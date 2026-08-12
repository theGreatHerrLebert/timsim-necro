"""timsim v2 as a necroflow DAG.

The three axes are a **cache model, not a cross product**, and that is exactly what a
content-addressed DAG is for:

    STRUCTURE   which molecules exist          computed ONCE, shared by every sample
    QUANTITY    how much of each, per sample   cheap; one node per sample
    MEASUREMENT how it is observed, per run    one node per run

necroflow's `DAG.add` deduplicates nodes by a content-addressed fingerprint, so the structure nodes
of N sample pipelines collapse to one set automatically. Nothing here asks for that; it falls out of
declaring the dependencies honestly. Adding a 20th sample re-runs `yield` and the simulator, and
re-runs nothing upstream of them.

That is the whole thesis of the redesign, expressed as a graph:

    proteome ─┬─ digest ─┬─ modify ── precursors ─┐
              │          │                        │
              │          └────────────────────────┼─── yield ── simulate(sample)  × N
              └─ design ──────────────────────────┘

Run:
    python timsim_flow.py --outdir /tmp/necro --dry-run
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
from types import SimpleNamespace

from necroflow import DAG, NodeType, Pipeline, command, output, resolve_command


# ── invalidation: a spec FILE is a dependency, not a string ──────────────────
#
# A node's fingerprint hashes its config *values*. A rule that takes `spec="design.toml"` therefore
# hashes the six characters of that filename — **not the mixture inside it.** Edit the fold change,
# re-run, and necroflow reports everything up to date and hands you the previous experiment's `.d`.
#
# That is exactly the failure this project exists to kill: a silently stale artifact that looks like
# a successful run. It is also invisible — nothing errors, the numbers are simply the old ones.
#
# So every artifact whose rule reads a spec file declares an `invalidator` that hashes the file's
# CONTENT. Change the file, change the token, and the node goes stale. Verified below, because a
# caching claim that has not been tested against an edit is not a caching claim.
def hashes_file(config_key: str):
    """Invalidate a node when the spec file named by `config_key` changes on disk."""

    def token(node) -> str:
        path = Path(node.config[config_key])
        return hashlib.sha256(path.read_bytes()).hexdigest()

    return token


def hashes_sciex_config(config_key: str):
    """Like `hashes_file`, but ALSO hashes the `.wiff` template the SCIEX config references.

    The legacy `timsim` CLI takes the template only *inside* the config TOML (there is no `--template`
    flag), so the template's CONTENT is a hidden dependency that hashing the TOML alone would miss —
    swapping the `.wiff` (or editing it in place) would silently reuse a stale render with the old SWATH
    windows / TOF calibration. `render_thermo` avoids this by hashing its template directly; do the same
    here by parsing the config for `template_path` and folding the file's bytes in. Fail graph
    construction if the referenced template is absent (a missing template is not a cache hit)."""
    import tomllib

    def token(node) -> str:
        cfg_path = Path(node.config[config_key])
        h = hashlib.sha256()
        h.update(cfg_path.read_bytes())
        data = tomllib.loads(cfg_path.read_text())
        tmpl = None
        for section in list(data.values()) + [data]:
            if isinstance(section, dict):
                tmpl = section.get("template_path") or section.get("astral_template_path")
                if tmpl:
                    break
        if not tmpl:
            raise ValueError(f"SCIEX config {cfg_path} declares no template_path to hash")
        tp = Path(tmpl)
        if not tp.exists():
            raise FileNotFoundError(f"SCIEX template {tp} (from {cfg_path}) does not exist")
        h.update(b"\x00wiff\x00")
        h.update(tp.read_bytes())
        return h.hexdigest()

    return token


def hashes_reference_d(config_key: str):
    """Invalidate a v2 Bruker render when its reference `.d` changes. The `.d` is a DIRECTORY (and often
    multi-GB), so — like `hashes_sciex_config` folds in the `.wiff` — hash the cheap `analysis.tdf`
    metadata DB inside it: that SQLite file carries the frame schedule, DIA windows and TOF/mobility
    calibration the render actually reads to place ions. Swap or edit the reference and the token moves;
    a missing reference is not a cache hit."""

    def token(node) -> str:
        ref = Path(node.config[config_key])
        tdf = ref / "analysis.tdf"
        if not tdf.exists():
            raise FileNotFoundError(f"reference .d {ref} has no analysis.tdf to hash")
        return hashlib.sha256(tdf.read_bytes()).hexdigest()

    return token


# ── artifacts ────────────────────────────────────────────────────────────────
#
# Each NodeType is a *typed artifact*, not a filename. The type is what lets necroflow check that a
# rule is fed the thing it asked for — the static half of the same discipline `timsim-schema`
# enforces at runtime. A stage that wants Peptides cannot be handed a Proteome, and the failure is at
# graph-construction time rather than three stages downstream.


class Proteome(NodeType):
    """The proteins. STRUCTURE — no amounts; abundance is a quantity and lives on its own axis."""

    filename = "proteome.parquet"
    invalidator = hashes_file("proteome_spec")


class Peptides(NodeType):
    """Distinct peptide sequences. STRUCTURE."""

    filename = "peptides.parquet"


class Occurrences(NodeType):
    """Where each peptide came from. The protein-inference answer key."""

    filename = "peptide_occurrences.parquet"


class CleavageSites(NodeType):
    """Every cleavage boundary. Needed to reconstruct the yield lattice."""

    filename = "cleavage_sites.parquet"


class Modforms(NodeType):
    """Modified species, with the fraction of molecules in each. STRUCTURE."""

    filename = "modforms.parquet"
    invalidator = hashes_file("mods_spec")


class Modifications(NodeType):
    """The modification spec as an artifact — read by BOTH modify and yield, so the two can
    never disagree about an occupancy."""

    filename = "modifications.parquet"
    invalidator = hashes_file("mods_spec")


class Precursors(NodeType):
    """The ion layer: m/z, isotope envelopes, charge fractions, ionisation propensity."""

    filename = "precursors.parquet"


class PrecursorCCS(NodeType):
    """Collision cross section per precursor. STRUCTURE — instrument-independent, predicted once.

    CCS is a property of the ion; 1/K0 is what an instrument measures from it. Keeping CCS here is
    what lets the same precursor space be measured on instrument A and instrument B (different gas)
    without recomputing anything upstream."""

    filename = "precursor_ccs.parquet"


class PeptideRT(NodeType):
    """Retention-time index per peptide. STRUCTURE — a hydrophobicity coordinate, gradient-independent,
    predicted once. Mapping it to seconds is a per-run measurement (the RT analog of CCS→1/K0)."""

    filename = "peptide_rt.parquet"


class Samples(NodeType):
    """The design: which samples exist, and the mixture that defines them."""

    filename = "samples.parquet"
    invalidator = hashes_file("design_spec")


class Runs(NodeType):
    """The runs, and which sample each measures."""

    filename = "runs.parquet"
    invalidator = hashes_file("design_spec")


class SampleRunMap(NodeType):
    """sample -> run. A sample measured twice is two runs, not two samples."""

    filename = "sample_run_map.parquet"
    invalidator = hashes_file("design_spec")


class ProteinQuantities(NodeType):
    """Protein amounts per sample, in amol. QUANTITY."""

    filename = "protein_quantities.parquet"
    invalidator = hashes_file("design_spec")


class PeptideQuantities(NodeType):
    """Peptide amounts per sample, in amol. QUANTITY — the digest applied to a mixture."""

    filename = "peptide_quantities.parquet"


class YieldReport(NodeType):
    """The `timsim-yield --report` accounting: missed-cleavage distribution, truncation/filter loss.

    A GUARANTEED co-output — yield always writes it when `--report` is passed, so it can be declared as
    one without ever poisoning cache classification by being intermittently absent. (Ported from the
    Sample Designer's pipeline, which wired it first.)"""

    filename = "yield_report.toml"


class RawData(NodeType):
    """A Bruker .d — the MEASUREMENT. One per run.

    A directory rather than a file: timsim names the `.d` itself from the experiment. The node's
    output is the directory that contains it.
    """

    filename = "raw"
    invalidator = hashes_file("config")


class FragmentPredictionInput(NodeType):
    """`(precursor_id, [UNIMOD]-annotated sequence, charge)` — the FROZEN input to the fragment model.
    Explicit (not hidden in the prediction rule) so the sequence + charge + MOD ENCODING are cached and
    inspectable: a modified precursor carries its modform's annotated sequence, so it fragments as
    modified. Built by the same `annotate()` the spectrum builder uses for m/z, so intensity and m/z
    agree on what the molecule is."""

    filename = "fragment_prediction_input.parquet"


class FragmentIntensities(NodeType):
    """Predicted fragment intensities — MEASUREMENT, and the instrument-DEPENDENT artifact. The
    fragment model (local timsTOF vs koina Orbitrap-HCD) is a node config value, so a different model
    is a different node; the materialised artifact is the auditable boundary for a network predictor."""

    filename = "fragment_intensities"


class IonSpectra(NodeType):
    """Instrument-independent MS1 isotopes + MS2 fragment spectra as `(m/z, intensity)`. One node bakes
    both (MS2 depends on the fragment model, MS1 does not — a later split can share MS1)."""

    filename = "ion_spectra"


class ThermoRawData(NodeType):
    """A Thermo `.raw` — the MEASUREMENT for a NO-IMS instrument (Orbitrap / Astral), authored into a
    real template. A DISTINCT type from Bruker [`RawData`] (`.d`) so a wrong-consumer is a type error.
    Restages when the template file changes."""

    filename = "data.raw"
    invalidator = hashes_file("template")


class ThermoTruth(NodeType):
    """The per-precursor answer key — a co-output of the render, so it is cached and invalidated EXACTLY
    with its `.raw` (never a drifting sidecar). A future search/score node consumes this by type."""

    filename = "truth.parquet"


class ThermoRunManifest(NodeType):
    """The auditable boundary for a render: renderer identity + version, template identity, fragment
    model, acquisition method, content-addressed input paths, and the render's own counts. Co-emitted
    with the `.raw` so a run is reproducible after the fact."""

    filename = "manifest.json"


class BrukerRawDataV2(NodeType):
    """A Bruker `.d` authored by the LEAN v2 projector (`timsim-render`) — the streaming, imspy-free render
    that places the same instrument-independent `ion_spectra` onto a reference `.d`'s acquisition grid. A
    DISTINCT type from the v1 [`RawData`] (`.d` via the monolithic `timsim` config seam): the v2 render is
    template-driven (a reference `.d`), so it restages on the reference's `analysis.tdf`, not a config file.
    NOTE: v2 DIA does not yet co-emit a truth answer key (only DDA does), so this path stops at the `.d` —
    phase-2 scoring stays on the Thermo path until `run_dia` writes a truth."""

    filename = "data.d"
    invalidator = hashes_reference_d("reference_d")


class BrukerTruthV2(NodeType):
    """The per-precursor DIA answer key co-emitted with the lean Bruker `.d` — the SAME 8-column schema as
    [`ThermoTruth`] (precursor_id, peptide_id, charge, mz, rt_seconds, abundance, has_ms2, in_any_window),
    so the instrument-agnostic `timsim_eval` score node consumes it unchanged. A co-output of the render,
    so it is cached and invalidated EXACTLY with its `.d` (never a drifting sidecar)."""

    filename = "truth.parquet"


class BrukerDdaData(NodeType):
    """A Bruker DDA-PASEF `.d` authored by `timsim-render --dda`: MS1 surveys + top-N precursor selection
    with dynamic exclusion + band-limited MS2. A DISTINCT type from the DIA [`BrukerRawDataV2`] so a
    wrong-consumer is a type error (DDA is searched by Sage, DIA by DiaNN). Restages on the reference `.d`."""

    filename = "data.d"
    invalidator = hashes_reference_d("reference_d")


class DdaTruth(NodeType):
    """The DDA answer key co-emitted with the `.d` (`--dda-truth`) — one row per SELECTION EVENT, tying
    each MS2 to the true precursor (ms2_frame, scan_begin/end, tdf_precursor_id, precursor_id, peptide_id,
    charge, mono_mz, rt_seconds, …). The `timsim_eval` DDA scorer maps Sage's `scannr` to `tdf_precursor_id`
    to score identifications against the precursors DDA actually fragmented. Cached with its `.d`."""

    filename = "dda_truth.parquet"


class SageReport(NodeType):
    """A Sage database search of a DDA `.d` — the SEARCH half of the DDA phase 2. A directory node (Sage
    writes `results.sage.tsv` + quant alongside). Restages when the search FASTA changes."""

    filename = "sage"


class SciexMzmlData(NodeType):
    """A SCIEX ZenoTOF SWATH run authored into open **mzML** by the LEAN v2 projector
    (`timsim-render-sciex`) — feature space → synthesised SWATH schedule → mzML via `timsim-core`'s
    mzdata writer, imspy-free and with NO `.wiff`/`sciexwiff` dependency. A single mzML FILE (not the v1
    directory). No config file: the SWATH params live in the command, so necroflow's fingerprint already
    covers them — no custom invalidator. Distinct from Thermo/Bruker so a wrong consumer is a type error."""

    filename = "sciex.mzML"


class SciexTruthV2(NodeType):
    """The per-precursor SWATH answer key co-emitted with the lean SCIEX mzML — the SAME 8-column schema
    as [`ThermoTruth`]/[`BrukerTruthV2`], so the instrument-agnostic `timsim_eval` scorer consumes it
    unchanged. Cached and invalidated with its mzML."""

    filename = "truth.parquet"


class DiannReport(NodeType):
    """A DiaNN library-free search of the rendered `.raw` — the SEARCH half of phase 2. A directory node
    (DiaNN emits report.parquet + stats + the predicted lib alongside). Restages when the search FASTA
    changes; the `.raw` it searches is an input node, so a different render is a different search."""

    filename = "diann"
    invalidator = hashes_file("search_fasta")


class ScoreMetrics(NodeType):
    """The SCORE half of phase 2: the DiaNN report scored against the render's answer key — hierarchical
    recall (all → present → in-window → has-frags → detectable), FDP, and a recall-vs-abundance-decile
    curve. This is the number the whole simulate→search→score run exists to produce, content-addressed
    to the exact `.raw` + truth + search DB that produced it."""

    filename = "metrics.json"


# ── rules ────────────────────────────────────────────────────────────────────


BIN = os.environ.get("TIMSIM_BIN", "target/release")
# Phase-2 search: DiaNN reads Thermo .raw natively only with the .NET 8 runtime (DOTNET_ROOT + on PATH).
DIANN = os.environ.get("TIMSIM_DIANN", "/home/administrator/dia-nn/diann-2.5.0/diann-linux")
DOTNET = os.environ.get("DOTNET_ROOT", os.path.expanduser("~/.dotnet"))
# Sage (DDA database search of Bruker `.d`, natively). Config supplies enzyme/mods/tolerances; the FASTA is
# overridden per run with `-f`. Built from lazear/sage with the local `.d` read patch.
SAGE = os.environ.get("TIMSIM_SAGE", "/home/administrator/Documents/promotion/rust/sage/target/release/sage")
SAGE_CONFIG = os.environ.get("TIMSIM_SAGE_CONFIG", "/scratch/timsim-demo/SAGEBench/configs/sage-smoke.json")


# Noise is passed to the render as INDIVIDUAL value placeholders (see the render nodes) — never as one
# multi-token string, which necroflow would shell-quote into a single broken arg. A1 signal-m/z (ppm/seed,
# 0 = off) is always passed; A2 / spike (presence flags) select a render VARIANT.
def _a1_kwargs(cfg) -> dict:
    return {
        "noise_mz_ppm": getattr(cfg, "noise_mz_ppm", 0.0),
        "noise_frag_ppm": getattr(cfg, "noise_frag_ppm", 0.0),
        "noise_seed": getattr(cfg, "noise_seed", 0),
    }


def _a2_kwargs(cfg) -> dict:
    return {
        "noise_precursor_frames": getattr(cfg, "noise_precursor_frames", 5),
        "noise_fragment_frames": getattr(cfg, "noise_fragment_frames", 5),
        "noise_intensity_max": getattr(cfg, "noise_intensity_max", 150000.0),
        "noise_precursor_fraction": getattr(cfg, "noise_precursor_fraction", 0.2),
        "noise_fragment_fraction": getattr(cfg, "noise_fragment_fraction", 0.2),
    }


def do_render(cfg, sample_id, P, control=False):
    """Pick + call the right Bruker render node for cfg's noise/spike mode. `control=True` renders the
    background-only variant (`--noise-only`) for the FDP-subtraction control. Returns (raw, truth)."""
    inputs = (P.precursors, P.rt, P.ion_spectra, P.ccs, P.peptide_quantities)
    common = dict(reference_d=cfg.reference_d, sample_id=sample_id,
                  intensity_scale=cfg.intensity_scale, **_a1_kwargs(cfg),
                  # Named explicitly rather than left to the function defaults: necroflow requires
                  # every declared input to be passed, and — more importantly — these are what the
                  # command-string fingerprint sees. See the note on _RENDER_HEAD.
                  peak_shape=cfg.peak_shape, cycle_seconds=cfg.cycle_seconds,
                  mobility_std_target=cfg.mobility_std_target, n_sigma=cfg.n_sigma,
                  ion_count_noise=str(cfg.ion_count_noise).lower(), instrument_cv=cfg.instrument_cv,
                  run_rt_sd=cfg.run_rt_sd, run_im_sd=cfg.run_im_sd,
                  run_intensity_cv=cfg.run_intensity_cv, target_p=cfg.target_p)
    if getattr(cfg, "spike_into", None):
        fn = render_spike_control if control else render_spike
        return fn(P, *inputs, spike_into=cfg.spike_into, **common)
    if getattr(cfg, "noise_real_data", False):
        fn = render_a2_control if control else render_a2
        return fn(P, *inputs, **_a2_kwargs(cfg), **common)
    return render(P, *inputs, **common)  # noiseless / A1 (no control)


@command(f"{BIN}/timsim-proteome --spec {{proteome_spec}} --out {{proteome}}")
def proteome(proteome_spec: str):
    """FASTAs -> proteins. STRUCTURE.

    Takes a multi-source spec rather than one FASTA, because that is what a real experiment is: HYE
    is three organisms, and the organism is a DECLARED column. v1 recovers it by substring-matching
    "HUMAN"/"YEAST"/"ECOLI" in the FASTA header, and peptides shared between two organisms silently
    become "Unknown" and get dropped.
    """
    proteome = output(Proteome)
    return proteome


@command(
    f"{BIN}/timsim-digest --proteome {{proteome}} "
    "--out-peptides {peptides} --out-occurrences {occurrences} --out-cleavage-sites {cleavage_sites} "
    "--max-missed-cleavages {max_missed_cleavages} --min-length {min_length} --max-length {max_length} "
    "--max-peptides {max_peptides} --seed {seed}"
)
def digest(proteome: Proteome, max_missed_cleavages: int, min_length: int, max_length: int,
           max_peptides: int, seed: int):
    """Proteins -> peptides. STRUCTURE, so it is computed once for every sample in the design.

    Three co-outputs of one call: they are one computation and necroflow treats them as such.
    `max_peptides=0` keeps the full analytic digest; a positive value samples that many (seeded) for a
    tractable run on a large proteome while keeping the full FASTA as the search space.
    """
    peptides = output(Peptides)
    occurrences = output(Occurrences)
    cleavage_sites = output(CleavageSites)
    return peptides, occurrences, cleavage_sites


@command(
    f"{BIN}/timsim-modify --peptides {{peptides}} --mods {{mods_spec}} "
    "--out-modforms {modforms} --out-modifications {modifications} --floor {floor}"
)
def modify(peptides: Peptides, mods_spec: str, floor: float):
    """Peptides -> modforms, driven by per-site OCCUPANCY rather than variable-mod combinatorics.

    Emits the modification spec alongside the modforms, because `yield` needs the same occupancies to
    know which cleavage sites are blocked. One artifact, two consumers, no flag to disagree about.
    """
    modforms = output(Modforms)
    modifications = output(Modifications)
    return modforms, modifications


@command(
    f"{BIN}/timsim-precursors --peptides {{peptides}} --modforms {{modforms}} "
    "--out {precursors} --charge-model {charge_model} --seed {seed} "
    # v1 drops charge states below a share of the charged mass (its dia-PASEF config uses 0.25);
    # v2 kept every state above zero, which renders precursors no instrument would record.
    # 0 = keep everything (v2's historical behaviour).
    "--min-charge-contrib {min_charge_contrib}"
)
def precursors(peptides: Peptides, modforms: Modforms, charge_model: str, seed: int,
               min_charge_contrib: float = 0.0):
    """Modforms -> ions. STRUCTURE: m/z and isotope envelopes are properties of the molecule, so
    they are shared by every sample too."""
    precursors = output(Precursors)
    return precursors


# ram: MEASURED, not guessed — 4.4 GB peak on a full-proteome HeLa run (9,010,877 precursors,
# 2026-08-01) after the chunked-tokenizer fix. Before that fix this stage peaked at 24.9 GB and was
# SIGKILLed. See the note above `fragments` on why these declarations matter.
@command(
    "timsim-ccs --precursors {precursors} --peptides {peptides} --out {precursor_ccs}",
    ram="6Gi",
)
def ccs(precursors: Precursors, peptides: Peptides):
    """Precursors -> CCS. STRUCTURE, and the one Python tool in the structure axis (the deep model
    is the standing exception). Runs once on the full precursor space and is shared by every sample
    and every simulated instrument."""
    precursor_ccs = output(PrecursorCCS)
    return precursor_ccs


@command("timsim-rt --peptides {peptides} --out {peptide_rt}")
def rt(peptides: Peptides):
    """Peptides -> RT index. STRUCTURE; deep model (Chronologer by default). Shared across every
    sample and every gradient."""
    peptide_rt = output(PeptideRT)
    return peptide_rt


@command(
    f"{BIN}/timsim-design --proteome {{proteome}} --spec {{design_spec}} "
    "--out-samples {samples} --out-runs {runs} --out-sample-run-map {sample_run_map} "
    "--out-protein-quantities {protein_quantities}"
)
def design(proteome: Proteome, design_spec: str):
    """The mixture. QUANTITY — this is where an A/B experiment is *declared* rather than recovered
    from a filename afterwards."""
    samples = output(Samples)
    runs = output(Runs)
    sample_run_map = output(SampleRunMap)
    protein_quantities = output(ProteinQuantities)
    return samples, runs, sample_run_map, protein_quantities


@command(
    f"{BIN}/timsim-yield --proteome {{proteome}} --occurrences {{occurrences}} "
    "--cleavage-sites {cleavage_sites} --protein-quantities {protein_quantities} "
    "--modifications {modifications} "
    "--digestion-efficiency {digestion_efficiency} --out {peptide_quantities} "
    "--report {yield_report}"
)
def peptide_yield(
    proteome: Proteome,
    occurrences: Occurrences,
    cleavage_sites: CleavageSites,
    protein_quantities: ProteinQuantities,
    modifications: Modifications,
    digestion_efficiency: float,
):
    """Structure x mixture -> peptide amounts. QUANTITY: cheap, and re-run whenever the design
    changes without touching the digest.

    Takes `modifications` so a blocking mod (acetyl-K, GG-K, TMT-K) actually stops the protease —
    the missed cleavage it forces is how the experiment localises the site.
    """
    peptide_quantities = output(PeptideQuantities)
    yield_report = output(YieldReport)
    return peptide_quantities, yield_report


@command(
    "mkdir -p {raw} && timsim {config} --save-path {raw} "
    "--v2-proteome {proteome} --v2-peptides {peptides} --v2-occurrences {occurrences} "
    "--v2-peptide-quantities {peptide_quantities} --v2-precursors {precursors} "
    "--v2-ccs {precursor_ccs} --v2-rt {peptide_rt} "
    "--v2-sample {sample_id} --seed {seed}",
    # Declared, because it was learned the hard way: the adapter loads the whole precursor table
    # into pandas, and a PTM-enriched design makes that table enormous — a phospho occupancy of 0.30
    # produced a **3.1 GB** precursors.parquet. Two of these ran in parallel and the kernel killed
    # one (exit 137). The blow-up is real chemistry (an enrichment genuinely has that many species),
    # so the answer is to tell the scheduler the truth about the cost rather than to hide it.
    threads=2,
    ram="8Gi",
)
def simulate(
    proteome: Proteome,
    peptides: Peptides,
    occurrences: Occurrences,
    peptide_quantities: PeptideQuantities,
    precursors: Precursors,
    precursor_ccs: PrecursorCCS,
    peptide_rt: PeptideRT,
    config: str,
    sample_id: str,
    seed: int,
):
    """MEASUREMENT: v1's LC, ion mobility, fragmentation and acquisition, driven from v2 artifacts.

    This is the strangler seam. One node per sample, and the only thing that changes between them is
    `sample_id` — so this fans out N ways while everything above it is computed once.
    """
    raw = output(RawData)
    return raw


# ── the measurement/render branch (Thermo, no-IMS) ───────────────────────────
# These replace the hand-fired steps: predict fragments (choosable model) -> assemble spectra ->
# author a real Thermo .raw. They reuse the SAME feature-space nodes, so the device/method matrix is a
# fan-out over (template, frag_model) with the feature space computed once.


@command(
    f"{BIN}/timsim-frag-input --precursors {{precursors}} --peptides {{peptides}} "
    "--modforms {modforms} --modifications {modifications} --out {fragment_prediction_input}"
)
def frag_input(precursors: Precursors, peptides: Peptides, modforms: Modforms, modifications: Modifications):
    """Precursors + modforms -> the frozen fragment-prediction input: `(precursor_id, [UNIMOD]-annotated
    sequence, charge)`. STRUCTURE (no CE/model), so it is shared by every fragment model. Annotates each
    precursor's modform, so a MODIFIED precursor fragments as modified — this is the correctness fix over
    the old bare-sequence join, which predicted every modform identically."""
    fragment_prediction_input = output(FragmentPredictionInput)
    return fragment_prediction_input


# ram: MEASURED at 15.4 GB peak on a full-proteome HeLa run (9,010,877 precursors -> 462,789,095
# fragment rows, 2026-08-01). This is the largest consumer in the pipeline and it was previously
# UNDECLARED, so necroflow scheduled it concurrently with `ccs` (both started in the same second on
# 2026-07-30). The combined footprint drove the box into deep swap, which silently corrupted 4
# floating-point values in the downstream `ion_spectra` artifact — see
# PhantomBENCH/simulations/v2/work/nodes/spectra/2a6187cf*/PROVENANCE-NOTE.md. Declaring real numbers
# is what stops the scheduler co-locating two stages that cannot both fit.
# NOTE: the historical note "streamed to ~1 GB" is wrong at this scale; treat 16Gi as the floor.
@command(
    "timsim-fragments --precursors {fragment_prediction_input} "
    "--collision-energy {collision_energy} --model {frag_model} --out {fragment_intensities}",
    ram="16Gi",
)
def fragments(fragment_prediction_input: FragmentPredictionInput, collision_energy: float, frag_model: str):
    """Annotated input -> predicted fragment intensities. MEASUREMENT, instrument-DEPENDENT: `frag_model`
    is "" (local timsTOF) or "koina:Prosit_2020_intensity_HCD" (Orbitrap-HCD), a config value — so the
    timsTOF-vs-HCD split is exactly two nodes and N renders sharing a model predict fragments once."""
    fragment_intensities = output(FragmentIntensities)
    return fragment_intensities


# ram: MEASURED at 4.86 GB peak on a full-proteome HeLa run (2026-08-01) with the chunked/streaming
# implementation (timsim-cli bb56ff2). The PREVIOUS implementation held the precursor->fragment map,
# every generated spectrum, and the Arrow copy simultaneously and peaked at 24.4 GB on the same input
# — do not restore a value sized for that version.
@command(
    f"{BIN}/timsim-spectra --precursors {{precursors}} --peptides {{peptides}} "
    "--modforms {modforms} --modifications {modifications} "
    "--fragment-intensities {fragment_intensities} --out {ion_spectra}",
    ram="6Gi",
)
def spectra(
    precursors: Precursors,
    peptides: Peptides,
    modforms: Modforms,
    modifications: Modifications,
    fragment_intensities: FragmentIntensities,
):
    """Precursors + fragments -> instrument-independent MS1 isotope + MS2 fragment spectra."""
    ion_spectra = output(IonSpectra)
    return ion_spectra


@command(
    f"{BIN}/timsim-render-thermo --precursors {{precursors}} "
    "--peptide-rt {peptide_rt} --ion-spectra {ion_spectra} --peptide-quantities {peptide_quantities} "
    "--sample {sample_id} --template {template} --intensity-scale {intensity_scale} "
    "--frag-model {frag_model} --method {method} --expected-ce {collision_energy} "
    "--out {data_raw} --thermo-truth {truth} --manifest {manifest}",
    threads=2,
    ram="8Gi",
)
def render_thermo(
    precursors: Precursors,
    peptide_rt: PeptideRT,
    ion_spectra: IonSpectra,
    peptide_quantities: PeptideQuantities,
    template: str,
    intensity_scale: float,
    sample_id: str,
    frag_model: str,
    method: str,
    collision_energy: float,
):
    """MEASUREMENT: author the feature space into a real Thermo `.raw` template (no-IMS). One node per
    sample (via `peptide_quantities` + `sample_id`); restages when the template changes. Three co-outputs
    of one command: the `.raw`, its answer key, and a durable run manifest — one computation, so the
    answer key and audit trail can never drift from the data."""
    data_raw = output(ThermoRawData)
    truth = output(ThermoTruth)
    manifest = output(ThermoRunManifest)
    return data_raw, truth, manifest


# ── the measurement/render branch (Bruker, WITH ion mobility — the lean v2 projector) ──
# The imspy-free counterpart of the monolithic v1 `simulate` seam: it reuses the SAME
# frag_input -> fragments -> spectra feature-space nodes as the Thermo path, then PROJECTS the
# instrument-independent spectra onto a reference `.d`'s DIA schedule. Bruker has ion mobility, so unlike
# the Thermo render it consumes `precursor_ccs` (CCS -> 1/K0, Mason-Schamp) — physical mobility a search
# engine needs. No config TOML: a reference `.d` (`--reference-d`) supplies the acquisition grid.


# The render command is assembled from single-token fragments. necroflow shell-quotes EACH `{placeholder}`
# (shlex.quote), so a multi-token string in one placeholder would be glued into a single arg — the flag NAMES
# must be literal in the template and only VALUES are placeholders. A1 signal-m/z is always passed (ppm 0 =
# off = byte-identical), so noiseless and A1 share one node; A2 / spike (presence flags) get node variants.
_RENDER_HEAD = (
    f"{BIN}/timsim-render --precursors {{precursors}} --peptide-rt {{peptide_rt}} "
    "--ion-spectra {ion_spectra} --precursor-ccs {precursor_ccs} "
    "--peptide-quantities {peptide_quantities} --sample {sample_id} "
    "--reference-d {reference_d} --dia --intensity-scale {intensity_scale} "
    # EXPLICIT even though these are the binary's own defaults. necroflow fingerprints on the COMMAND
    # STRING, so a default that changes behaviour without changing the command is invisible to the
    # cache: a render made by an older binary looks current and is silently reused. All four changed
    # on 2026-08-11 (peak shape went per-peptide, the clock started inheriting from the reference,
    # per-ion mobility widths arrived) and each moves the rendered SIGNAL. Naming them here is what
    # makes the cache correct rather than merely current -- the artifact provenance stamp detects a
    # stale result after the fact, it does not prevent one being produced and reused.
    "--peak-shape {peak_shape} --cycle-seconds {cycle_seconds} "
    "--mobility-std-target {mobility_std_target} --n-sigma {n_sigma} "
    # A3 + run-to-run measurement variation. A3 is counting statistics on the ANALYTE signal
    # (distinct from --noise-real-data, which is background); the run-* terms displace a peak
    # per run, keyed on (sample, precursor), so technical replicates stop being bit-identical.
    "--ion-count-noise {ion_count_noise} --instrument-cv {instrument_cv} "
    "--run-rt-sd {run_rt_sd} --run-im-sd {run_im_sd} "
    "--run-intensity-cv {run_intensity_cv} --target-p {target_p} "
    "--noise-mz-ppm {noise_mz_ppm} --noise-frag-ppm {noise_frag_ppm} --noise-seed {noise_seed}"
)
_RENDER_A2 = (
    " --noise-real-data --noise-precursor-frames {noise_precursor_frames} "
    "--noise-fragment-frames {noise_fragment_frames} --noise-intensity-max {noise_intensity_max} "
    "--noise-precursor-fraction {noise_precursor_fraction} --noise-fragment-fraction {noise_fragment_fraction}"
)
_RENDER_TAIL = " --out {raw} --truth {truth}"


@command(_RENDER_HEAD + _RENDER_TAIL, threads=2, ram="8Gi")
def render(
    precursors: Precursors, peptide_rt: PeptideRT, ion_spectra: IonSpectra, precursor_ccs: PrecursorCCS,
    peptide_quantities: PeptideQuantities, reference_d: str, sample_id: str, intensity_scale: float,
    noise_mz_ppm: float = 0.0, noise_frag_ppm: float = 0.0, noise_seed: int = 0,
    peak_shape: str = "per-peptide", cycle_seconds: float = 0.0,
    mobility_std_target: float = 0.009, n_sigma: float = 3.0,
    ion_count_noise: str = "false", instrument_cv: float = 0.0,
    run_rt_sd: float = 0.0, run_im_sd: float = 0.0,
    run_intensity_cv: float = 0.0, target_p: float = 0.0,
):
    """MEASUREMENT (Bruker): the lean v2 projector places `ion_spectra` onto the reference `.d`'s DIA grid.
    A1 signal-m/z noise is always wired (`--noise-mz-ppm/-frag-ppm`; 0 = off, byte-identical). One node per
    sample; co-emits the per-precursor answer key (`--truth`)."""
    raw = output(BrukerRawDataV2)
    truth = output(BrukerTruthV2)
    return raw, truth


@command(_RENDER_HEAD + _RENDER_A2 + _RENDER_TAIL, threads=2, ram="8Gi")
def render_a2(
    precursors: Precursors, peptide_rt: PeptideRT, ion_spectra: IonSpectra, precursor_ccs: PrecursorCCS,
    peptide_quantities: PeptideQuantities, reference_d: str, sample_id: str, intensity_scale: float,
    noise_mz_ppm: float, noise_frag_ppm: float, noise_seed: int,
    noise_precursor_frames: int, noise_fragment_frames: int, noise_intensity_max: float,
    noise_precursor_fraction: float, noise_fragment_fraction: float,
    peak_shape: str = "per-peptide", cycle_seconds: float = 0.0,
    mobility_std_target: float = 0.009, n_sigma: float = 3.0,
    ion_count_noise: str = "false", instrument_cv: float = 0.0,
    run_rt_sd: float = 0.0, run_im_sd: float = 0.0,
    run_intensity_cv: float = 0.0, target_p: float = 0.0,
):
    """render + A2 real-data background sampled from the reference `.d` (the v1 DIA recipe with A1)."""
    raw = output(BrukerRawDataV2)
    truth = output(BrukerTruthV2)
    return raw, truth


@command(_RENDER_HEAD + _RENDER_A2 + " --noise-only" + _RENDER_TAIL, threads=2, ram="8Gi")
def render_a2_control(
    precursors: Precursors, peptide_rt: PeptideRT, ion_spectra: IonSpectra, precursor_ccs: PrecursorCCS,
    peptide_quantities: PeptideQuantities, reference_d: str, sample_id: str, intensity_scale: float,
    noise_mz_ppm: float, noise_frag_ppm: float, noise_seed: int,
    noise_precursor_frames: int, noise_fragment_frames: int, noise_intensity_max: float,
    noise_precursor_fraction: float, noise_fragment_fraction: float,
    peak_shape: str = "per-peptide", cycle_seconds: float = 0.0,
    mobility_std_target: float = 0.009, n_sigma: float = 3.0,
    ion_count_noise: str = "false", instrument_cv: float = 0.0,
    run_rt_sd: float = 0.0, run_im_sd: float = 0.0,
    run_intensity_cv: float = 0.0, target_p: float = 0.0,
):
    """A2 background-ONLY control (`--noise-only`): the real-data background alone, same seed — searched, its
    IDs subtracted from FDP (score_bruker_bg)."""
    raw = output(BrukerRawDataV2)
    truth = output(BrukerTruthV2)
    return raw, truth


@command(_RENDER_HEAD + " --spike-into {spike_into}" + _RENDER_TAIL, threads=2, ram="8Gi")
def render_spike(
    precursors: Precursors, peptide_rt: PeptideRT, ion_spectra: IonSpectra, precursor_ccs: PrecursorCCS,
    peptide_quantities: PeptideQuantities, reference_d: str, sample_id: str, intensity_scale: float,
    noise_mz_ppm: float, noise_frag_ppm: float, noise_seed: int, spike_into: str,
    peak_shape: str = "per-peptide", cycle_seconds: float = 0.0,
    mobility_std_target: float = 0.009, n_sigma: float = 3.0,
    ion_count_noise: str = "false", instrument_cv: float = 0.0,
    run_rt_sd: float = 0.0, run_im_sd: float = 0.0,
    run_intensity_cv: float = 0.0, target_p: float = 0.0,
):
    """Spike-into-real: overlay the synthetic signal additively onto a real `.d` (`--spike-into`)."""
    raw = output(BrukerRawDataV2)
    truth = output(BrukerTruthV2)
    return raw, truth


@command(_RENDER_HEAD + " --spike-into {spike_into} --noise-only" + _RENDER_TAIL, threads=2, ram="8Gi")
def render_spike_control(
    precursors: Precursors, peptide_rt: PeptideRT, ion_spectra: IonSpectra, precursor_ccs: PrecursorCCS,
    peptide_quantities: PeptideQuantities, reference_d: str, sample_id: str, intensity_scale: float,
    noise_mz_ppm: float, noise_frag_ppm: float, noise_seed: int, spike_into: str,
    peak_shape: str = "per-peptide", cycle_seconds: float = 0.0,
    mobility_std_target: float = 0.009, n_sigma: float = 3.0,
    ion_count_noise: str = "false", instrument_cv: float = 0.0,
    run_rt_sd: float = 0.0, run_im_sd: float = 0.0,
    run_intensity_cv: float = 0.0, target_p: float = 0.0,
):
    """Spike background control (`--spike-into X --noise-only`): a re-encoded copy of X, no synthetic —
    searched, its IDs subtracted from FDP."""
    raw = output(BrukerRawDataV2)
    truth = output(BrukerTruthV2)
    return raw, truth


# ── Bruker DDA-PASEF → `.d` (top-N selection + dynamic exclusion; searched by Sage, not DiaNN) ──
# The DDA counterpart of the DIA render: same feature-space chain (frag_input → fragments → spectra) +
# CCS, then `timsim-render --dda` synthesises MS1 surveys + top-N precursor selection with dynamic
# exclusion + band-limited MS2, and co-emits a per-SELECTION-EVENT answer key.


@command(
    f"{BIN}/timsim-render --dda --precursors {{precursors}} --peptide-rt {{peptide_rt}} "
    "--ion-spectra {ion_spectra} --precursor-ccs {precursor_ccs} --peptide-quantities {peptide_quantities} "
    "--sample {sample_id} --reference-d {reference_d} --intensity-scale {intensity_scale} "
    # Same fingerprint reasoning as _RENDER_HEAD: these are the binary's defaults, named here
    # so a behaviour change in them invalidates the cache instead of hiding in it. The DDA path
    # took the same peak-shape / clock / per-ion-mobility changes as DIA.
    "--peak-shape {peak_shape} --cycle-seconds {cycle_seconds} "
    "--mobility-std-target {mobility_std_target} --n-sigma {n_sigma} "
    "--precursors-every {precursors_every} --max-precursors {max_precursors} --exclusion-width {exclusion_width} "
    "--out {raw} --dda-truth {truth}",
    threads=2,
    ram="8Gi",
)
def render_dda(
    precursors: Precursors,
    peptide_rt: PeptideRT,
    ion_spectra: IonSpectra,
    precursor_ccs: PrecursorCCS,
    peptide_quantities: PeptideQuantities,
    reference_d: str,
    sample_id: str,
    intensity_scale: float,
    precursors_every: int,
    max_precursors: int,
    exclusion_width: int,
    peak_shape: str = "per-peptide", cycle_seconds: float = 0.0,
    mobility_std_target: float = 0.009, n_sigma: float = 3.0,
):
    """MEASUREMENT (Bruker DDA-PASEF): MS1 surveys every `precursors_every` frames, top-N (`max_precursors`)
    precursor selection with `exclusion_width`-frame dynamic exclusion, band-limited MS2 on the reference's
    scan geometry. Co-emits `--dda-truth` (one row per selection event) so a Sage search of the `.d` scores
    against exactly the precursors DDA fragmented. Restages when the reference `.d` changes."""
    raw = output(BrukerDdaData)
    truth = output(DdaTruth)
    return raw, truth


# ── SCIEX ZenoTOF SWATH → open mzML (no-IMS, LEAN v2 projector, synthesised schedule) ──
# The imspy-free counterpart of the v1 `timsim` build-from-`.wiff`: it reuses the SAME
# frag_input → fragments → spectra feature-space nodes as the Thermo/Bruker paths, then projects the
# instrument-independent spectra onto a SYNTHESISED SWATH schedule and writes open mzML — no `.wiff`, no
# `sciexwiff`/`sciex-io` (legally clean; native `.wiff` is a separate rustims-local satellite).


@command(
    f"{BIN}/timsim-render-sciex --precursors {{precursors}} --peptide-rt {{peptide_rt}} "
    "--ion-spectra {ion_spectra} --peptide-quantities {peptide_quantities} --sample {sample_id} "
    "--gradient-length-s {gradient_length_s} --cycle-time-s {cycle_time_s} "
    "--mz-min {mz_min} --mz-max {mz_max} --window-width {window_width} "
    "--collision-energy {collision_energy} --intensity-scale {intensity_scale} "
    "--frag-model {frag_model} --out {mzml} --truth {truth}",
    threads=2,
    ram="8Gi",
)
def render_sciex(
    precursors: Precursors,
    peptide_rt: PeptideRT,
    ion_spectra: IonSpectra,
    peptide_quantities: PeptideQuantities,
    sample_id: str,
    gradient_length_s: float,
    cycle_time_s: float,
    mz_min: float,
    mz_max: float,
    window_width: float,
    collision_energy: float,
    intensity_scale: float,
    frag_model: str,
):
    """MEASUREMENT (SCIEX ZenoTOF SWATH, no-IMS): the lean v2 projector places the instrument-independent
    `ion_spectra` onto a synthesised fixed-width SWATH schedule and writes open **mzML** — imspy-free, no
    `.wiff`. One node per sample (via `peptide_quantities` + `sample_id`); the SWATH params are in the
    command so the fingerprint covers them (no template file). Co-emits the per-precursor answer key
    (`--truth`) so a DiaNN search of the mzML closes search→score like the Thermo/Bruker paths."""
    mzml = output(SciexMzmlData)
    truth = output(SciexTruthV2)
    return mzml, truth


# ── phase 2: search + score (close simulate -> search -> score) ──────────────


@command(
    f"mkdir -p {{diann}} && DOTNET_ROOT={DOTNET} PATH={DOTNET}:$PATH {DIANN} "
    "--f {data_raw} --fasta {search_fasta} --out {diann}/report.parquet "
    "--fasta-search --predictor --gen-spec-lib --qvalue {qvalue} --threads {search_threads} "
    "--met-excision --cut 'K*,R*' --missed-cleavages {max_missed_cleavages} "
    "--min-pep-len {min_length} --max-pep-len {max_length} --var-mods 1 --unimod35 "
    "--reanalyse --relaxed-prot-inf",
    threads=16,
    ram="32Gi",
)
def search(
    data_raw: ThermoRawData,
    search_fasta: str,
    qvalue: float,
    search_threads: int,
    max_missed_cleavages: int,
    min_length: int,
    max_length: int,
):
    """SEARCH: DiaNN library-free over the rendered `.raw` (predict a spectral library from the FASTA,
    then search). Reads `.raw` natively via the .NET 8 runtime. The FASTA is a content-hashed dependency;
    a different render (`.raw` input) or a different DB is a different search."""
    diann = output(DiannReport)
    return diann


@command(
    "python -m timsim_eval.v2_thermo_eval "
    "--report {diann}/report.parquet --truth {truth} --peptides {peptides} "
    "--fdr {qvalue} --out {metrics}"
)
def score(diann: DiannReport, truth: ThermoTruth, peptides: Peptides, qvalue: float):
    """SCORE: the DiaNN report against the render's answer key. Hierarchical recall + FDP + recall-by-
    abundance-decile, content-addressed to the exact `.raw`/truth/DB that produced it — so the number
    can never be attributed to the wrong run."""
    metrics = output(ScoreMetrics)
    return metrics


# ── phase 2 for the lean Bruker `.d` (dia-PASEF; DiaNN reads `.d` NATIVELY — no .NET) ──


@command(
    f"mkdir -p {{diann}} && {DIANN} "
    "--f {data_d} --fasta {search_fasta} --out {diann}/report.parquet "
    "--fasta-search --predictor --gen-spec-lib --qvalue {qvalue} --threads {search_threads} "
    "--met-excision --cut 'K*,R*' --missed-cleavages {max_missed_cleavages} "
    "--min-pep-len {min_length} --max-pep-len {max_length} --var-mods 1 --unimod35",
    threads=16,
    ram="32Gi",
)
def search_bruker(
    data_d: BrukerRawDataV2,
    search_fasta: str,
    qvalue: float,
    search_threads: int,
    max_missed_cleavages: int,
    min_length: int,
    max_length: int,
):
    """SEARCH (Bruker dia-PASEF): DiaNN library-free over the rendered `.d`. Unlike the Thermo `.raw`, a
    Bruker `.d` is DiaNN's NATIVE input on Linux — no .NET runtime. The FASTA is a content-hashed
    dependency; a different render or DB is a different search."""
    diann = output(DiannReport)
    return diann



# Same as `search_bruker`, but reusing an already-built predicted library instead of
# regenerating it. The library is a deterministic function of the FASTA and the digest
# settings, NOT of the data -- for the 20,535-protein human proteome it costs ~42 min and
# 2.1 GB, so a 15-arm ramp that rebuilds it per arm burns ~11 h producing identical bytes.
@command(
    f"mkdir -p {{diann}} && {DIANN} "
    "--f {data_d} --fasta {search_fasta} --out {diann}/report.parquet "
    "--lib {speclib} --qvalue {qvalue} --threads {search_threads} "
    "--met-excision --cut 'K*,R*' --missed-cleavages {max_missed_cleavages} "
    "--min-pep-len {min_length} --max-pep-len {max_length} --var-mods 1 --unimod35",
    threads=16,
    ram="32Gi",
)
def search_bruker_lib(
    data_d: BrukerRawDataV2,
    search_fasta: str,
    speclib: str,
    qvalue: float,
    search_threads: int,
    max_missed_cleavages: int,
    min_length: int,
    max_length: int,
):
    """SEARCH (Bruker dia-PASEF): DiaNN library-free over the rendered `.d`. Unlike the Thermo `.raw`, a
    Bruker `.d` is DiaNN's NATIVE input on Linux — no .NET runtime. The FASTA is a content-hashed
    dependency; a different render or DB is a different search."""
    diann = output(DiannReport)
    return diann


@command(
    "python -m timsim_eval.v2_thermo_eval "
    "--report {diann}/report.parquet --truth {truth} --peptides {peptides} "
    "--fdr {qvalue} --out {metrics}"
)
def score_bruker(diann: DiannReport, truth: BrukerTruthV2, peptides: Peptides, qvalue: float):
    """SCORE (Bruker): identical to `score`, but keyed on the Bruker answer key [`BrukerTruthV2`]. The
    `timsim_eval` scorer is instrument-agnostic — the truth schema is the same 8 columns — so the Bruker
    `.d` closes search→score exactly like the Thermo path."""
    metrics = output(ScoreMetrics)
    return metrics


@command(
    "python -m timsim_eval.v2_thermo_eval "
    "--report {diann}/report.parquet --truth {truth} --peptides {peptides} "
    "--background-report {background}/report.parquet "
    "--fdr {qvalue} --out {metrics}"
)
def score_bruker_bg(
    diann: DiannReport,
    background: DiannReport,
    truth: BrukerTruthV2,
    peptides: Peptides,
    qvalue: float,
):
    """SCORE (Bruker, A2 real-data noise): like `score_bruker`, but subtracts the noise-only control's IDs
    from FDP. With A2 the reference blank's real peptides get identified by the search — they are NOT false
    positives against the synthetic answer key, so counting them would inflate FDP. `background` is the
    DiaNN report of the `--noise-only` control render (same seed, background alone); the scorer removes its
    IDs (`background − ground_truth`, never a real synthetic hit). See REALISM_PLAN.md."""
    metrics = output(ScoreMetrics)
    return metrics


# ── HYE quant (P0.2): a JOINT DiaNN run over both conditions' `.d`, then per-organism fold-change ──


@command(
    f"mkdir -p {{diann}} && {DIANN} "
    "--f {data_a} --f {data_b} --fasta {search_fasta} --out {diann}/report.parquet "
    "--fasta-search --predictor --gen-spec-lib --qvalue {qvalue} --threads {search_threads} "
    "--met-excision --cut 'K*,R*' --missed-cleavages {max_missed_cleavages} "
    "--min-pep-len {min_length} --max-pep-len {max_length} --var-mods 1 --unimod35",
    threads=16,
    ram="32Gi",
)
def search_bruker_joint(
    data_a: BrukerRawDataV2,
    data_b: BrukerRawDataV2,
    search_fasta: str,
    qvalue: float,
    search_threads: int,
    max_missed_cleavages: int,
    min_length: int,
    max_length: int,
):
    """SEARCH (Bruker, JOINT): the two HYE conditions' `.d` in ONE DiaNN run (`--f A --f B`), so DiaNN does
    cross-run normalization + MaxLFQ as a single experiment — the quant benchmark needs joint quant, which
    a merge of two independent searches cannot recover. `data_a` is condition A (reference, Run.Index 0),
    `data_b` is B (Run.Index 1). See HYE_QUANT.md."""
    diann = output(DiannReport)
    return diann


@command(
    "python -m timsim_eval.v2_quant_eval "
    "--report {diann}/report.parquet --peptides {peptides} --occurrences {occurrences} "
    "--proteome {proteome} --design {design_spec} --run-col Run.Index --run-a 0 --run-b 1 "
    "--fdr {qvalue} --delta {delta} --out {metrics}",
)
def quant(
    diann: DiannReport,
    peptides: Peptides,
    occurrences: Occurrences,
    proteome: Proteome,
    design_spec: str,
    qvalue: float,
    delta: float,
):
    """SCORE (HYE quant): per-organism log2 fold-change recovery from the joint report. Organism-unique
    peptide-level `log2(B/A)` vs the design's expected ratios; median residual + MAD + %correct, three
    normalization views, eligibility + detection tables. `design` (design.toml) supplies the expected
    ratios (content-hashed so editing the mixture restages). See HYE_QUANT.md."""
    metrics = output(ScoreMetrics)
    return metrics


# ── phospho site-localization (P1.3): DiaNN with --monitor-mod, then FLR ──


@command(
    f"mkdir -p {{diann}} && {DIANN} "
    "--f {data_d} --fasta {search_fasta} --out {diann}/report.parquet "
    "--fasta-search --predictor --gen-spec-lib --qvalue {qvalue} --threads {search_threads} "
    "--met-excision --cut 'K*,R*' --missed-cleavages {max_missed_cleavages} "
    "--min-pep-len {min_length} --max-pep-len {max_length} "
    "--var-mod 'UniMod:21,79.966331,STY' --monitor-mod 'UniMod:21' --var-mods 2",
    threads=16,
    ram="32Gi",
)
def search_bruker_phospho(
    data_d: BrukerRawDataV2,
    search_fasta: str,
    qvalue: float,
    search_threads: int,
    max_missed_cleavages: int,
    min_length: int,
    max_length: int,
):
    """SEARCH (Bruker, phospho): DiaNN with variable Phospho on STY and `--monitor-mod` — which enables
    per-site LOCALIZATION (PTM.Site.Confidence + Site.Occupancy.Probabilities in the report). See
    PHOSPHO_FLR.md."""
    diann = output(DiannReport)
    return diann


@command(
    "python -m timsim_eval.v2_flr_eval "
    "--report {diann}/report.parquet --truth {truth} --precursors {precursors} "
    "--modforms {modforms} --peptides {peptides} --fdr {qvalue} --target-flr {flr_target} --out {metrics}",
)
def score_flr(
    diann: DiannReport,
    truth: BrukerTruthV2,
    precursors: Precursors,
    modforms: Modforms,
    peptides: Peptides,
    qvalue: float,
    flr_target: float,
):
    """SCORE (phospho FLR): empirical false-localization rate vs simulator truth. The true site comes from
    the render's modforms (`precursor→modform→mod_positions`); DiaNN's localized site + `PTM.Site.Confidence`
    give the FLR(τ) / recall(τ) curve over the isolated single-isomer eligible set. See PHOSPHO_FLR.md."""
    metrics = output(ScoreMetrics)
    return metrics


# ── phase 2 for the lean SCIEX mzML (DiaNN reads open mzML NATIVELY — no .NET) ──


@command(
    f"mkdir -p {{diann}} && {DIANN} "
    "--f {mzml} --fasta {search_fasta} --out {diann}/report.parquet "
    "--fasta-search --predictor --gen-spec-lib --qvalue {qvalue} --threads {search_threads} "
    "--met-excision --cut 'K*,R*' --missed-cleavages {max_missed_cleavages} "
    "--min-pep-len {min_length} --max-pep-len {max_length} --var-mods 1 --unimod35",
    threads=16,
    ram="32Gi",
)
def search_sciex(
    mzml: SciexMzmlData,
    search_fasta: str,
    qvalue: float,
    search_threads: int,
    max_missed_cleavages: int,
    min_length: int,
    max_length: int,
):
    """SEARCH (SCIEX SWATH): DiaNN library-free over the rendered open **mzML** — DiaNN's native open
    input, no .NET, no vendor SDK. The FASTA is a content-hashed dependency."""
    diann = output(DiannReport)
    return diann


@command(
    "python -m timsim_eval.v2_thermo_eval "
    "--report {diann}/report.parquet --truth {truth} --peptides {peptides} "
    "--fdr {qvalue} --out {metrics}"
)
def score_sciex(diann: DiannReport, truth: SciexTruthV2, peptides: Peptides, qvalue: float):
    """SCORE (SCIEX): the same instrument-agnostic `timsim_eval` scorer, keyed on the SWATH answer key
    [`SciexTruthV2`] — the SCIEX mzML closes search→score exactly like Thermo/Bruker."""
    metrics = output(ScoreMetrics)
    return metrics


# ── phase 2 for Bruker DDA-PASEF (Sage database search of the `.d`; NOT DiaNN, which is DIA-only) ──


@command(
    f"mkdir -p {{sage}} && {SAGE} -f {{search_fasta}} -o {{sage}} {SAGE_CONFIG} {{data_d}}",
    threads=16,
    ram="32Gi",
)
def search_dda(data_d: BrukerDdaData, search_fasta: str):
    """SEARCH (Bruker DDA): Sage database search over the rendered `.d` (native `.d` reader). The Sage
    config (enzyme/mods/tolerances) is `TIMSIM_SAGE_CONFIG`; the FASTA is overridden per run with `-f`, so
    the same content-hashed dependency logic applies (a different render or DB is a different search)."""
    sage = output(SageReport)
    return sage


@command(
    "python -m timsim_eval.v2_dda_eval "
    "--sage {sage}/results.sage.tsv --truth {truth} --peptides {peptides} --fdr {qvalue} --out {metrics}"
)
def score_dda(sage: SageReport, truth: DdaTruth, peptides: Peptides, qvalue: float):
    """SCORE (Bruker DDA): map Sage's PSMs (`scannr` → `tdf_precursor_id`) onto the selection-event answer
    key. DDA recall is CONDITIONAL — over the precursors DDA actually fragmented (top-N), not all
    precursors — because DDA only selects a subset per cycle."""
    metrics = output(ScoreMetrics)
    return metrics


# ── the pipeline ─────────────────────────────────────────────────────────────


def timsim_pipeline(P: Pipeline, cfg, sample_id: str) -> None:
    """One sample, end to end.

    Every node above `peptide_yield` depends only on things that do NOT vary with the sample, so N
    calls to this function produce N pipelines whose structure nodes share a fingerprint — and
    necroflow collapses them. The cache model is not implemented here; it is *implied* by declaring
    the dependencies honestly.
    """
    P.proteome = proteome(P, proteome_spec=cfg.proteome_spec)

    P.peptides, P.occurrences, P.cleavage_sites = digest(P, 
        P.proteome,
        max_missed_cleavages=cfg.max_missed_cleavages,
        min_length=cfg.min_length,
        max_length=cfg.max_length,
        max_peptides=cfg.max_peptides,
        seed=cfg.seed,
    )

    P.modforms, P.modifications = modify(P, P.peptides, mods_spec=cfg.mods, floor=cfg.floor)

    P.precursors = precursors(P, 
        P.peptides, P.modforms, charge_model=cfg.charge_model, seed=cfg.seed
    )
    P.ccs = ccs(P, P.precursors, P.peptides)
    P.rt = rt(P, P.peptides)

    # The seed lives INSIDE the design spec, not on the command line — it is a property of the
    # experiment, and a flag would let the artifact and the caller disagree about it.
    P.samples, P.runs, P.sample_run_map, P.protein_quantities = design(P, 
        P.proteome, design_spec=cfg.design_spec
    )

    P.peptide_quantities, P.yield_report = peptide_yield(P, 
        P.proteome,
        P.occurrences,
        P.cleavage_sites,
        P.protein_quantities,
        P.modifications,
        digestion_efficiency=cfg.digestion_efficiency,
    )

    # The fan-out. `sample_id` is the ONLY thing that differs between samples, so this node's
    # fingerprint differs and everything upstream of it does not.
    P.raw = simulate(P, 
        P.proteome,
        P.peptides,
        P.occurrences,
        P.peptide_quantities,
        P.precursors,
        P.ccs,
        P.rt,
        config=cfg.timsim_config,
        sample_id=sample_id,
        seed=cfg.seed,
    )


def timsim_thermo_pipeline(P: Pipeline, cfg, sample_id: str) -> None:
    """One sample, end to end, but the MEASUREMENT is a Thermo `.raw` authored into a template (no-IMS
    Orbitrap / Astral) via explicit `fragments -> spectra -> render_thermo` nodes.

    The feature-space nodes are IDENTICAL to `timsim_pipeline`'s (same config values -> same
    fingerprint), so necroflow collapses them: request both a Bruker and a Thermo pipeline and the whole
    structure axis is computed once. CCS is omitted (no ion mobility). The instrument/method matrix is a
    fan-out over `cfg.template` × `cfg.frag_model`.
    """
    P.proteome = proteome(P, proteome_spec=cfg.proteome_spec)
    P.peptides, P.occurrences, P.cleavage_sites = digest(P, 
        P.proteome,
        max_missed_cleavages=cfg.max_missed_cleavages,
        min_length=cfg.min_length,
        max_length=cfg.max_length,
        max_peptides=cfg.max_peptides,
        seed=cfg.seed,
    )
    P.modforms, P.modifications = modify(P, P.peptides, mods_spec=cfg.mods, floor=cfg.floor)
    P.precursors = precursors(P, P.peptides, P.modforms, charge_model=cfg.charge_model, seed=cfg.seed,
                             min_charge_contrib=cfg.min_charge_contrib)
    P.rt = rt(P, P.peptides)
    P.samples, P.runs, P.sample_run_map, P.protein_quantities = design(P, 
        P.proteome, design_spec=cfg.design_spec
    )
    P.peptide_quantities, P.yield_report = peptide_yield(P, 
        P.proteome,
        P.occurrences,
        P.cleavage_sites,
        P.protein_quantities,
        P.modifications,
        digestion_efficiency=cfg.digestion_efficiency,
    )
    # ── measurement branch (was hand-fired) ──
    P.fragment_prediction_input = frag_input(P, 
        P.precursors, P.peptides, P.modforms, P.modifications
    )
    P.fragment_intensities = fragments(P, 
        P.fragment_prediction_input, collision_energy=cfg.collision_energy, frag_model=cfg.frag_model
    )
    P.ion_spectra = spectra(P, 
        P.precursors, P.peptides, P.modforms, P.modifications, P.fragment_intensities
    )
    P.raw, P.truth, P.manifest = render_thermo(P, 
        P.precursors,
        P.rt,
        P.ion_spectra,
        P.peptide_quantities,
        template=cfg.template,
        intensity_scale=cfg.intensity_scale,
        sample_id=sample_id,
        frag_model=cfg.frag_model,
        method=getattr(cfg, "method", "DIA"),
        collision_energy=cfg.collision_energy,
    )
    # ── phase 2 (opt-in): search the .raw + score against the answer key ──
    if getattr(cfg, "search_fasta", None):
        P.diann = search(P, 
            P.raw,
            search_fasta=cfg.search_fasta,
            qvalue=cfg.qvalue,
            search_threads=cfg.search_threads,
            max_missed_cleavages=cfg.max_missed_cleavages,
            min_length=cfg.min_length,
            max_length=cfg.max_length,
        )
        P.score = score(P, P.diann, P.truth, P.peptides, qvalue=cfg.qvalue)


def timsim_bruker_v2_pipeline(P: Pipeline, cfg, sample_id: str) -> None:
    """One sample to a Bruker `.d`, but via the LEAN v2 projector (`timsim-render`) instead of the
    monolithic v1 `simulate` seam — imspy-free, same `frag_input -> fragments -> spectra` chain as the
    Thermo path plus CCS (Bruker has ion mobility). The feature-space nodes are IDENTICAL to the other
    pipelines (same config -> same fingerprint), so requesting a Bruker-v2, a Thermo and a SCIEX run
    collapses the whole structure axis to one computation. Opt-in via `--bruker-reference <ref.d>`; the
    v1 `timsim_pipeline` stays the default (it owns DDA + the DIA truth output v2 does not yet emit)."""
    P.proteome = proteome(P, proteome_spec=cfg.proteome_spec)
    P.peptides, P.occurrences, P.cleavage_sites = digest(P, 
        P.proteome,
        max_missed_cleavages=cfg.max_missed_cleavages,
        min_length=cfg.min_length,
        max_length=cfg.max_length,
        max_peptides=cfg.max_peptides,
        seed=cfg.seed,
    )
    P.modforms, P.modifications = modify(P, P.peptides, mods_spec=cfg.mods, floor=cfg.floor)
    P.precursors = precursors(P, P.peptides, P.modforms, charge_model=cfg.charge_model, seed=cfg.seed,
                             min_charge_contrib=cfg.min_charge_contrib)
    P.ccs = ccs(P, P.precursors, P.peptides)
    P.rt = rt(P, P.peptides)
    P.samples, P.runs, P.sample_run_map, P.protein_quantities = design(P, 
        P.proteome, design_spec=cfg.design_spec
    )
    P.peptide_quantities, P.yield_report = peptide_yield(P, 
        P.proteome,
        P.occurrences,
        P.cleavage_sites,
        P.protein_quantities,
        P.modifications,
        digestion_efficiency=cfg.digestion_efficiency,
    )
    # ── measurement branch: same feature-space nodes as Thermo, then the v2 Bruker projector ──
    P.fragment_prediction_input = frag_input(P, 
        P.precursors, P.peptides, P.modforms, P.modifications
    )
    P.fragment_intensities = fragments(P, 
        P.fragment_prediction_input, collision_energy=cfg.collision_energy, frag_model=cfg.frag_model
    )
    P.ion_spectra = spectra(P, 
        P.precursors, P.peptides, P.modforms, P.modifications, P.fragment_intensities
    )
    P.raw, P.truth = do_render(cfg, sample_id, P)
    # ── phase 2 (opt-in): DiaNN-search the .d natively + score against the answer key ──
    if getattr(cfg, "search_fasta", None) and getattr(cfg, "phospho", False):
        # Phospho site-localization: DiaNN with --monitor-mod (localization), scored by FLR vs the render's
        # true modform sites — NOT the FDP path. Use with --mods mods_phospho.toml. See PHOSPHO_FLR.md.
        P.diann = search_bruker_phospho(P, 
            P.raw,
            search_fasta=cfg.search_fasta,
            qvalue=cfg.qvalue,
            search_threads=cfg.search_threads,
            max_missed_cleavages=cfg.max_missed_cleavages,
            min_length=cfg.min_length,
            max_length=cfg.max_length,
        )
        P.score = score_flr(P, 
            P.diann, P.truth, P.precursors, P.modforms, P.peptides,
            qvalue=cfg.qvalue, flr_target=cfg.flr_target,
        )
    elif getattr(cfg, "search_fasta", None):
        # A prebuilt library routes to the sibling rule; see search_bruker_lib for why this is a
        # separate rule rather than a conditional flag.
        _search = search_bruker_lib if getattr(cfg, "speclib", None) else search_bruker
        _lib_kw = {"speclib": cfg.speclib} if getattr(cfg, "speclib", None) else {}
        P.diann = _search(P, 
            P.raw,
            search_fasta=cfg.search_fasta,
            qvalue=cfg.qvalue,
            search_threads=cfg.search_threads,
            max_missed_cleavages=cfg.max_missed_cleavages,
            min_length=cfg.min_length,
            max_length=cfg.max_length,
            **_lib_kw
        )
        if getattr(cfg, "noise_real_data", False) or getattr(cfg, "spike_into", None):
            # Real-data background (A2 sampled blank, or spike-into a real run): the real peptides get
            # identified by the search and would inflate FDP. Render a background-only control (same seed,
            # no synthetic), search it, and subtract its IDs. The control reuses the render/search nodes —
            # a distinct command (adds `--noise-only`) → a separate content-addressed node.
            P.raw_control, P.truth_control = do_render(cfg, sample_id, P, control=True)
            P.diann_control = search_bruker(P, 
                P.raw_control,
                search_fasta=cfg.search_fasta,
                qvalue=cfg.qvalue,
                search_threads=cfg.search_threads,
                max_missed_cleavages=cfg.max_missed_cleavages,
                min_length=cfg.min_length,
                max_length=cfg.max_length,
            )
            P.score = score_bruker_bg(P, 
                P.diann, P.diann_control, P.truth, P.peptides, qvalue=cfg.qvalue
            )
        else:
            P.score = score_bruker(P, P.diann, P.truth, P.peptides, qvalue=cfg.qvalue)


def timsim_hye_quant_pipeline(P: Pipeline, cfg, sample_ids) -> None:
    """HYE quant (P0.2): render TWO conditions on the Bruker projector, search them in ONE joint DiaNN run,
    and score per-organism fold-change recovery. The feature-space nodes are shared with every other
    pipeline (same fingerprint); only the two renders differ by `sample_id` (each selects its condition's
    per-sample amounts from the shared `peptide_quantities`). `sample_ids` = [condition-A run, condition-B
    run] (e.g. ["A_R1", "B_R1"]); A is the reference. Opt-in via `--quant`. See HYE_QUANT.md."""
    sid_a, sid_b = sample_ids
    P.proteome = proteome(P, proteome_spec=cfg.proteome_spec)
    P.peptides, P.occurrences, P.cleavage_sites = digest(P, 
        P.proteome,
        max_missed_cleavages=cfg.max_missed_cleavages,
        min_length=cfg.min_length,
        max_length=cfg.max_length,
        max_peptides=cfg.max_peptides,
        seed=cfg.seed,
    )
    P.modforms, P.modifications = modify(P, P.peptides, mods_spec=cfg.mods, floor=cfg.floor)
    P.precursors = precursors(P, P.peptides, P.modforms, charge_model=cfg.charge_model, seed=cfg.seed,
                             min_charge_contrib=cfg.min_charge_contrib)
    P.ccs = ccs(P, P.precursors, P.peptides)
    P.rt = rt(P, P.peptides)
    P.samples, P.runs, P.sample_run_map, P.protein_quantities = design(P, 
        P.proteome, design_spec=cfg.design_spec
    )
    P.peptide_quantities, P.yield_report = peptide_yield(P, 
        P.proteome, P.occurrences, P.cleavage_sites, P.protein_quantities, P.modifications,
        digestion_efficiency=cfg.digestion_efficiency,
    )
    P.fragment_prediction_input = frag_input(P, P.precursors, P.peptides, P.modforms, P.modifications)
    P.fragment_intensities = fragments(P, 
        P.fragment_prediction_input, collision_energy=cfg.collision_energy, frag_model=cfg.frag_model
    )
    P.ion_spectra = spectra(P, P.precursors, P.peptides, P.modforms, P.modifications, P.fragment_intensities)

    def _render(sid):
        raw, _truth = do_render(cfg, sid, P)
        return raw

    P.raw_a = _render(sid_a)
    P.raw_b = _render(sid_b)
    if not getattr(cfg, "search_fasta", None):
        raise SystemExit("--quant requires --search-fasta (the joint DiaNN search is the quant measurement)")
    P.diann = search_bruker_joint(P, 
        P.raw_a, P.raw_b,
        search_fasta=cfg.search_fasta, qvalue=cfg.qvalue, search_threads=cfg.search_threads,
        max_missed_cleavages=cfg.max_missed_cleavages, min_length=cfg.min_length, max_length=cfg.max_length,
    )
    P.quant = quant(P, 
        P.diann, P.peptides, P.occurrences, P.proteome,
        design_spec=cfg.design_spec, qvalue=cfg.qvalue, delta=cfg.quant_delta,
    )


def timsim_bruker_dda_pipeline(P: Pipeline, cfg, sample_id: str) -> None:
    """One sample to a Bruker DDA-PASEF `.d` — identical feature-space chain to the Bruker DIA pipeline
    (so requesting DIA and DDA collapses the structure axis to one computation), but the measurement is
    `timsim-render --dda` (top-N selection) and phase-2 is Sage → `v2_dda_eval` (DiaNN is DIA-only).
    Opt-in via `--bruker-dda <ref.d>`."""
    P.proteome = proteome(P, proteome_spec=cfg.proteome_spec)
    P.peptides, P.occurrences, P.cleavage_sites = digest(P, 
        P.proteome,
        max_missed_cleavages=cfg.max_missed_cleavages,
        min_length=cfg.min_length,
        max_length=cfg.max_length,
        max_peptides=cfg.max_peptides,
        seed=cfg.seed,
    )
    P.modforms, P.modifications = modify(P, P.peptides, mods_spec=cfg.mods, floor=cfg.floor)
    P.precursors = precursors(P, P.peptides, P.modforms, charge_model=cfg.charge_model, seed=cfg.seed,
                             min_charge_contrib=cfg.min_charge_contrib)
    P.ccs = ccs(P, P.precursors, P.peptides)
    P.rt = rt(P, P.peptides)
    P.samples, P.runs, P.sample_run_map, P.protein_quantities = design(P, 
        P.proteome, design_spec=cfg.design_spec
    )
    P.peptide_quantities, P.yield_report = peptide_yield(P, 
        P.proteome,
        P.occurrences,
        P.cleavage_sites,
        P.protein_quantities,
        P.modifications,
        digestion_efficiency=cfg.digestion_efficiency,
    )
    # ── measurement branch: same feature-space nodes as Bruker DIA, then the DDA-PASEF projector ──
    P.fragment_prediction_input = frag_input(P, 
        P.precursors, P.peptides, P.modforms, P.modifications
    )
    P.fragment_intensities = fragments(P, 
        P.fragment_prediction_input, collision_energy=cfg.collision_energy, frag_model=cfg.frag_model
    )
    P.ion_spectra = spectra(P, 
        P.precursors, P.peptides, P.modforms, P.modifications, P.fragment_intensities
    )
    P.raw, P.truth = render_dda(P, 
        P.precursors,
        P.rt,
        P.ion_spectra,
        P.ccs,
        P.peptide_quantities,
        reference_d=cfg.reference_d,
        sample_id=sample_id,
        intensity_scale=cfg.intensity_scale,
        precursors_every=cfg.dda_precursors_every,
        max_precursors=cfg.dda_max_precursors,
        exclusion_width=cfg.dda_exclusion_width,
    )
    # ── phase 2 (opt-in): Sage-search the .d + score against the selection-event answer key ──
    if getattr(cfg, "search_fasta", None):
        P.sage = search_dda(P, P.raw, search_fasta=cfg.search_fasta)
        P.score = score_dda(P, P.sage, P.truth, P.peptides, qvalue=cfg.qvalue)


def timsim_sciex_pipeline(P: Pipeline, cfg, sample_id: str) -> None:
    """One sample to a SCIEX ZenoTOF SWATH **mzML** via the LEAN v2 projector (`timsim-render-sciex`) —
    imspy-free, no `.wiff`. Reuses the SAME `frag_input → fragments → spectra` feature-space chain as the
    Thermo/Bruker pipelines (so requesting all three collapses the structure to one computation), then
    projects onto a synthesised SWATH schedule. No CCS (SCIEX has no ion mobility)."""
    P.proteome = proteome(P, proteome_spec=cfg.proteome_spec)
    P.peptides, P.occurrences, P.cleavage_sites = digest(P, 
        P.proteome,
        max_missed_cleavages=cfg.max_missed_cleavages,
        min_length=cfg.min_length,
        max_length=cfg.max_length,
        max_peptides=cfg.max_peptides,
        seed=cfg.seed,
    )
    P.modforms, P.modifications = modify(P, P.peptides, mods_spec=cfg.mods, floor=cfg.floor)
    P.precursors = precursors(P, P.peptides, P.modforms, charge_model=cfg.charge_model, seed=cfg.seed,
                             min_charge_contrib=cfg.min_charge_contrib)
    P.rt = rt(P, P.peptides)
    P.samples, P.runs, P.sample_run_map, P.protein_quantities = design(P, 
        P.proteome, design_spec=cfg.design_spec
    )
    P.peptide_quantities, P.yield_report = peptide_yield(P, 
        P.proteome, P.occurrences, P.cleavage_sites, P.protein_quantities, P.modifications,
        digestion_efficiency=cfg.digestion_efficiency,
    )
    # ── measurement branch: same feature-space nodes as Thermo/Bruker, then the SWATH mzML projector ──
    P.fragment_prediction_input = frag_input(P, 
        P.precursors, P.peptides, P.modforms, P.modifications
    )
    P.fragment_intensities = fragments(P, 
        P.fragment_prediction_input, collision_energy=cfg.collision_energy, frag_model=cfg.frag_model
    )
    P.ion_spectra = spectra(P, 
        P.precursors, P.peptides, P.modforms, P.modifications, P.fragment_intensities
    )
    P.mzml, P.truth = render_sciex(P, 
        P.precursors,
        P.rt,
        P.ion_spectra,
        P.peptide_quantities,
        sample_id=sample_id,
        gradient_length_s=cfg.gradient_length_s,
        cycle_time_s=cfg.cycle_time_s,
        mz_min=cfg.mz_min,
        mz_max=cfg.mz_max,
        window_width=cfg.window_width,
        collision_energy=cfg.collision_energy,
        intensity_scale=cfg.intensity_scale,
        frag_model=cfg.frag_model,
    )
    # ── phase 2 (opt-in): DiaNN-search the mzML natively + score against the answer key ──
    if getattr(cfg, "search_fasta", None):
        P.diann = search_sciex(P, 
            P.mzml,
            search_fasta=cfg.search_fasta,
            qvalue=cfg.qvalue,
            search_threads=cfg.search_threads,
            max_missed_cleavages=cfg.max_missed_cleavages,
            min_length=cfg.min_length,
            max_length=cfg.max_length,
        )
        P.score = score_sciex(P, P.diann, P.truth, P.peptides, qvalue=cfg.qvalue)


def _parser() -> argparse.ArgumentParser:
    """The full CLI surface. Shared by `main()` and the job-TOML entry point `job()`, so a TOML-driven
    run (and any GUI on top of it) gets exactly the same knobs and defaults as the command line."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="/tmp/necro-timsim")
    ap.add_argument("--proteome-spec", default="hye.toml", help="multi-FASTA proteome spec")
    ap.add_argument("--mods", default="mods.toml", help="modification spec (e.g. mods_basic.toml for a light HeLa run)")
    ap.add_argument("--design-spec", default="design.toml", help="experiment design spec (e.g. design_hela.toml for a single-organism run)")
    # --- realism knobs (all default OFF, so nothing changes unless asked) -------------------
    ap.add_argument("--ion-count-noise", action="store_true",
                    help="Poisson counting statistics on the ANALYTE signal (distinct from "
                         "--noise-real-data, which injects real BACKGROUND). Roughens traces near "
                         "the floor; leaves bright peaks smooth, as counting does.")
    ap.add_argument("--instrument-cv", type=float, default=0.0,
                    help="common-mode spray/transmission fluctuation, as the TOTAL per-bin CV. "
                         "Intensity-independent, so unlike --ion-count-noise it roughens BRIGHT "
                         "peaks too (v2's form of v1 noise_frame/scan_abundance).")
    ap.add_argument("--run-rt-sd", type=float, default=0.0,
                    help="per-run RT displacement, SECONDS sd (v1 rt_variation_std)")
    ap.add_argument("--run-im-sd", type=float, default=0.0,
                    help="per-run mobility displacement, 1/K0 sd (v1 ion_mobility_variation_std)")
    ap.add_argument("--run-intensity-cv", type=float, default=0.0,
                    help="per-run intensity variation, CV (v1 intensity_variation_std)")
    ap.add_argument("--target-p", type=float, default=0.0,
                    help="v1's cumulative-probability peak truncation (0.999); 0 leaves --n-sigma in charge")
    ap.add_argument("--speclib", default=None,
                    help="reuse a prebuilt DIA-NN .speclib instead of regenerating it from the FASTA "
                         "on every run (~42 min / 2.1 GB for the human proteome). Deterministic from "
                         "the FASTA + digest settings, so reuse is equivalent.")
    ap.add_argument("--min-charge-contrib", type=float, default=0.0,
                    help="drop charge states below this share of a modform's charged mass (v1's "
                         "dia-PASEF config uses 0.25). NOTE this is on the STRUCTURE axis, so it "
                         "re-fingerprints precursors and everything downstream, including the "
                         "fragment predictions.")
    ap.add_argument("--samples", nargs="+", default=["A_R1", "B_R1"])
    ap.add_argument("--max-peptides", type=int, default=0,
                    help="cap the simulated peptides to this many (seeded sample; 0 = full analytic digest). "
                         "Keeps the full FASTA as the DiaNN search space — for a tractable run on a big proteome.")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--graph", help="write the DAG to this file")
    ap.add_argument("--thermo-template", help="build the Thermo .raw pipeline against this template")
    ap.add_argument("--bruker-reference", help="build the LEAN v2 Bruker .d pipeline (timsim-render) against "
                                               "this reference DIA .d — imspy-free, replaces the v1 `simulate` seam")
    ap.add_argument("--bruker-dda", help="build the Bruker DDA-PASEF .d pipeline (timsim-render --dda) against "
                                         "this reference .d — top-N selection, searched by Sage (not DiaNN)")
    ap.add_argument("--dda-precursors-every", type=int, default=10, help="DDA: MS1 survey every Nth frame")
    ap.add_argument("--dda-max-precursors", type=int, default=25, help="DDA: max precursors per MS2 (PASEF) frame")
    ap.add_argument("--dda-exclusion-width", type=int, default=25, help="DDA: dynamic-exclusion window (frames)")
    ap.add_argument("--sciex", action="store_true", help="build the LEAN v2 SCIEX ZenoTOF SWATH -> open mzML "
                                                         "pipeline (timsim-render-sciex) — imspy-free, no .wiff")
    ap.add_argument("--sciex-gradient-s", type=float, default=1800.0, help="SCIEX SWATH gradient length (s)")
    ap.add_argument("--sciex-cycle-s", type=float, default=3.0, help="SCIEX SWATH cycle time (s)")
    ap.add_argument("--sciex-mz-min", type=float, default=400.0, help="SCIEX SWATH window coverage min m/z")
    ap.add_argument("--sciex-mz-max", type=float, default=1200.0, help="SCIEX SWATH window coverage max m/z")
    ap.add_argument("--sciex-window-width", type=float, default=25.0, help="SCIEX SWATH isolation window width (Th)")
    ap.add_argument("--frag-model", default="", help="fragment model: '' (local timsTOF) or 'koina:Prosit_2020_intensity_HCD'")
    ap.add_argument("--collision-energy", type=float, default=25.0)
    ap.add_argument("--intensity-scale", type=float, default=5.0e5)
    # A1 signal-m/z noise (Bruker DIA render). ppm is v1's 3σ envelope; 6.5 == v1's real DIA config.
    # 0/0 (default) keeps the render byte-identical to the noiseless baseline. See REALISM_PLAN.md.
    ap.add_argument("--noise-mz-ppm", type=float, default=0.0,
                    help="A1: Gaussian m/z scatter on precursor (MS1) peaks, ppm 3σ envelope (v1 6.5). 0=off")
    ap.add_argument("--noise-frag-ppm", type=float, default=0.0,
                    help="A1: m/z scatter on fragment (MS2) peaks, ppm 3σ envelope. 0=off")
    ap.add_argument("--noise-mz-uniform", action="store_true",
                    help="A1: use v1's uniform m/z scatter (mz ± mz·ppm/1e6) instead of the default Gaussian")
    ap.add_argument("--noise-seed", type=int, default=0, help="seed for the (deterministic) noise draws")
    # A2 real-data background (Bruker DIA): sample real peaks from the reference .d (v1 add_real_data_noise).
    # v1's real DIA recipe runs A1 + A2 together; --noise-real-data + a nonzero --noise-mz-ppm reproduces it.
    ap.add_argument("--noise-real-data", action="store_true",
                    help="A2: inject real background peaks sampled from the reference .d onto the frames")
    ap.add_argument("--noise-precursor-frames", type=int, default=5,
                    help="A2: reference MS1 frames sampled per output precursor frame (v1 5)")
    ap.add_argument("--noise-fragment-frames", type=int, default=5,
                    help="A2: reference MS2 frames sampled per output fragment frame (v1 5)")
    ap.add_argument("--noise-intensity-max", type=float, default=150000.0,
                    help="A2: background intensity cap in absolute counts (v1 reference_noise_intensity_max)")
    ap.add_argument("--noise-precursor-fraction", type=float, default=0.2,
                    help="A2: keep probability per sampled MS1 background peak (v1 0.2)")
    ap.add_argument("--noise-fragment-fraction", type=float, default=0.2,
                    help="A2: keep probability per sampled MS2 background peak (v1 0.2)")
    ap.add_argument("--spike-into", default=None,
                    help="Spike-into-real (mode B): overlay the synthetic signal additively onto this REAL "
                         ".d (selects the Bruker DIA pipeline; the .d is also the reference geometry). Only "
                         "the spikes are labeled; the real run's IDs are subtracted from FDP via a control.")
    ap.add_argument("--search-fasta", default=None,
                    help="opt into phase 2: DiaNN-search the rendered .raw against this FASTA, then score "
                         "against the answer key. Omit to stop at the .raw.")
    ap.add_argument("--qvalue", type=float, default=0.01, help="DiaNN + scoring q-value / FDR threshold")
    ap.add_argument("--search-threads", type=int, default=16)
    ap.add_argument("--quant", action="store_true",
                    help="HYE quant (P0.2): render the two design conditions, JOINT-search them in one DiaNN "
                         "run, and score per-organism fold-change recovery. Needs exactly 2 samples "
                         "(--samples A_R1 B_R1), a HYE proteome (--proteome-spec hye.toml), and --search-fasta.")
    ap.add_argument("--quant-delta", type=float, default=0.5,
                    help="log2FC tolerance for the quant %%correct metric")
    ap.add_argument("--phospho", action="store_true",
                    help="Phospho site-localization (P1.3): DiaNN --monitor-mod search + FLR scoring vs the "
                         "render's true modform sites. Selects the Bruker DIA pipeline; use with "
                         "--mods mods_phospho.toml and --search-fasta.")
    ap.add_argument("--flr-target", type=float, default=0.01,
                    help="target false-localization rate for the FLR operating point")
    return ap


def build_cfg(a) -> SimpleNamespace:
    """Parsed args -> the config object every pipeline factory consumes."""
    return SimpleNamespace(
        proteome_spec=a.proteome_spec,
        max_missed_cleavages=2,
        min_length=7,
        max_length=30,
        max_peptides=a.max_peptides,
        mods=a.mods,
        floor=1e-3,
        charge_model="site-specific",
        design_spec=a.design_spec,
        digestion_efficiency=0.9,
        timsim_config="v1.toml",
        seed=41,
        template=a.thermo_template,
        frag_model=a.frag_model,
        collision_energy=a.collision_energy,
        intensity_scale=a.intensity_scale,
        noise_mz_ppm=a.noise_mz_ppm,
        peak_shape=getattr(a, "peak_shape", "per-peptide"),
        cycle_seconds=getattr(a, "cycle_seconds", 0.0),
        mobility_std_target=getattr(a, "mobility_std_target", 0.009),
        n_sigma=getattr(a, "n_sigma", 3.0),
        min_charge_contrib=getattr(a, "min_charge_contrib", 0.0),
        speclib=getattr(a, "speclib", None),
        # A bare flag, so it is "" (absent) or "--ion-count-noise " (present).
        ion_count_noise=bool(getattr(a, "ion_count_noise", False)),
        instrument_cv=getattr(a, "instrument_cv", 0.0),
        run_rt_sd=getattr(a, "run_rt_sd", 0.0),
        run_im_sd=getattr(a, "run_im_sd", 0.0),
        run_intensity_cv=getattr(a, "run_intensity_cv", 0.0),
        target_p=getattr(a, "target_p", 0.0),
        noise_frag_ppm=a.noise_frag_ppm,
        noise_mz_uniform=a.noise_mz_uniform,
        noise_seed=a.noise_seed,
        noise_real_data=a.noise_real_data,
        noise_precursor_frames=a.noise_precursor_frames,
        noise_fragment_frames=a.noise_fragment_frames,
        noise_intensity_max=a.noise_intensity_max,
        noise_precursor_fraction=a.noise_precursor_fraction,
        noise_fragment_fraction=a.noise_fragment_fraction,
        spike_into=a.spike_into,
        quant_delta=a.quant_delta,
        phospho=a.phospho,
        flr_target=a.flr_target,
        search_fasta=a.search_fasta,
        qvalue=a.qvalue,
        search_threads=a.search_threads,
        # spike-into implies the reference geometry comes from the spike target itself.
        reference_d=a.spike_into or a.bruker_reference or a.bruker_dda,
        # Bruker DDA-PASEF selection params
        dda_precursors_every=a.dda_precursors_every,
        dda_max_precursors=a.dda_max_precursors,
        dda_exclusion_width=a.dda_exclusion_width,
        # SCIEX SWATH schedule (synthesised — no .wiff)
        gradient_length_s=a.sciex_gradient_s,
        cycle_time_s=a.sciex_cycle_s,
        mz_min=a.sciex_mz_min,
        mz_max=a.sciex_mz_max,
        window_width=a.sciex_window_width,
    )


def select_build(a, ap: argparse.ArgumentParser | None = None):
    """Pick the pipeline factory implied by the flags (and validate coupled ones)."""
    def fail(msg):
        ap.error(msg) if ap else (_ for _ in ()).throw(SystemExit(f"error: {msg}"))
    if a.sciex:
        return timsim_sciex_pipeline
    if a.thermo_template:
        return timsim_thermo_pipeline
    if a.bruker_dda:
        return timsim_bruker_dda_pipeline
    if a.bruker_reference or a.spike_into or a.phospho:
        if a.phospho:
            if not (a.bruker_reference or a.spike_into):
                fail("--phospho needs a reference .d (--bruker-reference <dia.d>) to render onto")
            if not a.search_fasta:
                fail("--phospho needs --search-fasta (the DiaNN --monitor-mod search is the FLR measurement)")
        return timsim_bruker_v2_pipeline
    return timsim_pipeline


def job(P: Pipeline, config: dict) -> None:
    """necroflow **job-TOML entry point** — the config-file face of this flow, and the target any GUI
    should build against::

        ".pipeline" = "flow/timsim_flow.py:job"
        ".requests" = ["score"]
        bruker_reference = "/data/ref.d"
        proteome_spec    = "hela.toml"
        noise_mz_ppm     = 6.5
        noise_real_data  = true

    Keys are the CLI flags with dashes as underscores (`--noise-mz-ppm` → `noise_mz_ppm`); everything not
    given keeps its CLI default, so a TOML run and a command line run are the same run. `sample` picks the
    design sample (default: the first of `samples`). Unknown keys are rejected rather than silently ignored
    — a typo'd knob would otherwise produce a plausible-looking run of the wrong experiment."""
    a = _parser().parse_args([])
    unknown = [k for k in config if k not in vars(a) and k != "sample"]
    if unknown:
        raise SystemExit(f"unknown config key(s) {unknown}; valid keys: {sorted(vars(a))}")
    for k, v in config.items():
        if k != "sample":
            setattr(a, k, v)
    build = select_build(a)
    sample = config.get("sample", a.samples[0])
    build(P, build_cfg(a), sample)


def main() -> None:
    ap = _parser()
    a = ap.parse_args()
    cfg = build_cfg(a)
    build = select_build(a, ap)
    dag = DAG(a.outdir)
    if a.quant:
        # HYE quant: ONE cross-sample pipeline (two renders → joint search → fold-change), not a per-sample
        # loop. Require exactly two samples resolving to the two DESIGN conditions, with the design's
        # reference condition FIRST (it becomes DiaNN run 0 = A, the fold-change denominator). This guards
        # the reversed-samples / two-replicates-of-one-condition mistakes the scorer can't detect.
        if len(a.samples) != 2:
            ap.error(f"--quant needs exactly 2 samples (the two conditions), got {a.samples}")
        import tomllib
        _design = tomllib.load(open(a.design_spec, "rb"))
        _conds = [c["name"] for c in _design["condition"]]
        _ref = _design.get("design", {}).get("reference", _conds[0])
        if len(_conds) != 2:
            ap.error(f"--quant needs a 2-condition design; {a.design_spec} has conditions {_conds}")

        def _cond_of(sid):
            m = [c for c in _conds if sid == c or sid.startswith(c + "_")]
            return m[0] if len(m) == 1 else None
        _sc = [_cond_of(s) for s in a.samples]
        if None in _sc or set(_sc) != set(_conds):
            ap.error(f"--quant samples {a.samples} must map 1:1 to the two design conditions {_conds} "
                     f"(got {_sc}); name them e.g. {_conds[0]}_R1 {_conds[1]}_R1")
        if _sc[0] != _ref:
            ap.error(f"--quant: the FIRST sample must be the design reference condition {_ref!r} "
                     f"(the fold-change denominator); got {a.samples[0]} → {_sc[0]!r}. Reorder --samples.")
        P = Pipeline(dag)
        timsim_hye_quant_pipeline(P, cfg, a.samples)
        dag.require([P.quant])
        print(f"\n  HYE quant conditions   : {a.samples[0]} (ref) vs {a.samples[1]}")
        print(f"  nodes actually to run  : {len(dag.nodes)}   <- structure deduplicated by fingerprint")
        if a.graph:
            dag.save(a.graph)
            print(f"  -> {a.graph}")
        if a.dry_run:
            print("\n  resolved commands:")
            for node in dag.nodes:
                cmd = resolve_command(node)
                if cmd:
                    print(f"    {cmd}")
            return
        dag.execute()
        return
    n_across = 0
    for sid in a.samples:
        P = Pipeline(dag)
        build(P, cfg, sid)
        # For the Thermo branch the answer key + run manifest are first-class deliverables (co-outputs of
        # the render), so request them explicitly alongside the .raw.
        req = [P.mzml] if getattr(P, "mzml", None) is not None else [P.raw]
        if getattr(P, "truth", None) is not None:
            req.append(P.truth)
            # The Thermo render also co-emits a run manifest; the lean Bruker render does not.
            if getattr(P, "manifest", None) is not None:
                req.append(P.manifest)
        # Phase 2 (opt-in): the score is the terminal deliverable — requesting it pulls search + the .raw.
        if getattr(P, "score", None) is not None:
            req.append(P.score)
        dag.require(req)
        n_across += len(P.nodes)

    # The claim this whole exercise rests on, stated as a number rather than an argument. Nodes are interned
    # on the shared DAG, so `len(dag.nodes)` is the deduplicated count across every sample's pipeline.
    n_unique = len(dag.nodes)
    print()
    print(f"  samples requested      : {len(a.samples)}")
    print(f"  nodes across pipelines : {n_across}")
    print(f"  nodes actually to run  : {n_unique}   <- structure deduplicated by fingerprint")
    print(f"  saved                  : {n_across - n_unique} redundant stage executions")

    if a.graph:
        dag.save(a.graph)
        print(f"  -> {a.graph}")

    if a.dry_run:
        print()
        print("  resolved commands:")
        for node in dag.nodes:
            cmd = resolve_command(node)
            if cmd:
                print(f"    {cmd}")
        return

    dag.execute()


if __name__ == "__main__":
    main()
