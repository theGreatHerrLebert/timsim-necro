# timsim-necro

The **timsim v2 simulator**, as a lean orchestration repo: a [necroflow](https://github.com/MatteoLacki/necroflow)
DAG that assembles a synthetic timsTOF/Thermo/SCIEX DIA run — and **ingests only the pieces it actually needs**,
not the old `imspy` monorepo.

## Why this exists

The v1 simulator shipped as one monolith: to run a simulation you installed the entire `imspy` stack
(`imspy-core`, `imspy-predictors`, `imspy-simulation`, `imspy-dia`, `imspy-search`, `imspy-vis`) plus the whole
Rust tree. The prediction step alone dragged PyTorch, Koina, the search engine, and the v1 simulator into every
install. `timsim-necro` is the opposite: the DAG declares its dependencies, and each is a small, independently
versioned, composable package.

```
necroflow (typed DAG framework)                         ← orchestration
      │  timsim_flow.py wires the nodes
      ▼
┌─────────────── the simulator's steps ───────────────┐
│ structure/render (Rust bins):  timsim-proteome/digest/design/precursors/
│                                yield/modify/frag-input/spectra/  render-thermo (.raw) /
│                                render (Bruker .d DIA, v2) / render --dda (Bruker DDA-PASEF .d) /
│                                render-sciex (SWATH mzML, v2)
│ prediction (Python, LEAN):     timsim-ccs / timsim-rt / timsim-fragments
│                                   └─ timsim-predict → pepdl → mscorepy
│                                      (mscore + ms-chem pyo3 primitives; imspy-free)
│ search (external):             DiaNN on every DIA path; Sage on the DDA path (DiaNN is DIA-only)
│ score (Python, LEAN):          timsim-eval (parse report → compare to truth → metrics)
└──────────────────────────────────────────────────────┘
```

**No rustims, no imspy monorepo — everything is ingested from small federated repos.** The prediction step — the one that used to pull the heaviest dead weight — is now
[`timsim-predict`](https://github.com/theGreatHerrLebert/timsim-predict) →
[`pepdl`](https://github.com/theGreatHerrLebert/pepdl) →
[`mscorepy`](https://github.com/theGreatHerrLebert/mscore), with **zero imspy** in its closure. CCS reproduces
the old output byte-for-byte (40,509 precursors, 0 diff); RT is Chronologer (Searle Lab); fragments are Prosit/
local. Torch is optional (`[koina]` runs remote, torch-free; `[local]` adds torch + Chronologer + the on-device
intensity model).

## Quickstart

```bash
make predict-deps         # necroflow + timsim-predict → pepdl → mscorepy   (the lean prediction stack)
make rust-bins            # cargo install timsim-cli from git (crates.io deps only, no rustims)
```

**Two things the flow needs from you, and both are easy to miss.**

1. **`TIMSIM_BIN` must point at the directory holding the Rust binaries.** The flow resolves them as
   `$TIMSIM_BIN/timsim-*` and defaults to `target/release`, which is *not* where `make rust-bins` puts
   them: the Makefile runs `cargo install --root $(PREFIX)` with `PREFIX ?= $(HOME)/.local`, so the
   binaries land in **`$HOME/.local/bin`**. Export it or every node fails with "command not found".
2. **Run from `flow/configs/`.** The spec defaults (`--proteome-spec hye.toml`, `--mods mods.toml`,
   `--design-spec design.toml`) are *bare filenames* resolved against the current working directory, and
   they only exist in `flow/configs/`. Absolute paths work from anywhere; `ramp/ramp.py` sets
   `cwd=flow/configs` for exactly this reason. The Python prediction jobs (`timsim-ccs`, `timsim-rt`,
   `timsim-fragments`) are invoked by bare name, so they must be on `PATH`.

A complete, working invocation — a HeLa Bruker-DIA run with v1's noise recipe, searched and scored:

```bash
export TIMSIM_BIN="$HOME/.local/bin"
cd flow/configs
python ../timsim_flow.py \
    --bruker-reference /path/to/reference-dia.d \
    --proteome-spec hela_proteome.toml --mods mods_basic.toml --design-spec design_hela.toml \
    --samples A_R1 --max-peptides 5000 \
    --noise-mz-ppm 6.5 --noise-frag-ppm 6.5 --noise-real-data \
    --search-fasta hela_subset.fasta \
    --outdir /tmp/necro-hela --dry-run          # drop --dry-run to execute
```

`python ../timsim_flow.py --help` lists every knob. Omit `--search-fasta` to stop at the rendered `.d`.

### Or drive it from a job TOML

The same flow exposes a necroflow **job-TOML entry point** (`flow/timsim_flow.py:job`) — the config-file
face, and the target a GUI builds against. `flow/configs/job_example.toml` is a complete example: keys are
the CLI flags with dashes as underscores, anything omitted keeps its CLI default, and an unknown key is
rejected rather than silently ignored.

```bash
export TIMSIM_BIN="$HOME/.local/bin"
cd flow/configs
necroflow run job_example.toml --outdir /tmp/necro-job --dry-run
```

### The Python surface

`requirements.txt` is all of it — **three git dependencies** (`necroflow`, `timsim-predict`, `timsim-eval`),
nothing from imspy. **necroflow must be ≥ 0.0.4** (the flow uses the module-level rule API). A PEP 508
direct-URL requirement cannot carry a version specifier, so that floor is documented rather than enforced
by pip; the tracked branch is at or above it, and `pip install 'necroflow>=0.0.4'` from PyPI is the
enforceable alternative.

## What's fully lean vs. still coupled

**Fully federated ✅**
- **Prediction** — `timsim-predict` (git) → `pepdl` → `mscorepy`. Independently installable, imspy-free,
  validated end-to-end (CCS exact parity, RT Chronologer, fragments local).
- **Orchestration** — `necroflow` (git), imports nothing from this project's internals.

- **Rust protocol/render tools** — [`timsim-cli`](https://github.com/theGreatHerrLebert/timsim-cli), its own
  repo, `cargo install --git`-able, depending only on published crates. **No rustims.**

- **Eval / validation** — [`timsim-eval`](https://github.com/theGreatHerrLebert/timsim-eval), its own repo.
  The SCORE node (`timsim_eval.v2_thermo_eval`) parses the DiaNN report and compares it to the render's
  ground-truth manifest. Pure-Python, **imspy-free** on the DiaNN path — the last imspy touchpoint is cut.

**Render backends**
- **Thermo `.raw`** (`--thermo-template`) — lean: `frag_input → fragments → spectra → render-thermo`
  (timsim-cli). Co-emits the answer key + manifest → the phase-2 DiaNN `search`/`score` closes on it.
- **Bruker `.d`** (`--bruker-reference <ref.d>`) — **lean v2, fully closed**: the same feature-space chain
  plus CCS → `timsim-render`, a streaming imspy-free projector onto a reference `.d`'s DIA schedule, which
  **co-emits the per-precursor answer key** (same 8-column schema as Thermo). With `--search-fasta` the
  DAG appends `search_bruker` (DiaNN reads the `.d` **natively** — no .NET) + `score_bruker`, so Bruker
  closes structure → render → search → score just like Thermo. Verified end-to-end: a 60-protein run
  authors a valid 3000-frame DIA `.d` (177 MS1 + 2823 MS2), a 42,919-row truth, DiaNN searches it as
  Slice-PASEF, and the scorer reports a monotonic recall-by-abundance ladder.
- **Bruker DDA-PASEF `.d`** (`--bruker-dda <ref.d>`) — **lean v2**: same feature-space chain + CCS →
  `timsim-render --dda` (MS1 surveys + top-N precursor selection with dynamic exclusion + band-limited MS2),
  co-emitting a per-**selection-event** answer key. Searched by **Sage** (not DiaNN — which is DIA-only):
  `search_dda` (Sage reads the `.d` natively) + `score_dda` (`v2_dda_eval` maps Sage's PSMs to the
  fragmented precursors; recall is *conditional* on the top-N DDA selected). Verified end-to-end: 6,557
  correct PSMs, FDP 0.17%.
- **SCIEX mzML** (`--sciex`) — **lean v2**: the same feature-space chain → `timsim-render-sciex`, which
  projects onto a **synthesised SWATH schedule** and writes open **mzML** via `timsim-core` (mzdata) — no
  `.wiff`, no `sciexwiff`/`sciex-io` (legally clean). Co-emits the answer key; with `--search-fasta` the
  DAG appends `search_sciex` (DiaNN reads open mzML natively) + `score_sciex`, so SCIEX closes
  structure → render → search → score like the others. Verified end-to-end (fresh sim through necro):
  84.4% detectable recall, FDP 0.81%, monotonic ladder — on par with Thermo. (Native `.wiff` output is a
  separate rustims-local satellite reusing the validated `sciexwiff` writer.)

**Thermo, Bruker DIA, Bruker DDA, and SCIEX all close structure → render → search → score** on small,
independently-versioned federated repos — zero imspy, zero rustims. The one remaining v1 cord is the
**native SCIEX `.wiff`** writer (a rustims-local satellite by design, since `sciexwiff` is legal-held; the
open **mzML** SCIEX path is fully lean).

## Realism: noise and spike-into-real

A noiseless render makes FDP trivially small, so the realism track (design: `REALISM_PLAN.md`) ports v1's
**two distinct** noise axes — v1's real DIA recipe runs one of each, together:

- **A1 — signal m/z scatter** (`--noise-mz-ppm`, `--noise-frag-ppm`; ppm is a 3σ envelope; `--noise-seed`;
  `--noise-mz-uniform` for v1's uniform draw instead of the default Gaussian). Moves the *simulated* peaks.
  `0` — the default — keeps the render byte-identical to the noiseless baseline, so the reproducibility
  gate still holds.
- **A2 — real-data background** (`--noise-real-data`, plus `--noise-precursor-frames`,
  `--noise-fragment-frames`, `--noise-intensity-max`, `--noise-precursor-fraction`,
  `--noise-fragment-fraction`). Samples real peaks out of the reference `.d` and injects them onto the
  rendered frames — background belonging to no simulated ion.

v1 parity is **A1 + A2 at once**: `--noise-mz-ppm 6.5 --noise-frag-ppm 6.5 --noise-real-data`.

**The noise-only control.** With A2 the DAG automatically builds a *second* render of the background alone
(`--noise-only`, same seed) and searches that too. Real peptides in the reference blank do get identified,
and they are not false positives against a *synthetic* answer key — so `score_bruker_bg` subtracts the
control's IDs and reports a **background-subtracted FDP**. Both views come out of the one `metrics.json`:
`fdp_sub = false / diann_ids` (the reported, background-adjusted number) and
`fdp_raw = (false + background_subtracted) / (diann_ids + background_subtracted)` (the end-to-end one).

**Spike-into-real** (`--spike-into <real.d>`) — the strongest realism mode: overlay the synthetic
ground-truth signal *additively* onto a **real** experimental `.d`, which also supplies the reference
geometry. Only the spikes are labelled; the real run's own identifications are removed from FDP through
the same control render.

## Benchmarks built on top

- **HYE quant** (`--quant`; design: `HYE_QUANT.md`) — a Human/Yeast/E. coli mixture at two known
  per-organism dilutions. Renders both design conditions, **joint-searches them in one DiaNN run**
  (`search_bruker_joint`), and scores per-organism log2 fold-change recovery against the ratios *derived*
  from the mixture (H 0, Y −0.585, E +1.585): bias, spread, and the fraction within `--quant-delta`.
  Needs exactly two samples mapping 1:1 onto the two design conditions with the reference condition first,
  plus `--search-fasta`; the flow validates all of that before it builds the graph.
- **Phospho / FLR** (`--phospho`; design: `PHOSPHO_FLR.md`) — site-localization scoring. Both positional
  isomers are simulated *in one run* (co-eluting), DiaNN searches with `--monitor-mod`, and `score_flr`
  compares the localized site to the render's true `mod_positions` at a `--flr-target` operating point.
  Use with `--mods mods_phospho.toml`; `mods_phospho_reg.toml` adds per-residue occupancy variation, which
  is what makes the benchmark meaningful — at uniform occupancy the two isomers are equally abundant by
  construction, so localization is ambiguous and FLR measures nothing.
- **Complexity ramp** (`ramp/ramp.py`; design: `RAMP.md`) — the empirical-FDP-vs-peptide-density curve.
  Drives the flow once per (density level × noise condition), factorial over `noiseless / A1 / A1+A2`, and
  aggregates each run's `metrics.json`. Committed results live in `ramp/results/` (`ramp.json`, `ramp.md`,
  `ramp.svg`): three `--max-peptides` levels (5k / 20k / 80k → 29,250 / 116,279 / 466,842 rendered
  precursors) at q = 1%. Read them before repeating the experiment — FDP does *not* climb with density
  here (1.22–1.45% across the noiseless levels); what shows up instead is A2's background, lifting raw FDP
  to 4.49% at the sparsest level and washing out as density rises.
- **`--max-peptides N`** — the density knob used throughout: cap the simulation to a seeded sample of the
  analytic digest (`0` = no cap) while keeping the full FASTA as the search space. It is a seeded-shuffle
  *prefix*, so for a fixed seed the levels are nested — a clean superset ladder.

## Layout
```
flow/timsim_flow.py     the DAG (nodes, typed edges, command wiring) + the `job` TOML entry point
flow/configs/           run configs (design*.toml, hela*.toml, hye.toml, sciex.toml, mods*.toml, fasta,
                        job_example.toml) — and the working directory the flow's defaults expect
golden/                 the regression gate — a frozen 60-protein run scored two ways (sim realism +
                        tool benchmark) against its own answer key; `golden/run.sh`, see golden/README.md
ramp/                   the complexity-ramp driver (ramp.py) + its committed results (ramp/results/)
docs/v2-design/         the v2 design record: spec, plan, render + MS2 designs, and their codex reviews
requirements.txt        the lean Python dependency surface (necroflow + timsim-predict + timsim-eval)
Makefile                predict-deps / rust-bins / setup

PORTING_ROADMAP.md      what of v1 to port, what to drop, and why — the tiering that drives the rest
REALISM_PLAN.md         P0   — the noise model (A1 / A2) + spike-into-real
HYE_QUANT.md            P0.2 — HYE quant + fold-change eval
PHOSPHO_FLR.md          P1.3 — phospho site-localization + FLR
RAMP.md                 P1.4 — the complexity ramp and the FDP-vs-density curve
```
