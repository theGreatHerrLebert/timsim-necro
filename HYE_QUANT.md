# P0.2 — HYE quant + fold-change eval

The second-most-cited v1 axis and a genuinely NEW eval dimension: a Human/Yeast/E.coli mixed proteome at two
known per-organism dilutions, scored on whether the search engine **recovers the known log2 fold-changes**
(and how much cross-species interference leaks). Unlike recall/FDP (per-sample ID), this is a **cross-sample
quant** metric that generalizes to every instrument.

## What already exists (most of it)

- **Two conditions with dilutions** — `flow/configs/design.toml`: A = {H .65, Y .30, E .05}, B = {H .65,
  Y .20, E .15 ("rest")}. Expected B/A ratios: **HUMAN 1.0 (log2 0)**, **YEAST 0.667 (log2 −0.585)**,
  **ECOLI 3.0 (log2 +1.585)**.
- **Per-sample fan-out** — `--samples ["A_R1","B_R1"]`; the loop builds a pipeline per sample, each renders
  its own `.d` + searches (DiaNN). Structure nodes dedupe; only render/search/score differ per sample.
- **Organism** is a column in `proteome.parquet`; peptide→protein is in `occurrences.parquet`; DiaNN
  `Precursor.Quantity` is parsed to `intensity`. So peptide→organism→intensity is all reachable.

## Gaps (the actual work)

### 1. FASTA repoint (config)
`hye.toml`'s three `[[source]]` paths point at the cleared `SUBMISSION/zenodo/...` tree. Fix: split
`/scratch/timsim-demo/PlasmaBENCH/data/raw/fasta/ute_hye_30867.fasta` (31k seqs, UniProt `_HUMAN/_YEAST/
_ECOLI` headers) into `HUMAN.fasta`/`YEAST.fasta`/`ECOLI.fasta` under `flow/configs/hye/`, repoint the three
sources. (Keeping three tagged sources matches the existing `organism = "..."` declaration.)

### 2. Quant-accuracy scorer — `timsim_eval/v2_quant_eval.py` (NEW)
Inputs: the two conditions' DiaNN reports (A, B), `occurrences` + `proteome` (organism map), and the
expected per-organism ratios. Steps:
- **Organism map:** join peptide→protein (occurrences)→organism (proteome). A peptide shared across
  organisms is **ambiguous → excluded** (v1's "Unknown" drop; report the count).
- **Per-peptide (or per-protein) ratio:** for peptides quantified in BOTH conditions, `log2(intB / intA)`.
  Precursor→peptide roll-up first (sum or MaxLFQ-style; start with summed precursor intensity per peptide).
- **Per-organism metrics vs the design's expected log2FC:**
  - **bias** = median(observed log2FC) − expected (how far off-center);
  - **precision** = MAD/IQR of observed log2FC (spread — interference widens it);
  - **% correctly regulated** (sign + within a tolerance band);
  - **cross-species leakage:** HUMAN is the 1:1 anchor — its observed spread IS the interference floor;
    YEAST/ECOLI pulled toward 0 (1:1) = ratio compression from co-isolation. Report per-organism
    distributions + a leakage summary (e.g. fraction of yeast/ecoli peptides mis-centered toward human).
- Output JSON + a text summary (per-organism n, median log2FC, expected, bias, spread, %correct).

### 3. Cross-sample quant node in the flow (the wiring challenge)
Today each sample is an independent pipeline in `for sid in a.samples`. Quant needs ONE node depending on
**both** conditions' reports. Plan:
- Build the per-sample pipelines as now, **collect each `P.diann`** (+ shared `occurrences`, `proteome`).
- Add a `quant` node `r.quant(diann_a, diann_b, occurrences, proteome, expected, qvalue)` →
  `python -m timsim_eval.v2_quant_eval ...`. Both `diann_*` are typed `DiannReport` inputs from the two
  sample searches (distinct `.d` → distinct fingerprints), so the DAG orders it after both.
- Gate it on a `--quant` flag (or auto when ≥2 samples map to ≥2 conditions). Request the quant node as a
  terminal deliverable so it pulls both searches.
- **Conditions vs replicates:** first cut = 1 replicate/condition (A_R1 vs B_R1). Generalize later to
  average replicates per condition before the ratio (the node takes lists per condition).

### 4. Expected ratios (source of truth)
Parse `design.toml`'s conditions → per-organism mix, compute B/A. Simplest: a tiny helper reads the two
conditions' `mix` and emits `{HUMAN: log2(0.65/0.65), YEAST: log2(0.20/0.30), ECOLI: log2(0.15/0.05)}`.
Pass to the scorer (so the eval isn't hard-coded to one design). "rest" resolves to `1 − Σ others`.

## Review resolutions (codex)

- **JOINT DiaNN run over both `.d` (architectural — the big one).** Independent per-sample searches do NOT
  give cross-run quant (DiaNN's normalization + MaxLFQ + MBR operate over an *experiment*, not one run). A
  post-hoc node can't recover that. So the quant path searches **both conditions' `.d` in ONE DiaNN run**
  (`--f A.d --f B.d`) → one report with per-run quantities. This joint search IS the cross-sample join (and
  it supersedes the "post-hoc quant node over two reports" idea): render A/B → one `search_quant` node →
  `v2_quant_eval`. Record the DiaNN normalization setting/version in the output.
- **Normalization — do NOT human-anchor the primary score** (it forces HUMAN residual ≡ 0, launders a global
  error into yeast/ecoli, and makes leakage circular). Emit **three labelled views**: (1) raw
  `Precursor.Quantity` ratio, (2) **`Precursor.Normalised` ratio = the primary engine-performance result**,
  (3) human-anchored calibration as a *diagnostic only* (never the score). Note the HYE design is
  deliberately compositional, which violates DiaNN's "most peptides unchanged" normalization assumption —
  hence reporting all three.
- **Metric:** peptide-level `log2(B/A)` on organism-unique sequences is the **primary measurement-level**
  endpoint (exposes interference protein roll-up hides); add **`PG.MaxLFQ` protein-level as a SECONDARY**
  endpoint (don't replace). Summary = **median residual `median(log2FC − truth)` + MAD** (robust). `%correct`
  uses an explicit tolerance `|log2FC − truth| ≤ δ` (HUMAN has no sign, so band-not-sign). Peptides are not
  independent replicates — bootstrap CIs are descriptive only.
- **"Leakage" is ratio COMPRESSION, not proven interference.** Pulled-toward-0 can come from normalization,
  censoring, or nonlinearity too. Two options: (a) v1-style — call it **"cross-species-associated ratio
  compression"** and report per-organism residual distributions; (b) **the simulator upgrade** — because we
  own the truth, compute the *actual* fraction of each precursor's extracted signal attributable to another
  organism (co-isolation in the same window/RT) and correlate it with the FC residual. Ship (a) first; (b)
  is the honest-leakage upgrade (needs per-peak organism provenance through the render).
- **Coverage tables (report, don't hide).** Complete-case (quantified in BOTH) is the primary estimand, but
  it induces survivorship bias (low-abundance ECOLI absent in A ⇒ ECOLI FC looks tighter/better). Emit a
  per-organism **eligibility table** (theoretical / organism-unique / quantified-in-both / scored) and a
  **detection table** (A-only / B-only / both, by abundance). Do NOT impute extreme ratios into bias/MAD.
- **Gate on explicit `--quant`;** validate the samples resolve to **exactly two named conditions** (≥1 run
  each) — never auto-pick an arbitrary A/B. Design the interface for replicate aggregation from the start
  (average replicate quantities per condition before the ratio), even if the first cut runs 1 run/condition.
- **Sequence key + aggregation:** filter at the precursor q-value; treat zero/absent quantity as
  below-quantification (not 0); key on stripped peptide (+ charge/mod handled by the roll-up). Peptide-level
  quantity = summed **normalized** precursor quantities across charge/mod forms (state the rule); MaxLFQ only
  for the protein endpoint. Exclude ambiguous (>1 organism) peptides; no razor for the primary metric.
- **Validation upgrade:** besides no-noise, run a **no-interference counterfactual** (e.g. a single-organism
  or well-separated render) — the clean test that the compression metric tracks interference, not ordinary
  quant error. Note: in a pure sim there's no run-order/instrument drift, so batch confounding is N/A until
  we add drift.

## First-cut scope (what to build now)

1. FASTA repoint. 2. A joint `search_quant` node (`--f A.d --f B.d`) → one report. 3. `v2_quant_eval`:
organism-unique peptide-level `log2(B/A)` from `Precursor.Normalised`, complete-case, per-organism median
residual + MAD + `%correct(δ)`, the eligibility + detection tables, and the three normalization views;
MaxLFQ protein endpoint + true-sim-leakage deferred as documented upgrades. 4. `--quant` gate + 2-condition
validation. 5. expected ratios parsed from `design.toml`.

## Original open questions (now resolved above)

1. **Precursor→peptide→protein roll-up:** start with summed precursor `intensity` per peptide, ratio at the
   peptide level, then aggregate peptides→organism? Or protein-level (MaxLFQ)? v1 used protein-level
   fold-changes. Peptide-level is simpler and adequate for a first cut — is it faithful enough, or do we
   need protein roll-up to match v1's numbers?
2. **Missing-in-one-condition peptides:** a peptide seen in A but not B (or vice-versa) — drop from the
   ratio, or treat as a large/small ratio? v1 requires both. Drop + report the count seems right.
3. **Normalization:** DiaNN quantities may need median-normalization across runs before ratios (a global
   loading/normalization factor). HUMAN at 1:1 is the natural anchor — normalize so median HUMAN log2FC = 0?
   Or trust DiaNN's own normalization? This materially affects yeast/ecoli bias.
4. **Cross-sample node in necroflow:** is referencing two sample pipelines' `P.diann` in one added node
   supported cleanly, or does the per-sample `dag.add(P)` loop need restructuring so the quant node is
   added once with both handles?
5. **Ambiguous peptides:** exclude any peptide mapping to >1 organism. Confirm that's the right call vs
   protein-unique-peptide-only (razor). Start with organism-unique.
6. **Which quant column:** `Precursor.Quantity` vs `PG.MaxLFQ`/`Genes.MaxLFQ` from DiaNN. MaxLFQ is the
   quant-grade column; Precursor.Quantity is rawer. Start with what `parse_diann_report` already gives
   (Precursor.Quantity), note MaxLFQ as the upgrade.

## Validation
- FASTA repoint: `timsim-proteome` on `hye.toml` tags H/Y/E correctly (counts per organism sane).
- Scorer on a rendered A/B pair: HUMAN median log2FC ≈ 0, YEAST ≈ −0.585, ECOLI ≈ +1.585 (within the
  interference-driven spread); cross-species leakage reported. A no-noise render should recover the ratios
  tightly; A2/spike noise should widen them toward realism (ties back to the realism track).
- Deterministic given the same renders + reports.
