# P1.4 — HeLa complexity ramp + true-FDR-vs-density curve

The capstone of the realism track. Render HeLa at increasing peptide DENSITY, **with noise on**, search +
score each, and plot **true-FDR (empirical FDP vs the answer key) against density**. The point (roadmap):
as density rises the search engine makes more co-isolation errors, so true-FDR climbs **above the nominal
q-value** — an inflation that is "only meaningful after noise" (which we now have: A1 signal-m/z + A2
real-data background).

## Why it's mostly orchestration (machinery already exists)
- **Render at a density** = `--max-peptides N` (the digest cap; the direct density knob).
- **Noise on** = `--noise-mz-ppm 6.5 --noise-real-data` (the v1 DIA recipe). With A2 the flow already runs a
  noise-only control + subtracts background IDs (score_bruker_bg) — so the reported FDP is the *genuine*
  search error rate, not inflated by real blank-peptide IDs. That subtraction is CORRECT here: the FDR
  *inflation* we want is wrong matches on the SYNTHETIC signal (co-isolation), which are false synthetic IDs,
  NOT background IDs.
- **true-FDR** = the score node's `fdp` (= false / diann_ids vs the answer key); **recall** = its `hierarchy`
  (detectable level). Both already emitted per run.

So the ramp = run the flow at each level → collect each `metrics.json` → aggregate `(density, fdp, recall,
nominal_q)` → a table + the true-FDR-vs-density curve.

## Design
- **Driver** `ramp/ramp.py`: for each level, invoke the flow (venvA python, TIMSIM_BIN, `--bruker-reference`
  <blank.d>, `--mods mods_basic.toml`, `--noise-mz-ppm 6.5 --noise-real-data`, `--search-fasta`,
  `--max-peptides N`, `--samples A_R1`), then read the score `metrics.json`. Aggregate to `ramp/ramp.json`
  + a markdown table; emit a self-contained SVG/HTML curve (true-FDR + recall vs density, with the nominal
  q-value line).
- **Two tiers** (compute-bounded):
  - **Validation ramp** (do now): `hela_subset.fasta` (2500 proteins) → levels ~[2k, 8k, 20k, 40k]
    max-peptides. DiaNN on 2500 proteins is FAST, so the whole sweep is tractable and shows the trend.
  - **Full ramp** (documented long-run): full human proteome (`HUMAN.fasta`, 20420 proteins) → up to 250k.
    DiaNN lib-gen on 20k proteins is slow (hours) — a background sweep, not interactive.
- **Density axis:** report the ACTUAL rendered precursor / peptide count (from the truth), not just the
  `--max-peptides` request (sampling + charge/mods shift it).
- **Fixed:** q-value 1%, the same reference `.d`, the same noise seed — so only density varies.

## Review resolutions (codex)

- **Terminology: empirical FDP, not "true FDR"** (FDR is an expectation; a run realizes FDP). The claim is a
  **trend with confidence bounds separating from the nominal 1%**, NOT run-by-run monotonicity. Multiple
  render/noise **seeds** per density point → mean/median + bootstrap/binomial CI. First cut: 1 seed (compute),
  driver loops seeds.
- **Show BOTH FDP views** (subtraction hides part of the inflation A2 causes): **unsubtracted empirical FDP =
  the PRIMARY end-to-end calibration result**; background-adjusted (subtracted) FDP = a labelled secondary
  attribution. Both are derivable from the one metrics.json: `fdp_sub = false/diann_ids` (reported);
  `fdp_raw = (false + background_subtracted)/(diann_ids + background_subtracted)` — no scorer change.
- **Factorial attribution: A1-only vs A1+A2 × density.** A1-only = density/interference without chemical
  background; `(A1+A2) − (A1-only)` = A2's contribution. Two curves, not one.
- **Confounds — controllable here:** (a) `--max-peptides` is a **seeded-shuffle-prefix → NESTED** set for a
  FIXED seed (level N ⊂ level N'), so it's a clean superset ladder. (b) **No per-peptide dilution**: the
  design's load is per-PROTEIN (hockeystick over `n_proteins`), and `--max-peptides` subsamples the digest
  independently, so a peptide keeps its protein-derived amount. Set **`n_proteins` = the full proteome** so
  every sampled peptide's protein is expressed (no expression confound). Hold seed, reference `.d`,
  acquisition, q, abundance model fixed — only `--max-peptides` varies.
- **Density axis:** report the **actual rendered precursor count** AND the **detectable count** (present &
  in-window & has-frags & abundance>floor) AND a **co-elution proxy** (mean precursors per RT-bin from the
  truth) — global count alone is a weak proxy for DIA overlap.
- **"Only meaningful after noise" is too strong:** empirical FDP is informative noiselessly; noise
  strengthens the *realism* claim. Include a noiseless anchor per density (cheap, no A2 control).
- **Answer-key UNIT check:** the truth is keyed (sequence, charge) — confirm DiaNN reporting matches
  (precursor/charge, I→L); a unit mismatch manufactures apparent inflation. (Already the scorer's contract.)

## Original open questions (resolved above)
1. Is `--max-peptides` the right density knob, or should density be precursors (post-charge/mods)? Report
   actual precursor count on the x-axis regardless.
2. With A2 background-subtraction ON, is the reported FDP the right "true-FDR" for the inflation story (yes:
   it isolates synthetic-signal search errors from real blank IDs) — or should we ALSO show the un-subtracted
   FDP to see the background's contribution? Show both?
3. Does density-driven inflation need A2, or does A1 + density alone suffice (co-isolation is a density
   effect independent of chemical background)? Run A1+A2 (the recipe); consider an A1-only control curve to
   attribute the inflation.
4. Expected shape: true-FDR flat ≈ nominal at low density, rising above it as density climbs. Recall should
   fall (more interference). Is there a saturation/ceiling to watch for (u32 intensity, DiaNN's own FDR
   control kicking in)?
5. Compute: the A2 control double-renders + double-searches each level. Acceptable, or run A1-only for the
   ramp and A2 at one anchor point?

## Validation
- The driver reproduces the existing 5K point (FDP ~1.4%) when run at that density noiseless.
- The curve is monotone-ish: true-FDR rises with density under noise; recall falls.
- Deterministic given the seed; the aggregation is pure (re-runnable from the per-level metrics.json).
