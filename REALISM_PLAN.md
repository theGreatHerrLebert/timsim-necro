# P0 realism track — noise + spike-into-real

The gate to believable numbers. Our measured FDPs are clean (0.17–1.4%) for **two** reasons that must be
disentangled: (1) the render is **noiseless** — no chemical background to seed false matches; and (2) we
search with **DiaNN 2.5**, which is genuinely better-calibrated than the tools v1 caught. v1's FDR-inflation
headline was **tool-specific**, not universal: **Spectronaut** inflated true FDR on *modified* peptidoforms,
and **DiaNN 1.8** showed the large discrepancies — *neither has been tested here*. So a clean DiaNN-2.5 FDP
is not automatically wrong; it may be a real "modern tool is well-calibrated" result.

That means reproducing v1's finding needs **both** pieces, and they map onto the golden gate's two axes:
- **noise** (this track) → a realistic background, so an FDR claim is stress-tested rather than trivial;
- the **tool axis** (already built) → run DiaNN **1.8**/1.9/2.0/2.5 + **Spectronaut** on the *same* noised
  data and show *which* inflate. (Neither DiaNN 1.8 nor Spectronaut is installed — a dependency to add.)

Noise is still P0 regardless: even a well-calibrated engine can't be fairly stress-tested on a clean
background, and every sample-type benchmark below needs it.

**v1's two noise axes (don't conflate them).** v1 distinguishes noise on the *simulated signal* from noise
that is *background* (peaks tied to no real ion). The real v1 DIA recipe (`IT-DIA-HYE-B.toml`) turns on one
of each — signal-m/z **and** real-data-background — together:

| v1 config | axis | our track | v1 DIA default |
|---|---|---|---|
| `mz_noise_precursor/fragment` + `_ppm`, `mz_noise_uniform` | **signal** — m/z position scatter | **A1** ✅ | ON: Gaussian, 6.5 ppm |
| `detection_noise` (isotopes.py `add_detection_noise`) | **signal** — intensity shot noise | — (skip) | OFF (isotope-gen only) |
| `add_real_data_noise` (+ `reference_noise_intensity_max`, `num_*_noise_frames`) | **background** — real peaks from the ref `.d` | **A2** ✅ | ON |
| `baseline_shot_noise` (noise.py) | **background** — synthetic baseline | **A3** | unwired / legacy |
| `noise_frame/scan_abundance` | abundance | dropped | OFF |

So true v1 parity is **A1 + A2 on at once**; the intensity/abundance layers are v1-default-OFF (which is why
the plan skips them). This track has two capabilities, from "model the background" to "use a real one":

- **A. Noise** — synthesise a realistic background onto a fully-synthetic render (everything stays ground
  truth; noise just makes FDP/scoring real).
- **B. Spike-into-real** (`superimpose_on_reference` in v1) — overlay synthetic ground-truth precursors
  onto a **real experimental `.d`/`.raw`**. The real run is the background (real chemical noise, real
  interference, real dynamic range); **only the spikes are labeled**. The strongest realism there is.

Both are **opt-in and seeded** (deterministic per seed, not byte-identical to the noiseless baseline — the
byte-test stays the reproducibility gate for `--noise off`).

---

## A. Noise model

Three layers, in value/effort order. v1's machinery: `jobs/add_noise_from_real_data.py`, `noise.py`,
m/z-noise in `jobs/assemble_frames.py`.

**A1. m/z-ppm scatter (easy, all instruments). ✅ DONE — checked against v1 (`timsim-cli` 7105d73).**
Per-peak m/z jitter before m/z→tof, so a search engine sees a realistic non-degenerate mass-error
distribution to calibrate against. Matched to v1's `mscore::add_mz_noise_normal/_uniform`:
- Render flags: `--noise-mz-ppm <ppm>` (precursor) + `--noise-frag-ppm <ppm>` (fragment); `0` = off.
  `--noise-mz-uniform` selects v1's uniform mode; default is Gaussian (v1 default). `--noise-seed`.
- **v1 scale match:** the ppm value is a **3σ envelope** (v1 convention): Gaussian sd = `mz·ppm/1e6/3`;
  uniform = `mz ± mz·ppm/1e6`. So `--noise-mz-ppm 6.5` reproduces the real v1 DIA config. (v1's asymmetric
  `right_drag` tailing variant not ported.)
- Seed per `(precursor_id, is_frag, peak_index)` via successive splitmix64 avalanches (identity-keyed, like
  `survival`) so adding an ion doesn't reshuffle others — reproducible + stable under `--limit`.
  **Deliberate divergence:** v1 redraws m/z noise *per scan*; we draw once per (precursor, peak) — same
  marginal distribution, coherent across the elution (v2 projects each spectrum to tof once, then deposits).
- Applied in the projection closure; the noiseless path (`--noise-mz-ppm 0`) stays byte-identical (verified:
  hash unchanged vs the frozen baseline). Unit tests pin N(0,1)/U(−1,1) shape + MS1/MS2 key independence.

**A2. Real-data-noise injection (Bruker). ✅ DONE — checked against v1 (`timsim-cli` 4b0aaf4).** Sample
**actual background peaks from the reference `.d`** and add them per frame — window-group-aware for MS2.
The realistic chemical/electronic background; the piece that (with A1) reproduces v1's real DIA recipe.
- Flags: `--noise-real-data` + `--noise-{precursor,fragment}-frames` (v1 5), `--noise-intensity-max` (v1
  cap 150000, absolute counts), `--noise-{precursor,fragment}-fraction` (v1 0.2); reuse `--noise-seed`.
- **v1 match** (`add_real_data_noise_to_frames` DIA branch + mscore `filter_ranged`/`generate_random_sample`,
  reused verbatim): per output frame sample N reference frames of the matching type (MS1, or MS2 of the same
  DIA window group — classified via the schedule, pools from reference metadata), keep peaks with intensity
  in `[1, cap]`, downsample by the fraction, add **real detector counts** on top of the scaled synthetic
  signal (== v1's `frame + noise`; injected post-scale in `dedup_and_quantise`).
- **Deliberate divergence:** seeded (v1 is `thread_rng`) on `(output_frame, sample_slot, peak)` so it's
  reproducible for the gate — distributional/logic parity, not byte parity, with v1. Grid-safe; each unique
  reference frame decoded+filtered once (cached). Off ⇒ byte-identical to the noiseless baseline (verified).
- Bruker-only (needs the reference's real data); Thermo/SCIEX get A1 (ppm) only. Known cost: the cache is
  memory-heavy on dense references (~3GB in test) — streaming per-chunk is a future optimization.

**A3. Synthetic chemical/baseline (optional).** Poisson-count baseline peaks, uniform m/z, exponential
intensities (`noise.py`). Only if A1+A2 leave the background too clean, or for no-reference renders.

**Eval — background subtraction is REQUIRED for A2 FDP (not optional).** A2's background comes from a real
reference; even a "blank" carries real low-level peptides, and the search engine legitimately identifies
some. Those IDs are **real, not false positives** against our synthetic answer key — counting them inflates
measured FDP and would make noise look like it breaks FDR when it doesn't. So FDP under A2 must **subtract
background IDs**:
1. Render a **noise-only control**: `timsim-render --dia --noise-real-data --noise-only --noise-seed S`
   (same seed as the real run) — deposits only the A2 background, no synthetic signal.
2. Search both the real run and the control with the same engine.
3. Score with `--background-report <control_report>`: the scorer removes the control's IDs
   (`background − ground_truth`, so a real synthetic hit that also shows up isn't dropped) and reports the
   subtracted count. `fdp = false_after_subtraction / kept`.
Implemented in `timsim_eval` (`score_thermo_dia` — the flow's Bruker/Thermo/SCIEX scorer — and
`v2_eval.score`). This is the same "only labeled signal counts" principle as spike-into-real (mode B).
Recall is unaffected (background IDs are never in the truth). Everything is still ground truth; noise +
subtraction just make the measured FDP believable (toward v1's 3–5%).

---

## B. Spike-into-real-experiment

Render synthetic precursors **additively onto a real `.d`/`.raw`**, not a blank template. Output = the real
run's peaks **plus** the synthetic spikes on the same grid.

**Design (the new machinery vs. the current from-scratch render):**
- The current DIA render *writes a fresh `.d`* with only synthetic peaks. Spike mode instead **copies the
  real `.d`'s frames** and **adds** synthetic deposits per frame (matched to the real `.d`'s schedule +
  calibration, which we already read). Bruker: add onto real frames; Thermo: merge onto the `.raw`
  template's real peaks at `--spike-merge-ppm` (v1's `superimpose_merge_ppm`).
- Flag: `--spike-into <real.d>` (implies the reference geometry comes from that same real run).

**Truth + eval (the key difference):** only the **synthetic spikes** are labeled. The answer key lists the
spikes; the real-background IDs are **unlabeled** — they must NOT count as false positives. So `timsim_eval`
needs a **spike-recovery mode**: recall over the injected spikes, and background IDs excluded from FDP (or
reported separately as "background/ambient"). This is the PlasmaBENCH PYE pattern (spike human/yeast/E.coli
into real plasma, recover the known log2 ratios).

**Why it's the strongest mode:** the background is *real* — real co-elution, real noise, real dynamic
range — so a tool's spike recovery + quant accuracy is measured in a genuine matrix, with perfect ground
truth for the spikes. It's the most honest benchmark we can produce, and it needs no noise *model* (the
real run IS the noise).

---

## Sequencing + validation

1. **A1 (m/z-ppm) ✅ DONE + v1-matched.** Validated: `--noise-mz-ppm 0` byte-identical to the frozen
   baseline; nonzero shifts the mass-error distribution DiaNN calibrates against (Gaussian sd `mz·ppm/3e6`,
   reproducing v1). A1 alone won't reproduce v1's headline — that needs **A2 (below)** for the real
   background **and** the **tool axis** (Spectronaut + DiaNN 1.8 on the noised data; DiaNN 2.5 alone may
   stay clean, itself a finding). **Wired through the flow** (`--noise-mz-ppm/-frag-ppm/-mz-uniform/-seed`
   on the pipeline CLI → `render_noise_flags(cfg)` → the Bruker DIA `render` node; off ⇒ byte-identical
   command ⇒ caches unaffected). Bruker DIA only for now — Thermo/SCIEX/DDA render bins need the same noise
   closure to gain it.
2. **A2 (real-data noise) ✅ DONE + v1-matched + flow-wired.** Reuses v1's own mscore primitives
   (`filter_ranged`/`generate_random_sample`); off ⇒ byte-identical, on ⇒ deterministic. A1 + A2 together
   now reproduce v1's real DIA recipe. **FDP background-subtraction is implemented** (scorer
   `--background-report` + `--noise-only` control render — see "Eval" above). **Remaining:** (a) wire the
   noise-only control render + its search + `--background-report` into the necroflow DAG so the searched
   FDP validation runs end-to-end; (b) then confirm measured FDP moves toward v1's 3–5% on the golden
   gate's real fixture (needs DiaNN). Optional: memory-stream the A2 cache for full-scale runs.
3. **B (spike-into-real) ✅ DONE + flow-wired.** `timsim-render --spike-into <real.d>` copies every real
   frame + adds synthetic on top (v1 `superimpose_reference_frames`); reuses A2's deposit path. Control =
   `--spike-into X --noise-only` (re-encoded copy of X), searched, IDs subtracted via `--background-report`
   (spike-recovery == the A2 subtraction). Validated: control reproduces X EXACTLY (per-frame + total
   70,876,746==70,876,746); spike-full = real+synthetic; off byte-identical. A critical emit_all fix (every
   frame visited when background present) also corrected A2 under-depositing background only near signal.
   REMAINING (empirical): recover a known spike at a believable FDP in the real matrix (needs DiaNN).

Each lands behind the golden gate; seed makes every noise/spike render reproducible.

## Scope discipline (what we're NOT doing)

- No `abundance noise` (`noise_frame/scan_abundance`) as a separate stochastic layer — A1+A2 cover the
  realism a search engine cares about.
- No provenance/preview-video coupling.
- Spike mode reuses the existing render + a copy-real-frames path; it is NOT a second renderer.
