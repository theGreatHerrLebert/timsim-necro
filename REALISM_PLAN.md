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
background, and every sample-type benchmark below needs it. This track has two capabilities, from "model the
background" to "use a real one":

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

**A1. m/z-ppm scatter (easy, all instruments).** Per-peak Gaussian ppm jitter on m/z before m/z→tof, so a
search engine sees a realistic non-degenerate mass-error distribution to calibrate against.
- Render flag: `--noise-mz-ppm <ppm>` (precursor) + `--noise-frag-ppm <ppm>` (fragment); `0` = off.
- Seed per `(precursor_id, peak_index)` (identity-keyed, like `survival`) so adding an ion doesn't
  reshuffle others — reproducible + stable under `--limit`.
- Applied in the projection closure; the noiseless path (`--noise-mz-ppm 0`) stays byte-identical.

**A2. Real-data-noise injection (moderate, Bruker).** Sample **actual background peaks from the reference
`.d`** and add them per frame — window-group-aware for MS2. This is the realistic chemical/electronic
background. We already open the reference for calibration; extend to read its real `(scan, tof, intensity)`
per frame, sample a fraction (`--noise-sample-fraction`, cap `--noise-intensity-max`), deposit.
- Bruker-only (needs the reference's real data); Thermo/SCIEX get A1 (ppm) only for now.
- The v1 mechanism (`sample_precursor_signal` / `sample_fragment_signal`) is the reference.

**A3. Synthetic chemical/baseline (optional).** Poisson-count baseline peaks, uniform m/z, exponential
intensities (`noise.py`). Only if A1+A2 leave the background too clean, or for no-reference renders.

**Eval:** unchanged. Everything is still ground truth; noise just moves the measured FDP toward v1's real
3–5%. The *point* is that our recall/FDP become believable.

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

1. **A1 (m/z-ppm)** — small, all instruments, immediately makes FDP less fake. Validate: `--noise-mz-ppm 0`
   byte-identical to today; a nonzero value shifts the mass-error distribution DiaNN calibrates against. It
   should raise measured true-FDR in the *expected* direction, but the full **reproduction of v1's headline
   needs the tool axis too** — Spectronaut + DiaNN 1.8 on the noised data (that's where the inflation
   lived; DiaNN 2.5 alone may stay clean, which is itself a finding worth reporting).
2. **A2 (real-data noise)** — Bruker; validate the reference's sampled peaks land in-frame and FDP moves
   toward v1's 3–5%.
3. **B (spike-into-real)** — the additive-onto-real render + the spike-recovery eval mode. Validate: the
   real `.d`'s peaks are preserved and a known spike is recoverable in the real matrix.

Each lands behind the golden gate; seed makes every noise/spike render reproducible.

## Scope discipline (what we're NOT doing)

- No `abundance noise` (`noise_frame/scan_abundance`) as a separate stochastic layer — A1+A2 cover the
  realism a search engine cares about.
- No provenance/preview-video coupling.
- Spike mode reuses the existing render + a copy-real-frames path; it is NOT a second renderer.
