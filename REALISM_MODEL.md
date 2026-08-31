# timsim v2 realism: what we measured, what we changed, and the model it points to

Working document, 2026-08-28 → 08-31. Everything here is measured unless flagged otherwise; where a
claim was later falsified, both the claim and the correction are kept, because several were mine and
quoting only the corrections would hide how they were found.

---

## 0. Summary — and what it does NOT claim

**Established.** The protein abundance marginal (`hockey_stick`) contained a deterministic pathology:
`exp(-r/decay) + tail` puts **45% of proteins at one identical abundance** by construction (§4.1).
That is a property of the function, not an inference from any comparator, and it demonstrably
destabilises the intensity quantile ratio with detection depth. Replacing it with a log-normal
removes that pathology.

**Separately established, and a different kind of finding.** The ionisation term (`flyability`) ran
1.47× wider than the v1 it claimed parity with (§5.2). That is a **code/documentation defect**, now
corrected. It is NOT a demonstration that the biological distribution was wrong — we in fact ship
v2's wider value because it fits the comparator better (§5.3).

**NOT established, and previously overstated here.** That these were "the dominant realism defects",
or that the simulator is now "within ~40% of real". The end-to-end evidence is ONE timsTOF cohort
against ONE Bruker run, at mismatched load, AFTER tuning three coupled parameters. Critically, the
simulated cohort reports **180,933 precursors against real 49,867** (§5.4): matched top-N censoring
equalises the *reported-set* depth, not the generative sampling regime, so an apparent improvement in
conditional intensity summaries can be produced by simply making more species detectable. **The
current numbers are a tuned fit, not a validated one** — see §11 for the held-out protocol that would
change that.

The residual work is a joint calibration with uncertainty, which the parameter structure suggests
should be a hierarchical Bayesian model (§8) rather than the hand grid-search that got us here.

---

## 1. What was broken at the start

The v2 pipeline **could not render at all**. `timsim-rt` / `-ccs` / `-fragments` in the working venv
were stale shims pointing at `imspy_simulation.timsim.jobs.*`, a module never tracked in git.

**Fix:** `timsim-predict` is the intended provider — its `[project.scripts]` match the flow's command
templates verbatim, so this was a restore, not a predictor swap. Required building `mscorepy`
(pyo3), then `pepdl` and `timsim-predict` from local sources. `chronologer` and `koinapy` were
already present.

**Also broken and fixed along the way**
- `timsim-cli/target/` was deleted mid-session (freeing 48 GB); all 14 binaries rebuilt with
  `--features tdf,thermo,sciex`.
- The Thermo render path had never received the 2026-08-13 floor work (§3).
- SCIEX was not wired into the PhantomBENCH v2 driver at all (§6).

---

## 2. THE METHOD LESSON, because it caused four wrong conclusions

**Every wrong intermediate conclusion this session came from comparing two populations that had not
been through the same filter.** In order:

| # | wrong conclusion | cause | correction |
|---|---|---|---|
| 1 | "our render's p99/p50 is 15× real" | compared our RAW peaks against DIA-NN OUTPUT | checked: raw Astral 32.5 vs DIA-NN 32.3–38.8 — they agree, so the comparison was salvageable, but it was luck |
| 2 | "real p99/p50 is stable with depth" | compared DIFFERENT RUNS at different depths (confounds depth with instrument and load) | within-run top-N subsampling: real grows 2.5–2.9× over 20× depth, not flat |
| 3 | "intensity is sub-linear in load (load^0.9)" | compared the MEDIAN OF EACH LOAD's own detected set | matched precursors across both loads: ×352–400 for ×200 load — **linear, if anything steeper** |
| 4 | "the yield model alone exceeds real total spread, so the fix must go beyond abundance" | compared our UNCENSORED peptides against real DETECTED precursors | censored to matched depth: yield-only floor is sd 0.118, far BELOW real |

**Rule:** censor both sides to the same depth before comparing, and state the denominator. The
acceptance function should take a depth argument and refuse mismatched comparisons rather than
relying on the operator to remember.

**BUT THAT RULE IS TOO NARROW, and its narrowness has already bitten.** Depth is one conditioning
variable among many. The same error class covers mismatched **load** (§6 knowingly retains 200 ng vs
50 ng), **gradient**, **instrument**, **acquisition mode** (§7 switches comparator from DIA to a DDA
corpus), **search library / FASTA**, **FDR procedure**, **protein inference**, and **feature
definition**. And matched-depth censoring is itself a conditioning choice that introduces its own
bias: top-N conditions ON INTENSITY, which can erase real differences in detection probability and in
the population just below threshold (§5.4).

**Corollary:** acceptance must report matched-depth conditional summaries AND unconditional
quantities together — IDs, precursors/PG, detection curves, intensity-vs-load, near-threshold
distributions. A conditional summary alone can be satisfied by a simulator that is wrong in exactly
the way ours might still be.

---

## 3. Thermo render calibration (implemented, reviewed, not yet used for a cohort)

The floor work fixed only the Bruker path. `render_thermo` defaulted `--min-peak-intensity` to 1.0
with no inherit path, and carried `--intensity-scale 5e5` from a path whose axis is **ion counts** —
Orbitrap intensity is an arbitrary detector unit. Measured against a real Exploris template, that put
the rendered MS1 median ~12,000× below the template's median and ~25,000× below its reporting floor.

**Implemented**
- Floor inherited from the template **per MS level**. The levels genuinely differ: a stock Exploris
  DIA run censors MS1 at 25,760 and MS2 at 575.5. One shared floor censored fragments ~45× too hard,
  and a DIA search scores on fragments.
- Scale solved so the median of peaks that **survive the floor** matches the template's median — the
  two medians must be over the same population, since the template's is already above its own floor.
- `--calibrate-only` prints the calibration record and exits, so a cohort shares one frozen constant.
  Re-estimating per arm lets sample composition be absorbed into a compensating rescale, which in a
  cohort with a planted differential partially cancels the signal being measured.
- `template_level_stats` picks ONE domain per level by majority — mixing a profile baseline with
  centroid peak heights in one pool compares quantities ~4 orders apart.

**Result** (real Exploris template, matched): floor, p10 and p50 match on both levels. Residual is
upper-tail only — an abundance defect, not a renderer one (§4).

**The constant is a property of the (template, design, digest depth) TRIPLE**, not of the template
alone. A constant measured on a 3,000-peptide smoke run does not transfer to a full-digest cohort.

**Codex review** found 7 real issues, all fixed: unvalidated new args; a post-loop `lo == vals.len()`
clamp that would report a one-element "survivor distribution" while authoring nothing; the iteration
cap reported as success; the domain mixing above; SCIEX placed/tagged as Thermo; `lfq_axes.py`
duplicate reports and `--label` leakage. Its **highest-severity finding was wrong** — it claimed
`--calibrate-only` destroys `--out`; `ThermoRawWriter::from_template` only stores `out_path`, writes
are buffered to `finalize()`, verified in source and empirically with a sentinel file.

---

## 4. The abundance marginal — the largest defect we have identified

### 4.1 Diagnosis

`hockey_stick(rank) = exp(-rank/decay) + tail`, `decay=0.06, tail=1e-4`: an exponential head glued to
a **constant additive floor**. The exponential drops below the floor at rank 0.553, so **45% of
proteins sit on a flat shelf at one identical abundance**. `p50` lands on the shelf as soon as the
bottom half is included while `p99` stays in the exponential head, so the ratio explodes with depth.

Censored to the ~50k precursors a real run detects: **p99/p50 = 6688, depth-swing 743×**, against real
Bruker Ultra 2's **24.1** and **2.5×**.

The shelf also distorts the **recall-vs-load** curve — raising load pushes the whole shelf across the
detection threshold at once instead of gradually — which is the curve ramp-005 exists to measure.

### 4.2 What real abundance looks like

Real protein rank-abundance (DIA-NN `PG.MaxLFQ`, Astral 250 pg, 4,419 PGs), normalised to the top
protein, against the hockey stick:

| rank | real | hockey_stick |
|---|---|---|
| 0.01 | 0.0878 | 0.846 (10× too flat at the head) |
| 0.50 | 0.00124 | 0.000340 (3.6× too low mid) |
| 0.99 | 6.42e-05 | 0.000100 (real still decaying, ours flat) |

Real dynamic range 5.0 orders vs 4.0. Form is **approximately log-normal**: log10 mean 3.988, sd
0.665, skew +0.486, excess kurtosis +0.409, Shapiro-Wilk W = 0.986. Log-normal in FAMILY but not
exactly — more overall spread with a less extreme top than a log-normal at matching spread.

### 4.3 The fix, and the regression guard

`AbundanceProfile::LogNormal` already existed. Switched the design configs to it.

Added `the_abundance_marginal_keeps_its_quantile_ratio_stable_with_depth` (`timsim-chem/src/design.rs`),
bound 8×. **Verified it has teeth: passes on LogNormal, fails on HockeyStick at 1649×** — which
independently reproduces the numpy estimate of 1655× from the Rust implementation. The existing
`the_hockey_stick_curve_gives_a_realistic_dynamic_range` checks the RANGE and cannot see the SHAPE;
it is now annotated to say so.

### 4.4 End-to-end verification (σ=1.5, n_proteins=2400)

3v3 Bruker cohort re-rendered and searched with **byte-identical DIA-NN flags** to the hockeystick
cohort. Scored at matched depth (top-46,849 — the ID count itself changed, so "100%" is not
comparable without this):

| | 5% | 25% | 50% | 100% | swing | log10 sd |
|---|---|---|---|---|---|---|
| hockeystick | 15.7 | 55.6 | 116.8 | 365.7 | 23.3× | 0.941 |
| **lognormal σ=1.5** | 14.5 | 25.9 | 38.8 | **72.9** | **5.0×** | **0.582** |
| REAL Bruker Ultra 2 | 9.4 | 14.0 | 16.6 | 23.2 | 2.5× | 0.483 |

Defect reduced from ~300× (design) / ~15× (searched) to ~3×. Raw-domain cross-check: same TIC (+2%,
mass balance held) with **45% more peaks above the floor** — the shelf's peaks becoming detectable.

---

## 5. `n_proteins` and `flyability` — the second and third terms

### 5.1 `n_proteins`: coverage depth

Searched cohort gave **34.0 precursors per protein group** against real 4.4–11.9 (Bruker Ultra 2 5.4,
Astral 15 min 11.9, Astral 250 pg 4.4, ZenoTOF 6.5). Cause: 200 ng landing in ~12% of the proteome,
so every present protein is covered ~4× more deeply than reality.

Raised 2400 → **9000** (real Bruker reports 9,179 PGs). It is **not independent of σ**: spreading the
same load over ~4× more proteins divides the per-protein amount and changes what clears the floor.
Calibrating σ first and then changing `n_proteins` invalidates σ.

### 5.2 `flyability`: THE DOMINANT peptide-level term, and it was not v1-parity

`ionization_propensity` (drawn in `timsim-precursors`) has **log10 sd 0.880** at the binary default
`sigma = 1.0`, measured over 9,010,877 precursors. That is **larger than the abundance contribution
(0.245–0.412) and larger than real data's ENTIRE observed precursor spread (0.483)**.

`ionize.rs` claims "v1's parameters (median 1e-2, σ = 1)". **v1's `generate_normal_efficiency` uses
`std_log = 0.6`** in log10 space. `flyability_sigma` is in log10 units — measured: 0.6 → sd 0.597,
1.0 → sd 0.880. **v2 ran 1.47× wider than the v1 it claims to match.** The median (1e-2) does match.

**Every σ calibration before this measured `peptide_quantities`, which is UPSTREAM of flyability** —
tuning the smaller term while blind to the larger. Part of the "1.41× render+search overhead" was
just flyability applied at the precursors node.

**THE QUANTITY THE RENDERER CONSUMES** (`render.rs:996`):

    abundance = amount_amol × ionization_propensity × modform_fraction × charge_fraction

Calibrate against **that product**, never `amount_amol` alone.

### 5.3 Joint calibration on the product (n=9000, top-50k, design level)

| σ | flyability | sd | p99/p50 | swing |
|---|---|---|---|---|
| 1.10 | 0.6 | 0.341 | 15.5 | 2.1× |
| 1.40 | 0.6 | 0.371 | 20.2 | 2.3× |
| **1.60** | **0.6** | 0.392 | **24.0** | **2.3×** |
| 1.60 | 1.0 (old) | 0.424 | 31.5 | 3.9× |
| REAL Bruker Ultra 2 | – | 0.483 | 23.2 | 2.5× |

**SHIPPED SET: `n_proteins=9000, abundance σ=1.6, flyability_sigma=1.0`** — chosen on EVIDENCE, not
on v1 parity. See §5.4: the end-to-end cohort at σ_fly=1.0 gives searched log10 sd **0.437** against
real **0.483**, i.e. already slightly NARROW, so dropping to v1's 0.6 moves further from the data.
v1 was itself a model; parity with it is not a claim about reality. The misleading "v1's parameters"
doc comment in `ionize.rs` is corrected to state the real value and the measurement.

**σ_fly is NOT separately identifiable from the abundance σ** — intensity observes only their
product. The pair belongs in a joint fit (§8), and the grid search below found *a* good point on a
ridge, not *the* optimum.

### 5.4 END-TO-END RESULT (2026-08-31) — the shipped set, searched

3v3 timsTOF cohort at `n_proteins=9000, σ=1.6, σ_fly=1.0`, searched with DIA-NN flags identical to
both earlier cohorts, scored at matched depth (top-46,849, the smallest set):

| | n total | 5% | 25% | 50% | 100% | swing | log10 sd |
|---|---|---|---|---|---|---|---|
| hockeystick, n=2400 | 46,849 | 15.7 | 55.6 | 116.8 | 365.7 | 23.3× | 0.941 |
| lognormal σ=1.5, n=2400 | 91,057 | 14.5 | 25.9 | 38.8 | 72.9 | 5.0× | 0.582 |
| **lognormal σ=1.6, n=9000** | 180,933 | **9.8** | **17.1** | **22.4** | **32.4** | **3.3×** | **0.437** |
| REAL Bruker Ultra 2 | 49,867 | 9.4 | 14.0 | 16.6 | 23.2 | 2.5× | 0.483 |

- p99/p50 at full depth **365.7 → 72.9 → 32.4** (real 23.2): from 15.8× off to **1.4×**
- depth-swing **23.3× → 5.0× → 3.3×** (real 2.5×): from 9.3× off to **1.3×**
- log10 sd **0.941 → 0.582 → 0.437** (real 0.483): from 1.95× too wide to **0.90×, now slightly narrow**
- the 5% quantile is within 4% of real (9.8 vs 9.4)
- **8,639 protein groups** at 1% global FDR, against a 9,000-protein design and real's 9,179

**INTERNAL CONSISTENCY, NOT VALIDATION.** The design-level prediction for this config was sd 0.424,
p99/p50 31.5, swing 3.9×; measured searched 0.437, 32.4, 3.3× — within ~10%. This was previously
written up as "method validation". **It is not.** Both numbers are generated by the same fitted model
and evaluated with closely related statistics, so agreement is expected and carries no external
evidential weight. What it does show is narrower and still useful: measuring the FULL PRODUCT at the
precursors node (rather than `amount_amol` alone) makes the design→searched mapping stable, which is
why the earlier "overhead multiplier" was erratic.

**THE ID-COUNT CONFOUND — the main threat to everything in this table.** We report 180,933 precursors
against real 49,867. Matched top-N censoring conditions ON INTENSITY: it equalises how many reported
values enter the comparison, NOT the generative regime that produced them. A simulator with
implausibly many detectable species can therefore look *better* on these conditional summaries. The
confound is not confined to this row — the §4.4 switch also doubled IDs (46,849 → 91,057) at
identical `n_proteins`, so it is present in the "controlled" comparison too. Nothing in this document
separates "the shape is now right" from "more species are now detectable".

**Caveat on absolute counts:** this cohort is 200 ng against a 50 ng comparator (§6), and identifies
180,933 precursors vs real 49,867. Load-adjusted (precursors ∝ load^0.37) real at 200 ng would be
~84k, so we are still ~2× deep. A 50 ng arm removes that mismatch.

`flyability` / `flyability_median` / `flyability_sigma` are now settable from a job TOML (they were
on a binary default) — plumbed through all six pipelines and the CLI in `timsim_flow.py`.

---

## 6. Load → intensity

**Linear over the range tested.** Following the SAME precursors across a 200× load change gives
×352–400 signal (exponents 1.107–1.131 on the best-matched pairs). The apparent sub-linearity
(×78–131) was a censoring artifact of comparing medians of differently-censored sets.

**Do not add a compression exponent on this evidence — but do not claim linearity is established
either.** The result rests on TWO load levels and a selected matched-precursor set (precursors
detected at BOTH loads, which biases toward the bright end). It rules out neither saturation nor
instrument-dependent compression outside 250 pg – 50 ng, and the paired runs differ in date and
column, so the >1 exponents should not be read as genuine super-linearity.

Depth scaling: real ×6.8–7.4 precursors for ×200 load (exponent ~0.37); ours ×4.8 (0.30). Close, not
a priority.

Only two load levels exist in the LFQ set (250 pg, 50 ng), so an "exponent" is one ratio restated,
not a fit.

**Open denominator issue:** our cohort is 200 ng and the real Bruker comparator is 50 ng. Load-adjusted,
real at 200 ng would sit near ~9 precursors/PG rather than 5.4. **Render a 50 ng arm** to remove the
last mismatch rather than tuning against it.

---

## 7. RT and IM

| | ours | real Bruker Ultra 2 | verdict |
|---|---|---|---|
| FWHM | 4.31 s | 1.93 s @5 min | **ok** — real FWHM scales ~`gradient^0.39` (3× gradient → 1.53×), so 40 min predicts ~3.67 s; ours 1.17× wide |
| RT residual sd | 0.98 s | 0.52 s | 1.9× wide, but IQR only 1.24× → heavy **tails**, core ok |
| IM residual sd | 0.00135 | 0.00317 | ours narrower — but see below, this is largely by construction |
| IM p5–p95 | 0.688–1.236 | 0.858–1.209 | ours extends too low |

**Do NOT gradient-scale FWHM linearly.** It goes as ~`gradient^0.39`, measured across 5/11/15 min.

**IM and RT WIDTHS are checked against a modern DDA corpus — this supersedes the old-instrument check
below. Note the scope: peak WIDTHS only. Residual LOCATION structure (systematic vs random), the skew
sign convention, and whether a DDA-derived width transfers to the DIA rendering target all remain
open, so this is not a general validation of the mobility/RT model.**

`theGreatHerrLebert/timstof-dda-pasef-cc0` (HF, CC0, "Claudius timsTOF DDA-PASEF PSM corpus",
10M-100M rows; tier1 6.45 GB / tier3 10.03 GB). Local: `/media/hd02/data/timstof-dda-pasef-cc0/`.
It carries **directly fitted peak shapes with quality gates** — `ms1_{im,rt}_{apex,fwhm,sigma,skew,
r2,snr}` plus `{im,rt}_width_reliable` (gates on r2>=0.8 & snr>=20) — so the widths are measured, not
inferred from PSM scatter. It also carries `ms1_iso_0..4` (isotope intensities) and a tier3 table of
per-fragment intensity by type/ordinal/charge.

Measured on the validation split (3,265,602 PSMs), **reliable fits only**:

| | n (reliable) | p25 | p50 | p75 | mean |
|---|---|---|---|---|---|
| IM sigma (1/K0) | 941,095 (28.8%) | 0.00575 | 0.00791 | 0.01218 | **0.00962** |
| RT sigma (s) | 822,663 (25.2%) | 2.47 | **3.25** | 4.65 | 3.85 |
| RT FWHM (s) | " | 3.79 | 5.50 | 7.00 | |

- **`V1_MOBILITY_STD_TARGET = 0.009` targets the population MEAN. Measured mean 0.00962 — within 7%.**
  The constant is right, now on current-generation hardware. The distribution is right-skewed
  (median 0.0079 < mean 0.0096), so targeting the mean is the correct choice.
- **`--sigma-seconds` default 3.0 vs measured median 3.25** — within 8%. Our searched cohort's FWHM
  (4.31 s) is NARROWER than this reference's 5.50 s, reversing the direction inferred from LFQ's
  5-min DIA runs. This DDA corpus is the better comparator: fitted peaks, not gradient-extrapolated.
- Only ~25-29% of fits pass the quality gate — using all PSMs would give a different, worse answer.
- **Peak skew is non-zero in BOTH dimensions**: RT median -0.437, IM median -0.956. Our renderer
  models asymmetry in RT only (EMG) and none in mobility. **The fitter's skew sign convention must be
  checked before concluding direction** — EMG produces tailing; whether -0.437 means fronting depends
  on their definition.

**Location is NOT something to extract.** A peptide's mobility is a function of its sequence — that is
what the CCS deep model is for. What IS extractable is the **location variance**: `observed -
predicted`, and critically **whether it is systematic or random**. That distinction matters: a search
engine can calibrate out systematic offsets (DIA-NN does per-run IM calibration) but not random
scatter, so simulating bias as noise would make the search artificially harder. **Not yet measured.**

### How the model works today (unchanged, and correct)

The CCS deep model predicts a per-ion `ccs_std` (learned shape), rescaled by a constant gain so the
population mean hits `V1_MOBILITY_STD_TARGET` (`render.rs:210-234`). v2 improved on v1 by using a
fixed denominator rather than a per-run mean — v1's version made an ion's width depend on which other
ions were in the run (4.1% swing across subsets).

### Superseded: the old-instrument check

`/media/hd02/data/raw/dda/ccs` (787,419 PSMs) gave within-run repeat-PSM spread rising monotonically
with sampling depth (n=2: 0.00184 -> n=5-6: 0.0084-0.0094 -> n=10+: 0.01611), consistent with PASEF
apex-bias and bracketing 0.009. **But those runs are 2019-era `TIMS1`/`TIMS2`/`tims03`**, so they do
not speak to current hardware. Superseded by the corpus above. Also: all three FragPipe searches
there are **yeast** (zero HeLa PSMs), and overlap with our human digest was only 479 (sequence,
charge) pairs — too few for residual-structure work regardless.

## 8. The hierarchical model this is pointing at

The generative chain is already a hierarchy; we have been fitting it by hand, sequentially, on
summary statistics:

    protein abundance   ~ LogNormal(μ, σ)              population level
    peptide | protein   = deterministic digestion physics (P(cut)·P(cut)·Π(1-P(cut)))
    flyability          ~ LogNormal(1e-2, σ_fly)       per-peptide random effect
    charge, modform     = multiplicative fractions
    RT_obs   | peptide  = Predicted.RT + ε_rt          (ε_rt heavy-tailed vs Gaussian)
    IM_obs   | precursor= Predicted.IM + ε_im          (systematic vs random UNKNOWN — §7)
    intensity| precursor= product of the above, then CENSORED at the instrument floor

**For:**
- The parameters are coupled (σ ↔ n_proteins ↔ flyability) and we discovered that the hard way, one
  at a time, with wrong conclusions in between. A joint fit gets the covariance.
- Censoring is the central difficulty and is exactly what a likelihood handles cleanly. All four
  errors in §2 were censoring/matching errors.
- Three observation channels (intensity, RT, IM) on the same latent peptide give more constraints
  than three separate fits — intensity alone cannot separate abundance from flyability because it
  observes only their **product**, but RT and IM are independent channels on the same latent.
- Instrument is a natural grouping level (Bruker / Astral / ZenoTOF share peptide latents, differ in
  dispersion) — a partial-pooling problem.
- We would get uncertainty, not point estimates. We are shipping these into a benchmark and cannot
  currently answer "σ = 1.6 ± what?".

**Against / limits:**
- **Identifiability.** abundance × flyability is a product and intensity observes the product; from a
  single run they are identifiable only up to it. What breaks it: flyability is a per-peptide property
  shared across runs, abundance varies by condition. Repeated measures plus the planted design
  separate them — but only if the model is specified to exploit that.
- **Circularity.** Fitting the simulator to LFQ data and then benchmarking tools on the simulator puts
  some of the data inside the answer. Keep the fitting data and the benchmarking claim distinct.
- **The forward model is not differentiable and costs hours** (Rust renderer + DIA-NN). Inference must
  be on the **generative sub-model** fitted to DIA-NN outputs, not the full pipeline.

**Scope proposal:** protein → precursor abundance with an explicit censoring likelihood, fitted
against real DIA-NN reports, renderer and search kept out of the inference loop.

---

## 9. Infrastructure defects found

**necroflow re-runs a STALE node IN PLACE.** The node directory keeps the same hash and path, and the
old artifact is overwritten. Downstream fingerprints reference the parent by **path**, so in a LATER
invocation they stay `up_to_date` against data whose content changed underneath them.

This silently produced a **mixed cohort**: 1 arm on the new config, 5 on the old, with `necroflow`
reporting "done: 6 arms" and printing a `linked` line for every one. Only comparing file mtimes and
sizes caught it.

Same family as the documented binary-fingerprint gap, which was fixed for BINARIES by stamping
`# tools:<sha>` into the command template. **Data parents are still identified by path, not content.**

- **Detect:** `stat -c%y` every arm after any config change, before searching. Split mtimes/sizes =
  mixed cohort. Do not trust "done: N completed" or `linked` lines.
- **Workaround:** `rm -rf` the stale node dirs, or use a **fresh `--work-root`** after a config change.
- **Real fix:** stamp the parent's output content hash into the dependent's fingerprint, or make node
  dirs content-addressed on output.
- **Any prior result produced by re-running into an existing work root after a config change is
  suspect**, not just this one.

---

## 10. External data

**LFQ Benchmark Gen Beta** — PXD070049 (raw, ~2.9 TB) + Zenodo 17936657 (16 GB processed).
Orbitrap Astral, ZenoTOF 7600+/8600, timsTOF Ultra/Ultra 2; 2,340 runs; human/E. coli/yeast at
**exact published ratios**, giving an external answer key our own planted design cannot provide:

| organism | Cond A | Cond B | Cond C |
|---|---|---|---|
| E. coli | 0.05 | 0.20 | 0.32 |
| Human | 0.65 | 0.65 | 0.65 |
| Yeast | 0.30 | 0.15 | 0.03 |

→ B vs A: E. coli **+2.000**, human **0.000**, yeast **−1.000**; C vs A: +2.678 / 0 / −3.322.

**PRIDE has essentially no processed results** (2 `.mzid`, 3 `.mgf`, 1 FASTA). The search output is on
**Zenodo**, processed with AlphaDIA, DIA-NN, FragPipe, PEAKS and Spectronaut.

**DENOMINATOR WARNING:** 5/11/15/30-min gradients at 50 ng and 250 pg. Our cohorts are ~40 min at
200 ng. Match gradient and load before using as a comparator; the timsTOF Ultra arms are closest to
our Bruker work.

Local: `/media/hd02/data/lfqbench-beta/` — SDRF, HYE FASTA, the Zenodo zip, `extracted/`, and
`raw-onepertype/` (one DIA run per instrument, ~8 GB, for building analysis scripts locally before
slurm-submitting at scale on MOGON).

**Tooling:** `analysis/lfq_axes.py` scores ID density, RT width and the scale-invariant p99/p50 from
any DIA-NN report — including our own, via `--report/--label`, so our output sits on identical
footing to five real instruments.

---

## 11. What would actually validate this (the held-out protocol)

Everything in §4–§6 is a **tuned fit against one comparator**. The following is the minimum that
would turn it into a validated one. It is written down before running so the acceptance criteria
cannot drift to match whatever comes out.

**Freeze first, then test.** Freeze `n_proteins`, abundance σ, σ_fly, AND the depth-matching rule.
No re-tuning against any held-out set; a failure is a result, not a prompt to adjust.

**Held-out axes, in order of value**
1. **50 ng arm** — removes the load mismatch that currently contaminates every absolute comparison.
2. **Other Bruker batches and gradient lengths** — tests that the fit is not specific to one run.
3. **Astral and ZenoTOF** — tests that it is not specific to one vendor.

**Quantities NOT used in fitting** (the fit used p99/p50, depth-swing, log10 sd only):
- precursors per protein group, and IDs — **unconditional**, which is where the ID-count confound
  (§5.4) would show itself
- missingness structure across replicates
- between-replicate variance
- load-response curves (recall vs load)
- fragment concentration (§12 reference: median 13 fragments, top-1 at 19.6%)
- RT and IM residual structure, including whether it is systematic or random
- **HYE fold-change recovery against the published ratios** ([[§10]]: E. coli +2.000, human 0.000,
  yeast −1.000) — the one axis with EXTERNAL ground truth, and therefore the strongest single test

**Report the ridge, not a point.** σ and σ_fly are not separately identifiable (§5.3). Report the
profile of acceptable combinations with uncertainty, not one selected pair.

**Disjoint calibration and evaluation data.** Fit on one subset of the LFQ Benchmark, evaluate on
held-out instruments and runs. Without this the simulator can reproduce DIA-NN's own reporting
artifacts and then be used to judge DIA-NN-like tools (§8, circularity).

**A model that matches p99/p50 but fails reproducibility or organism-ratio recovery has not earned
the word "realistic".**

---

## 12. State and sequencing

**Done — but see §0 and §11 for what "done" does and does not mean here**
- **abundance + coverage + ionisation jointly TUNED and run end-to-end** (§5.4): conditional
  intensity summaries land within ~40% of one Bruker comparator. This is a fit, not a validation:
  three coupled parameters, one instrument, one run, and an unresolved ID-count confound.
- pipeline restored end-to-end (render → DIA-NN → score)
- abundance shape fixed and verified end-to-end; regression guard with teeth
- flyability corrected to v1 parity; exposed to job configs
- Thermo floor/scale calibration implemented, codex-reviewed, bugs fixed
- SCIEX arm wired (own link root, tag and glob hint)
- all four repos committed and pushed on branches; `main` untouched — **operational status, not
  scientific verification.** And §9's stale-parent fingerprint defect means the provenance of ANY
  result produced by re-running into an existing work root is an open audit question.

**Queued**
1. ~~**timsTOF 3v3** — DONE (§5.4). Shipped set verified end-to-end.~~
2. **50 ng arm** to remove the load denominator mismatch (§6).
3. **Thermo arm** — never rendered a cohort; needs its frozen constant against the cohort's own
   feature space.
4. **SCIEX arm** — wired, never executed once; shares the `fragments` node with Thermo.
5. **PhantomBENCH core** last — its DE/GO layer sits on all three arms.

**Open**
- **fragment allocation** — the one genuinely unexamined axis. Reference measured from the HF corpus
  (tier3, one row group, 2,774 fragments): **median 13 matched fragments per precursor** (p25 8,
  p75 17, p95 22), **top-1 carries 19.6%** of the precursor's MS2 signal, **top-3 carry 46.4%**,
  y:b ≈ 1.33, singly-charged dominant. Compare our `prospect-local` predictions against this: a model
  that concentrates signal into too few fragments makes DIA scoring easier than reality, the same way
  the abundance shelf did for precursors. A proper pass is chunked over the 49.2M-row split.
- **isotope allocation — CLOSED, not a defect.** Real envelopes agree with theory at
  `isotope_cosim` median **0.9796**, and our renderer computes envelopes from theory (mscore). Real
  medians (z=2) 0.365/0.314/0.192/0.113 for M+0..M+3, shifting rightward with charge as physics
  requires. Note the corpus extracts at most 4 isotopes (M+4 identically 0), so truncate to match
  before comparing.
- **CCS location residual — systematic vs random.** Still open and still the interesting question: a
  search engine can calibrate out a systematic offset (DIA-NN does per-run IM calibration) but not
  random scatter, so simulating bias as noise makes the search artificially harder. The HF corpus
  supersedes the planned HeLa search of the 2019 `Raw_HeLa_Trp/` data — it is a join against 3.27M
  modern PSMs instead.
- **peak skew** — real peaks are asymmetric in BOTH dimensions (RT median −0.437, IM median −0.956).
  We model asymmetry in RT only (EMG) and none in mobility. Check the fitter's sign convention before
  concluding direction.
- **the hierarchical model (§8)** — and note §5.3/§5.4: σ and σ_fly are not separately identifiable,
  so the shipped values are *a* point on a ridge. That is the concrete argument for the joint fit.

**A governance question for PhantomBENCH.** We have built a detector for our own fabrication:
p99/p50 depth-stability, log-normality and kurtosis all separate our synthetic data from real. Against
**DIA-NN alone** the fabrication still passes — the search engine compresses and masks most of the
defect, which is why raw-level inspection was necessary to see its true size. But a reviewer applying
these three statistics would catch it. That is a stronger result than the original thesis, not a
weaker one, and it is produced by the honest sibling work rather than by the fabrication — but it is
a different paper from the one the repo currently describes.
