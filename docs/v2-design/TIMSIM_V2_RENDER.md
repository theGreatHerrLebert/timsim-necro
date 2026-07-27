# timsim v2 — the streaming frame render

Companion to `TIMSIM_V2_SPEC.md` / `TIMSIM_V2_PLAN.md`. This is the design for the **render**: the
last, heaviest measurement stage, which turns the instrument-independent feature space into a real
Bruker `.d`. It is the S4 chunk of the strangler plan, and the single largest remaining piece.

**Goal:** render a full run in **bounded memory (target: peak RSS < 1 GB)** and faster than v1, on a
6-year-old 32 GB machine — so that render cost scales with the *elution window*, not the *run length*.

---

## 1. The problem with v1's render

v1's frame assembly (`assemble_frames.py` → `rustdf/src/sim/lazy_builder.rs`) already reaches for the
right idea — `load_data_for_frame_range(frame_min, frame_max)` loads only peptides eluting in a frame
window. But the *how* is expensive:

- **SQLite-backed.** `read_peptides_for_frame_range_with_source`, `read_ions_for_peptides_with_source`
  and friends are string-column queries against the `.d`'s SQLite, re-issued every batch. Slow, and
  it is the same by-column coupling the schema work exists to kill.
- **Batch-based, not a true sweep.** It loads a *block* of frames at once, so RAM scales with batch
  size, and a peptide whose elution spans a batch boundary is loaded twice (or more).
- **Materialization.** Per batch it builds `frame_to_abundances` maps and holds the batch's peptides,
  ions, and fragment vectors together.

The result is RAM that grows with how much you load at once, and throughput bounded by SQLite.

---

## 2. The insight: the data is a 1-D temporal sweep

A peptide affects only a **contiguous run of frames** — its elution window (~5–10 s) is tiny against
the gradient (~30–120 min). Its m/z, isotope envelope, CCS, and fragment intensities are *fixed* over
that window (they are feature-space facts, not per-frame). So at any retention time, only a small
**active set** of peptides contributes, and the whole render is a sweep-line over RT.

This is the "rolling pattern": process frames in RT order, keep an active set that peptides enter and
leave exactly once, accumulate each frame's peaks, emit, and drop.

---

## 3. Design: the streaming sweep-line render

### 3.1 Inputs (the Parquet feature space, streamed — no SQLite)

| artifact | provides |
|---|---|
| `precursors` | precursor_id, peptide_id, modform, charge, m/z, isotope_intensity (envelope), charge_fraction |
| `precursor_ccs` | mobility apex (CCS → 1/K₀ at the run's gas) + width (ccs_std) |
| `peptide_rt` | rt_index (apex) + rt_sigma_hat, rt_k_hat (the elution **peak shape**, normalized) |
| `fragment_intensities` | the MS2 payload: per-precursor (ion_type, ordinal, frag_charge, intensity) |
| `peptide_quantities` | abundance per peptide, per sample |

### 3.2 Precompute (per precursor, once)

Given the run's instrument config (gradient length + shape, drift gas/temp, isolation scheme,
frame cycle, scan↔mobility and tof↔m/z calibrations):

- **Elution window** `[frame_start, frame_end]`: map rt_index → apex seconds → apex frame; map
  sigma_hat/k_hat → absolute EMG (σ, k) via the gradient band; truncate the window at a probability
  target `target_p` (as v1 does) so the tail is bounded.
- **Mobility peak** `(scan_center, scan_width)`: CCS → 1/K₀ → scan_center; ccs_std → scan_width.
- **m/z → tof** for the precursor (MS1) and for each fragment (MS2).

### 3.3 Terminology: frame vs duty cycle

A timsTOF **frame** is one TIMS ramp (~100 ms), producing many **scans** (the mobility axis, ~1 per
µs of ramp). A **duty cycle** is one MS1 frame followed by *N* MS2 frames. The acquisition schedule
assigns each frame a role (MS1, or MS2 with a window group). The sweep advances over **frames**; RT
integration and window assignment are per frame. (v1's ~36k count is frames, not cycles.)

### 3.4 The sweep — per-ion → per-scan → frame, mirroring v1

The render's inner structure is **not** "one peak = one (frame, scan, tof)". v1 builds a frame by, for
each active ion, emitting a per-scan spectrum, then layering the scans into the frame
(`build_precursor_frame` → `from_tims_spectra_filtered`). The streaming render keeps that structure
and only changes *what feeds it* (Parquet, swept) and *when it writes* (incrementally). Per frame `f`:

1. **Advance the active set** (a min-heap keyed on `frame_end`): push precursors with
   `frame_start ≤ f`, pop those with `frame_end < f`.
2. **Role from the schedule:** MS1, or MS2 with a window group `g`.
3. **For each active precursor, for each scan it occupies** (the scans its mobility peak covers):
   - the scan's weight is `mobility_weight(scan)`; the temporal weight is `elution_weight` integrated
     over *this event's exposure interval* `[start, end]` (§6), not a frame index.
   - **MS1:** emit the precursor's isotope envelope into the scan spectrum, scaled by
     `abundance · charge_fraction · elution_weight(exposure) · mobility_weight(scan)`.
   - **MS2:** emit its **transmitted** fragments, weighted by `elution_weight · mobility_weight` and by
     transmission — where transmission is applied by **exactly one** of v1's modes (§3.5), never two.
4. **Layer** the per-scan spectra into the frame; add noise / superimpose the reference frame (§6);
   emit the frame to TDF (§5); drop it.

### 3.5 The transmission is NOT `m/z ∈ W` — it is the diagonal, and v1 has it

The naive "MS2 if m/z is in the isolation window" is wrong for diaPASEF: **the quadrupole moves during
the TIMS ramp**, so the isolation window is a function of *scan* (mobility), and its edges are a soft
roll-off, not a hard cut. v1 already models this precisely: `IonTransmission::apply_transmission(frame_id,
scan_id, mz)` and `TimsTransmissionDIA::get_setting(window_group, scan_id)` give a **per-scan window
center/width** (the diagonal) with a `k`-parameterized edge, and `is_transmitted(frame, scan, mz)`.
The streaming render **reuses this verbatim**; it does not reinvent transmission. So the factorization is:

```
MS1:  I(frame, scan, tof) = abundance · charge_fraction
                          · elution_weight(exposure)     # emg integral over the EVENT's exposure interval
                          · mobility_weight(scan)        # cdf_normal range   (mscore::utility)
                          · isotope_envelope(tof)        # normalized; no quad on MS1

MS2:  I(frame, scan, tof) = abundance · charge_fraction
                          · elution_weight(exposure) · mobility_weight(scan)
                          · TRANSMISSION(mode; mz, scan, frame)   # ONE of the modes below
                          · fragment_intensity(tof)               # relative, base-peak normalized
```

**Transmission is one mode, never two — the equation must not double-count `Q`.** v1 offers three
*mutually exclusive* modes, and the render picks one:

- **`PrecursorScaling`** — scale the fragment payload by a single precursor transmission factor
  `Q(mz, scan, frame)` (`apply_transmission`, the soft diagonal).
- **`PerFragment`** — derive the fragment isotope distribution from the *set of transmitted precursor
  isotopes* (`calculate_transmission_dependent_fragment_ion_isotope_distribution`). This already folds
  transmission in, so it must **not** also be multiplied by `Q`.
- **`None`** — gate only (a fragment is present iff its precursor is transmitted at all), no scaling.

Multiplying a `PerFragment` payload by `Q` applies transmission twice; the earlier draft did exactly
that. The render selects the mode per the instrument capability and applies transmission once.

**Soft vs threshold, stated honestly.** `apply_transmission` is a *soft* curve (the `k`-edge), but the
`PerFragment` helper takes a **thresholded** `HashSet` of transmitted isotope indices, not soft
probabilities. So `PerFragment` fragment chemistry is a threshold approximation of the soft edge. If
soft-edge fidelity is needed for MS2, add a probability-weighted isotope-distribution variant;
otherwise this approximation is documented, not hidden.

**Collision energy** is baked into the fragment payload (intensities predicted at the window group's
CE — how `fragment_intensities` is produced), so the render selects the fragment set matching the
event's CE rather than carrying a separate CE factor. Detector **saturation / dynamic range** is a
known simplification, added as a post-accumulation per-scan clip (§9).

**Two separate conservation checks — MS1 and MS2 are different acquisitions, do not sum them.**

- *MS1 oracle (exact):* `charge_fraction` splits the peptide amount across charge states *once* before
  isotope emission, so for a fixed charge, `Σ` over `(frame, scan, tof)` of the emitted **MS1**
  intensity equals `abundance · charge_fraction` up to the *measured* truncation (elution `target_p`,
  mobility range, isotope depth). The elution, mobility, and isotope factors are exact CDF integrals
  over a partition of each axis, so this is exact and is the render's primary oracle.
- *MS2 is not a partition of abundance.* Fragment intensities are **relative** (base-peak normalized),
  and MS2 is a separate acquisition from MS1 — summing MS1 and MS2 does **not** return the peptide
  amount, and claiming so would be wrong. MS2 is validated separately: against the chosen transmission
  mode (a `PerFragment` spectrum must match the transmitted-isotope model; a `PrecursorScaling` one
  must scale linearly with `Q`) and against the golden micro-fixtures (§8), not against a conservation
  sum.

---

### 3.6 The sweep requires frame_start-sorted input — a sorted artifact, not a sort

The active-set bound holds only if precursors **arrive in nondecreasing `frame_start`**. Parquet is
not RT-sorted, and materializing the whole precursor set to sort it defeats the purpose. So a
`timsim-precompute` step (§3.2) emits a **`render_precursors` artifact sorted by `frame_start`** —
carrying, per precursor, the elution window, mobility peak, tof, and a pointer to its fragment set —
produced with an external/partitioned sort so the sort itself is bounded. The sweep is then pure
sequential I/O + accumulation over that artifact. (Fragment intensities can stay in their own table,
read by the sweep for active precursors only, or be denormalised into the sorted artifact — §9.)

## 4. Memory analysis (why 32 GB is headroom — with the real pressures named)

For a 60-min gradient, ~36k frames (100 ms ramp), ~100k peptides, ~7 s elution width:

- **Active set** ≈ `(7 s / 3600 s) · 100k` ≈ **~200 precursors** eluting at once; 10× margin → 2000,
  ~100k active fragment-ions → **megabytes**. Bounded by **elution width, not run length** — a 3-hour
  gradient costs the same as a 30-minute one.

But the active heap is **not the only pressure**, and the earlier "a few thousand sparse bins" was
optimistic. A dense diaPASEF frame co-isolates many precursors × fragments × mobility bins; the frame
buffer, hash-map overhead, duplicate `(scan, tof)` aggregation, the reference-frame data (§6), and the
parallel-compression queue can each dominate. The design bounds **all** of them:

- **Bounded producer/compressor channels** — a fixed-depth queue between the sweep and the (parallel)
  compress-and-write stage, so a slow writer applies backpressure instead of growing an unbounded
  buffer.
- **Per-frame peak cap + overflow policy** — a hard ceiling on `(scan, tof)` entries per frame with a
  defined, deterministic policy on overflow (drop lowest-intensity, logged), never silent OOM.
- **Peak-RSS telemetry** — the prototype reports peak RSS; the target is **< 1 GB** on a realistic run.
- **A stress fixture** — concentrated RT (everything co-eluting), broad elution tails, and a
  narrow-gradient/enrichment case, to prove the bound under pathology, not just the happy path.
- **Cap/spill (if ever needed) must preserve deterministic reduction order** and must not silently
  alter intensities — determinism is a correctness property here, not a nicety.

---

## 5. The incremental TDF writer — append-only is confirmed feasible

This was flagged as the load-bearing risk; the code review *retired* it. A `.d` is `analysis.tdf`
(SQLite metadata: frames, scans, calibration) + `analysis.tdf_bin` (compressed binary frame payloads).
`tdf_bin` is a **concatenation of self-delimiting per-frame blocks** — the current encoder emits a
4-byte block length, then a 4-byte scan count, then the compressed data — and `Frames.TimsId` is
simply the **byte offset** of each frame, written as the current file position. So the writer is:
*append the block → insert the `Frames` row with `TimsId = position` in one transaction.* **No global
binary index, no final offset pass.** (AlphaTims reads a frame by seeking to `TimsId` and reading its
length header — independent confirmation the format is offset-addressed, not index-addressed.)

What still must be gotten right, and tested against the target reader:

- **`(scan, tof)` must be sorted and deduplicated per scan** before delta-coding — v1's writer already
  notes unsorted TOFs produce invalid delta coding. The per-scan accumulation buffer sorts on emit.
- **Compression type must match `GlobalMetadata.TimsCompressionType`** (v1 uses zstd); it is per-frame
  and streamable.
- **End-of-run SQLite finalisation** — global aggregates (run/frame counts, mz/mobility ranges) are
  SQL `UPDATE`s after the sweep; they do **not** require retaining frame payloads.
- **Crash / partial-output policy** — an interrupted render must leave a `.d` that is *explicitly*
  invalid or resumable, never one that *looks* complete (e.g. write a `finalised` marker last, or
  render to a temp path and atomically rename on success).
- **Reader interop is a test, not an assumption** — the vendor SDK, AlphaTims, and DIA-NN must open
  and round-trip the emitted `.d` (§8).

---

## 6. Noise, reference frames, DIA vs DDA

- **Noise / reference superimposition.** v1 can superimpose simulated signal onto a *real* blank `.d`'s
  frames (§8.6 of the plan — real peptides leak in). In a streaming render this means, per frame,
  reading the reference frame's peaks and adding them. Streamable if the reference is read frame-by-
  frame in lockstep; a memory risk if it is loaded whole.
- **The sweep consumes an `AcquisitionEvent` stream, not a hardcoded frame role — and each event
  carries its TIMING.** This one seam makes DDA *not a rewrite* and makes non-IMS timing correct. An
  event carries `{ start_time, end_time/exposure, role (MS1 | MS2), isolation geometry, collision
  energy, encoding }`. Critically, **`elution_weight` integrates the RT profile over the event's
  exposure interval `[start, end]`, not over a frame index.** This is the real IMS/non-IMS unifier: a
  timsTOF frame is a ~100 ms exposure; a fast Astral MS2 is a much shorter one; uneven MS1/MS2 cadence
  or a nonuniform template schedule must change total signal *through the exposure integral*, not
  through frame-index discretization (which would leak or inflate signal purely from the schedule).
  The sweep, active-precursor store, and writer are identical regardless of how events are produced.
- **DIA (first target)** supplies the event stream **statically** — a fixed, periodic schedule known
  up front.
- **DDA (follow-up)** supplies it from an **online controller**: after each MS1 frame is rendered and
  centroided, a selector consumes that (noisy, realized) MS1 evidence + a dynamic-exclusion list +
  timing, and *emits* the next MS2 events. Because selection depends on the *realized* MS1 signal, a
  two-pass "render all MS1, then pick, then render MS2" would be *less faithful* (it would select on a
  cleaner signal than the instrument sees). The online-controller-behind-the-event-stream design is
  both more faithful and keeps the streaming architecture. DDA is deferred, but the event-stream seam
  is designed in now so deferring it costs nothing.

---

## 6b. Instrument as a first-class abstraction (non-IMS is not a special case)

Today, non-IMS instruments are bolted on as a **parallel path**: `astral_dispatch.rs` +
`AstralAcquisitionBuilder` + `ThermoRawWriter` render a separate DB to a `.raw`, with capability flags
force-set (isotope transmission off, NCE-per-window). The Bruker IMS path is "the render" and
Thermo/Sciex are exceptions grafted on. That is the hackiness to remove.

**And it is deeper than a writer fork.** The current non-IMS builders **fabricate a 451-scan mobility
grid** and marginalize it away later — so "non-IMS has no mobility" is *not* true upstream; it is a
pseudo-mobility representation carried and then collapsed. v2 must make `n_scans = 1` **real at the
source** (the sweep emits a single scan for a non-IMS instrument), not a writer-side reinterpretation
of a fabricated grid. Otherwise the special-casing simply moves; it does not disappear.

The streaming sweep makes unification natural, because **ion mobility is just one axis that may have
one bin or many.** One render kernel, parameterised by an `Instrument`:

| capability | Bruker timsTOF | Thermo (Orbitrap/Astral) | Sciex |
|---|---|---|---|
| **mobility axis** | many scans (Gaussian per ion) | **1 scan** (collapsed) | 1 scan |
| **transmission `Q`** | scan-dependent **diagonal** (diaPASEF) | **static m/z window** (scan-independent) | static window |
| **CE unit** | eV | NCE | (instrument-specific) |
| **isotope transmission** | quad-dependent | none | none |
| **output writer** | TDF `.d` | Thermo `.raw` / mzML | mzML |
| **calibration** | tof↔m/z, scan↔mobility | tof/peak↔m/z, no mobility | ↔m/z, no mobility |

The sweep (§3.4) is **identical** across all of them. What the `Instrument` supplies:

1. **Mobility axis size.** IMS → the ion spreads over its scans via `mobility_weight`. Non-IMS → **one
   scan**, `mobility_weight ≡ 1`, and the whole ion lands in that single scan. Not a code fork — a
   parameter (`n_scans = 1`) that makes the mobility factor degenerate.
2. **Transmission `Q(mz, scan, frame)`.** Bruker → the existing diagonal (`mscore::quadrupole`).
   Non-IMS → a static window `Q(mz)` (no scan term), which is a *special case of the same interface*,
   not a different pipeline. So `m/z ∈ W` — the thing that was *wrong* for diaPASEF — is *correct* for
   Thermo/Sciex, and both are the one `Q` interface at different settings.
3. **CE unit + fragment model.** The model registry keyed on the instrument (the `fragment_predictor_capability`
   guard already refuses feeding NCE to an eV model) picks the CE encoding and fragment model.
4. **Isotope transmission mode.** A capability flag consumed by the payload step, not a forked builder.
5. **An `AcquisitionLayout` + a pluggable writer — not just a serializer.** A writer alone is too thin,
   because vendor formats carry *format physics*, not just bytes: Thermo `.raw` authoring is
   **template-slot-bound** (it uses each template scan's profile grid and calibration, and can
   defer/repack/clear overflowing payloads); mzML must declare **profile vs centroid** correctly
   (real non-IMS DIA can use profile MS2, and v1's mzML currently declares everything centroided); the
   mzML activation field records a CE value whose conventional unit is eV even when it is an NCE; and
   scan-event metadata (RT units, MS level continuity, precursor/activation fields, polarity/analyzer,
   packet/grid constraints) all matter for whether DIA-NN/Spectronaut accept it. So the `Instrument`
   supplies an **`AcquisitionLayout`** — profile/centroid mode, calibration, per-event scan metadata,
   exposure timing, encoding — that the writer serializes. v1 already has `AcquisitionWriter` /
   `ThermoRawWriter` / the mzML writer to consolidate behind this.

The payoff: **cross-instrument becomes the same feature space measured by a different `Instrument`** —
the exact pattern CCS/RT already proved (one CCS, different gas → different 1/K₀). One HYE feature
space rendered as timsTOF-diaPASEF *and* Astral-DIA *and* a Sciex run (Sciex output is **mzML**, not a
`.wiff`), by swapping the `Instrument`, not by three code paths. The non-IMS "hack" disappears because
non-IMS is the `Instrument` with `n_scans = 1` (real, not fabricated), a static-window `Q`, and its
`AcquisitionLayout`.

**Out of scope (flagged, not designed):** genuinely *multiplexed / coded* isolation — scanning-SWATH
or overlapping-multiplexed acquisition where several isolations are summed into one spectrum and must
be deconvolved — is not a static `Q` nor a diagonal; it needs a transmission/mixing matrix and
deconvolution-aware metadata. Ordinary overlapping (non-multiplexed) Astral windows are fine: an
event-indexed static `Q` gives each event its own signal, so overlap needs no special geometry.

## 7. Sequencing (prove the win before the heavy port)

1. **Sweep-line prototype + memory benchmark.** ✅ **DONE** — `timsim-cli/src/bin/render_bench.rs`
   (`timsim-render-bench`). Reads the Parquet feature space, runs the active-set sweep, accumulates a
   real per-frame sparse `(scan, tof)` buffer from each active precursor's isotope envelope (the MS1
   path), emits frame statistics, writes no TDF. Run against the 11.9 M-precursor e2e artifact.

   **Measured — the central claim holds.** Fixing elution *density* and **doubling the run length**
   (3 000 → 6 000 frames, 30 k → 60 k precursors) left the sweep's working set flat:

   | metric | 30 k / 3 000 fr | 60 k / 6 000 fr |
   |---|---|---|
   | peak **active set** | 2 022 | 2 033 |
   | peak frame buffer (scan,tof) | 121 313 | 119 887 |
   | **sweep working set** (active heap + buffer) | **3.73 MB** | **3.69 MB** |
   | MS1 conservation (worst rel. err) | 1.15e-9 | 1.26e-9 |

   So the render's own memory is **independent of run length** — the active set (~2 k) is ~15× below
   the total even on this small slice, and would be ~200× below on the full run. **But read the claim
   precisely** (Codex flagged the prose as over-broad): what is proven is independence from run
   *duration* at fixed density. The working set is

   > `working set = f(local arrival rate, elution-width distribution, mobility width, ion/fragment
   > payload, acquisition geometry, peak density)` — **not** `f(run length)`.

   "Fixed density" deliberately holds the quantity that sets active-set size constant, so this is an
   independence result, not a universal few-MB *capacity* bound. Turning it into a capacity claim is
   what the **stress fixture (§4) is for, and it is therefore not optional**: a matrix over local
   density, elution width/tails, mobility width, isotope depth, and fragments-per-precursor, reporting
   **p50/p95/p99 and max** of local arrival rate / active precursors / unique bins / RSS — not just
   averages or one peak. Known breaking regimes to cover: a concentrated RT/enrichment burst (raises
   local arrival rate), broad or long-tailed elution (raises overlap at identical average density),
   gradient-edge compression of many `rt_index` into a short region, and the "density grows with run
   length" workload (more material per second on a longer run — flatness is then neither expected nor
   required). MS2 can be materially worse than MS1: fragment multiplicity and co-isolation, not
   precursor count, drive `(scan,tof)` occupancy.

   Conservation: the bench's old check was only a **consistency** check (Codex's sharpest catch) —
   `expected` built from the *same* truncated, clamped, frame/scan-partitioned integration the
   emission uses, so a fault duplicated in *both* paths would pass. That check is **removed from the
   bench** and replaced by an independent oracle, ✅ **DONE** in `timsim-cli/src/render.rs`
   (`cargo test -p timsim-cli --lib render`). The production sweep is now a shared function
   (`stream_render`) that both the bench and the oracle drive, so the tested code is the code that
   runs. The oracle is independent in *both* decomposition and code:
   - `reference_render` — an **ion-major** render (the sweep is frame-major) that discovers each ion's
     frame window by a direct `fs..=fe` loop (the sweep discovers it via heap enter/leave against a
     moving cursor), with no active set, no per-frame buffer, no input sort, and window bounds
     recomputed by a different expression. It shares only the pure Gaussian weight, which is
     unit-tested on its own. `sweep_render` vs `reference_render` is compared at **every bin**
     (`test sweep_matches_independent_reference_every_bin`), worst rel-diff < 1e-12.
   - Metamorphic tests: duplicate one ion → exactly 2× *its* bins and nothing else; permute input
     order → identical render; union of two overlapping-chunk renders == a single render; plus a
     lone-ion analytic mass check computed with no render code.
   - **Teeth verified by mutation**: flipping the sweep's leave-condition by one frame (ion lingers an
     extra frame) fails both the every-bin and the analytic test — the class of enter/leave off-by-one
     the check exists to catch.

   Still owed for MS2 (when fragments + transmission enter the render): the same independent-oracle
   treatment with an independent recomputation of diagonal `Q(scan, mz)` including soft-edge points,
   exercising all three transmission modes — a total-intensity check is insufficient. The `Ion` type
   the render consumes already generalises (MS1 isotopes and MS2 fragments are both just a locus + a
   set of `(tof, intensity)` peaks), so the MS2 oracle reuses this harness.

   **Two findings that shape the real writer:**
   - *Peak RSS (~967 MB) is dominated by the load, not the sweep.* The prototype materializes the whole
     `Vec<Precursor>`, builds an O(total) in-memory RT-index map (3.5 M entries), and lets Arrow
     decode full row groups of a 396 MB Parquet with no column projection. The sweep itself needs only
     ~4 MB. The production render must (a) project only the needed columns, (b) stream precursors
     **pre-sorted by `frame_start`** off disk rather than holding them all, and (c) not materialize the
     full RT map. Those three are the difference between "< 1 GB by a hair, for the wrong reason" and
     "flat few-MB working set regardless of run size." Fragment payloads must be **co-sorted /
     denormalised with (or partitioned identically to) the render artifact**, or the streamed read
     trades RAM pressure for random-I/O stalls.
   - *The `HashMap<(scan,tof)>` per-frame buffer is the throughput bottleneck* — ~10 M peaks/s, ~92
     frames/s. At the full run's real density (~330 precursors/frame) that is hours. The real render
     needs a **dense per-scan buffer** (a `Vec` indexed by tof, reused across frames) — the same shape
     v1's `from_tims_spectra_filtered` already uses. Two caveats (Codex): do **not** clear the full
     TOF range every frame — track touched-TOF lists/epochs so clearing is O(touched), not O(range);
     and **include the buffer's fixed allocation in RSS reporting**. It aligns with TDF's per-scan
     sorting/delta coding and conflicts with neither diagonal transmission (which only decides *which*
     contributions are added) nor DDA (which only changes the event stream online) — but it makes the
     buffer-cap/spill policy more consequential, since dense storage bounds *allocation* but not
     emitted peak count or compression time.
2. **Production-shaped throughput + RSS benchmark (still writer-free).** ✅ **DONE** —
   `timsim-render-bench --encode` (build `--features bench-encode`), driving the real Bruker TDF block
   encoder (`rustdf::data::utility::reconstruct_compressed_data`: sort-dedup + delta-TOF + zstd).
   Measured on the 11.9 M-precursor artifact, ~9.5 k active precursors/frame (a co-elution-heavy
   slice), on a 16-core / 31 GB box:

   > **Correction (post code-review).** An earlier version of this section reported flat MS1 at ~12 min
   > and called zstd "the MS2 bottleneck." A Codex review caught a quantisation bug — the flat encode
   > path quantised each contribution to `u32` *before* summing duplicates, so co-eluting sub-quantum
   > signal was floored to zero. Fixing it (sum in f64, then quantise, in the parallel worker) both
   > made the two accumulators produce **byte-identical output** and revealed the numbers below. The
   > true dominant cost is the **dedup** (billions of contributions → unique `(scan,tof)` bins), not
   > zstd; the old "encode bottleneck" was the dedup happening *inside* the encoder's `sort_dedup`.

   **The accumulator is the render-*pass* lever, but the dedup is unavoidable.** Swapping the per-frame
   `HashMap` for a flat push-`Vec` makes the *emission* pass **6.0× faster** (render-only: hashmap
   105.8 s vs flat 17.6 s at 150 k/3000; `Vec::push` vs `HashMap::entry`). But reducing every
   contribution to unique bins must happen somewhere — flat defers it to the (parallel) encode stage
   rather than fusing it into the render. Both accumulators are bin-identical (proven:
   `flat_accumulator_dedups_to_the_same_cube`) and produce identical encoded output.

   **Where the time actually goes (flat MS1, 150 k/3000, dedup+encode on 16 threads):**
   | stage | this slice | full 36k-run extrapolation |
   |---|---|---|
   | emission (push) | 42 s | ~8 min |
   | **dedup + quantise + zstd** (parallel) | 65 s | ~13 min |
   | end-to-end | 107 s | **~21 min** (16 threads) |

   So the dedup+encode back end is the larger half, and — the key point — it is **per-frame
   independent, hence parallel**. Moving dedup out of the serial render callback into the parallel
   encode worker is what makes this scale (an intermediate version that deduped serially was ~134 s and
   MS2 did not finish). Memory is the trade: the parallel worker buffers raw triples, so RSS rose to
   ~3 GB at `--encode-batch 64` (tunable down for dense MS2, up on a big-RAM box).

   **Parallelism headroom (for the 128-core / 512 GB target).** Both stages are per-unit independent:
   - *Dedup + quantise + zstd by frame* — ✅ measured parallel (`--encode-threads`). This is the larger
     half of the cost, and it scales near-linearly with cores (frames are independent). The bounded
     `--encode-batch` of raw triples is what feeds the workers; raise it on a big-RAM box.
   - *Emission (render) by RT chunk* — ✅ **DONE and measured** (`--render-chunks K`). The run splits
     into K contiguous frame ranges rendered concurrently (`stream_render_flat_range` on a rayon pool);
     ions are bucketed into the chunks their active window overlaps, and each chunk emits only its own
     frames so the pieces concatenate (no summing, no double-emit). Correctness is pinned by a *second*
     oracle invariant — `frame_range_partition_equals_whole` (frame-partition, distinct from the
     ion-partition `chunk_union_equals_whole`), every bin, worst < 1e-12. Measured on the 16-core box
     (**emission only** — the dedup is measured separately in the encode stage, not folded in here):

     | | serial (K=1) | K=16 | speedup | full-run emission |
     |---|---|---|---|---|
     | MS1 flat | 17.3 s | 3.26 s | **15.0× (94%)** | 3.5 min → **0.7 min** |
     | MS2 flat (26 pk/ion) | 12.9 s | 4.35 s | **14.7× (92%)** | 7.8 min → **2.6 min** |

     Two honest caveats the measurement exposed:
     - *Memory-bandwidth contention, not emission duplication.* Total peaks emitted is invariant in K
       (emission is partitioned, never redone), but `Σ per-chunk` wall grew (17 s → 49 s at K=16) as
       threads contend for RAM bandwidth. So per-core render throughput falls as cores rise — the 94%
       efficiency on 16 cores will be **sub-linear at 128**; project the 128-core number from
       bandwidth, not a naive ×128.
     - *Halo tax when chunks get small.* A boundary ion is bucketed into every chunk its ±n_sigma
       window touches; when chunk length drops **below** the elution window the bucketing dup-factor
       blows up (measured 3.6× at 62-frame chunks vs a 180-frame halo). Keep chunk length ≫ the
       elution window (for a 36 k-frame run, K=128 → 281-frame chunks vs ~180-frame halo is already
       marginal; fewer, larger chunks or a smaller halo is better). This inflates bucketing/sort only,
       not emission, but it is real overhead.
   - *512 GB dissolves the load confound* — the ~1 GB peak RSS is the O(total) load, not the sweep, so
     with that much RAM you simply hold the whole precursor+fragment set resident and skip
     streaming-from-disk entirely; each render chunk's active set is a few MB.

   **Projected full run on 128 cores / 512 GB: a few minutes wall (vs tens of minutes single-node)**,
   floored by the genuinely serial parts (Parquet load; in-order block append to `tdf_bin`, which is
   cheap) and by RAM bandwidth (see the render-by-chunk caveat above — do not assume ×128). So "prove
   the win, then port" is vindicated: the renderer is viable, the true cost (the **dedup**, not zstd)
   is understood, and **both halves scale near-linearly with cores** (emission 15×/16 by chunk;
   dedup+encode by frame). Still owed before/with the port: add real diagonal transmission (it only
   *removes* MS2 points, so these numbers are a conservative upper bound), and fuse the two parallel
   stages into one pipeline (measured separately today). The single most load-bearing correction the
   review produced: **quantise only after summing** — carry f64 through the dedup, never floor a raw
   contribution.
3. **Incremental TDF writer** in `rustdf`, driven by the proven sweep (§5). Fan-out follows the
   measured axes above: parallel render-by-chunk feeding the already-parallel block encoder, serial
   only at the final ordered append.
4. **Wire in + parity.** Replace the SQLite lazy builder on the v2 path; validate the streaming `.d`
   against v1's render on the same feature space (see §8).

---

## 8. Validation (aggregate tests are too forgiving — layer them)

TIC correlation + search recovery are aggregate and *forgiving*: a search can identify the same
peptides while isolation leakage, mobility placement, cofragmentation, intensity rank, and calibration
are all subtly wrong. Keep them, but they are the *top* of a layered suite, not the whole thing:

1. **Golden micro-fixtures (exact).** A handful of hand-built peptides, noise-free, with the *exact*
   expected per-bin `(frame, scan, tof, intensity)` after quantisation. This is the bit-level pin —
   catches an off-by-one in scan/tof placement or a transmission error that aggregates hide. (The same
   discipline as the fragment-decode pin test.)
2. **Conservation.** Per precursor per charge, emitted intensity sums to `abundance · charge_fraction`
   up to the measured truncation and the modelled `Q < 1`.
3. **Per-precursor distributional checks.** 2-D RT×mobility peak moments and widths, isotope ratios,
   fragment rank/intensity distributions, and **isolation-edge transmission curves** (does a precursor
   at the window edge get the soft-`Q` partial transmission, not full or zero?).
4. **Per-frame checks.** Peak count, base-peak, saturation/clipping behaviour, m/z and mobility
   marginal distributions, and MS2 interference / chimericity metrics (a naive window would show *too
   little* chimericity — this is where the diagonal `Q` earns its place).
5. **Reader interoperability — import/search behaviour, not just "it opens".** The vendor SDK,
   AlphaTims, and DIA-NN/Spectronaut must *import and search* the emitted TDF (and the `.raw`/mzML for
   non-IMS) and recover the ground truth — testing isolation offsets, RT units, MS-level continuity,
   precursor/activation metadata, and **whether profile MS2 is retained when the template uses it**,
   not merely that the file parses. Plus TimsId-offset / metadata consistency checks.
6. **Parity with v1** on the same feature space: per-frame TIC correlation + per-precursor peak-set
   agreement. Bit-exactness is *not* expected (different accumulation order).
7. **Search accuracy, stratified** — not just total IDs, but IDs sliced by abundance, charge, mobility,
   **window-edge proximity**, and coelution. A render that is only wrong at window edges or in dense
   coelution shows up here and nowhere in the aggregate number.
8. **Determinism.** Given the feature space + seed, the render is bit-identical run to run
   (identity-keyed noise; deterministic reduction order even under cap/spill).

---

## 9. Open questions (after review — the resolved ones removed)

Resolved across two review rounds: **TDF is append-only** (self-delimiting blocks, `TimsId` = offset,
§5); the **transmission diagonal already exists** and is reused as **one mode**, de-double-counted
(§3.5); the **precompute is a sorted artifact** (§3.6); **timing is per-event exposure**, not frame
index (§6); the **MS1/MS2 conservation split** is correct; ordinary **overlapping Astral windows** need
no special geometry (event-indexed static `Q`); and non-IMS `n_scans = 1` must be **real upstream**,
not a fabricated 451-scan grid collapsed at the writer (§6b). Still open:

1. **Fragment locality.** Denormalise fragment intensities into the sorted `render_precursors`
   artifact (pure sequential I/O, larger artifact) vs keep `fragment_intensities` separate and
   random-access it for active precursors (smaller artifact, random reads). The sweep only ever needs
   the active set's fragments, so a co-sorted/co-partitioned layout likely wins.
2. **Detector saturation / dynamic range.** Modelled as a post-accumulation per-scan clip, or skipped
   in v1? It affects high-abundance peaks and the realism of the top of the dynamic range.
3. **Reference-frame superimposition** (§6) streamed in lockstep — the reference `.d` must be read
   frame-by-frame, not loaded whole, and its frame grid must align with (or be resampled to) the
   simulated schedule.
4. **The cap/spill worst case.** What is the realistic peak active-set size under a real enrichment or
   a plasma-depletion-style narrow-RT band, and does the per-frame peak cap interact badly with it?
5. **Non-IMS `AcquisitionLayout` fidelity (§6b).** Does the layout carry enough (profile-vs-centroid,
   calibration, scan-event metadata, exposure timing) that DIA-NN/Spectronaut *search* the emitted
   `.raw`/mzML as they would a real file? The MS1/MS2 profile-retention path is the likeliest gap.
6. **`PerFragment` soft-edge fidelity (§3.5).** Is the thresholded transmitted-isotope selection good
   enough for MS2 benchmark realism, or is a probability-weighted isotope-distribution variant worth
   building? Decide with a fixture at a window edge.
7. **Multiplexed / coded acquisition.** Out of scope for v1 of the render (§6b), but if scanning-SWATH
   or overlapping-multiplexed methods become a target, the transmission model needs a mixing matrix +
   deconvolution metadata — a real extension, not a parameter tweak.
