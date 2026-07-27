# timsim v2 — MS2 / DIA render scope

The MS1 `.d` writer + reference-calibrated render is done and opens in DiaNN's reader (openTIMS/SDK,
sub-ppm m/z). This scopes adding **MS2 (DIA-PASEF fragment) frames**: fragments emitted into MS2 frames,
gated by the diagonal quadrupole transmission. Grounded in a full read of v1 (`rustdf/src/sim/dia.rs`,
`mscore/src/timstof/quadrupole.rs`, `acquisition.py`) — the headline is that **most of the physics
already exists in Rust; the new work is a scan-resolved sweep emitter and wiring.**

---

## 1. What MS2 adds over MS1

An MS1-only run places each precursor's isotopes into every frame of its elution window. DIA adds a
second frame type interleaved on a cycle:

- **Frame types.** `MsMsType = 0` (MS1/precursor) every `precursor_every` frames; `MsMsType = 9`
  (`FragmentDia`) otherwise. `ScanMode` stays 9 (DIA). Reference cycle: `precursor_every ≈ 17`.
- **The diagonal.** Each MS2 frame belongs to a **window group**; the group defines, *per mobility
  scan*, an isolation `(m/z, width)`. So an MS2 frame transmits *different* precursor m/z at different
  scans — a diagonal stripe in (mobility, m/z). A precursor at scan `s`, m/z `m` is fragmented in
  **every** MS2 frame whose group's window at scan `s` transmits `m` with non-negligible probability —
  which, because adjacent diagonal windows overlap and the logistic edge is soft, can be **more than
  one group**. The expected MS2 signal is the **sum over all transmitting frames/groups**. Co-isolation
  (multiple precursors transmitted at the same scan contributing to the same MS2 spectrum) is real and
  must occur — it is a feature of DIA, not to be suppressed.
- **Fragments, not isotopes.** In an MS2 frame a transmitted precursor deposits its **b/y fragment
  ions** (at the fragments' own m/z → TOF), scaled by the precursor abundance × elution × mobility ×
  **transmission**, not its precursor isotopes. Fragments are gated by the precursor's transmission,
  **never re-gated at their own (lower) m/z** — they sit below the isolation window.

So the render gains: (a) a frame-type layout, (b) a per-`(frame, scan, precursor-m/z)` transmission
gate, (c) a fragment peak set per precursor.

**First-version idealization (state it plainly).** v1 of the MS2 render is **b/y-only, complete
fragmentation, fixed collision energy** — no surviving precursor, no neutral-loss/internal ions, and a
constant fragmentation *yield* across windows/scans (CE-dependent yield and CE-dependent relative
intensities are deferred; see §4). This is enough for a valid, gating-correct DIA file; it is *not* yet
a claim of spectrally realistic MS2 for robust search scoring.

---

## 2. What we REUSE (already in Rust — do not rebuild)

- **The diagonal gate: `mscore::timstof::quadrupole::TimsTransmissionDIA`** (`quadrupole.rs:431`).
  Constructed from the window tables:
  `new(frame, frame_window_group, window_group, scan_start, scan_end, isolation_mz, isolation_width,
  k=Some(15))`. It **expands each group across every scan** in its range internally. Key methods:
  `apply_transmission(frame_id, scan_id, &Vec<mz>) -> Vec<f64>` (per-m/z probability, soft logistic
  edge), `is_precursor(frame_id) -> bool`, `get_setting(group, scan) -> Option<(mz, width)>`. This
  *is* the diagonal — reuse verbatim.
- **Fragment m/z from sequence** (`mscore/src/data/peptide.rs`): `calculate_product_ion_series(charge,
  FragmentType)` (`:319`) builds b/y ions; `PeptideProductIon::mz()` (`:195`) gives each fragment's
  m/z. So m/z is derivable from `(sequence, ion_type, ordinal, frag_charge)`.
- **The three scaling modes** (`IsotopeTransmissionMode`: `None` / `PrecursorScaling` / `PerFragment`)
  and the `fragment_series_spectrum` kernel (`dia.rs:38`) — reference for how transmission scales
  fragment intensity once gated.
- **DIA table readers** (`meta.rs`): `read_dia_ms_ms_info` → `DiaMsMisInfo{frame_id, window_group}`;
  `read_dia_ms_ms_windows` → `DiaMsMsWindow{window_group, scan_num_begin, scan_num_end, isolation_mz,
  isolation_width, collision_energy}`. These let us copy a real DIA schedule from a reference `.d` —
  exactly the pattern the calibration copy already uses.
- The MS1 sweep, oracle, parallel chunking, and the TDF writer (calibration copy) — all as-is.

**Reuse note:** v1's own MS2 render (`build_fragment_frame`, `dia.rs:575`) is frame-by-frame and reads
a *simulation* SQLite (`fragment_ions`, Prosit-174 blobs), not our Parquet feature space. We reuse its
*physics* (transmission gate, scaling modes, fragment m/z) but not its frame loop — v2 renders with the
sweep.

---

## 3. What we BUILD

1. **Spectra from the peptide ion via mscore (✅ DONE — `timsim-cli/src/spectrum.rs`).** *Superseded the
   original "materialize a `fragment_mz` column" plan* after two decisions:
   - **Reuse v1's peptide-ion path, don't hand-roll.** mscore already turns a peptide ion into an
     isotopic spectrum: `PeptideIon::calculate_isotopic_spectrum` (precursor) and
     `PeptideSequence` → `calculate_product_ion_series` → `generate_isotopic_spectrum` (fragments).
     `spectrum.rs` wraps both as `precursor_peaks` / `fragment_peaks`.
   - **Unify precursor AND fragment on this one path.** The MS1 render previously hand-placed precursor
     isotopes at `mono_mz + k·1.0033/z` from the `isotope_intensity` envelope — approximate and
     inconsistent. Both now come from mscore (exact isotope m/z from composition). This makes the
     render **sequence-driven** (it needs the annotated peptide sequence to build the ion) — which MS2
     required anyway — and retires the `1.0033/z` placement.
   - **Attach per-ion, never re-flatten.** v1 feeds mscore the Prosit flat-174; our decoded
     `fragment_intensities` `(ion_type, ordinal, frag_charge, intensity)` is attached directly onto the
     computed ion series (`n_ions[i]`=b_{i+1}, `c_ions[i]`=y_{i+1}), mirroring
     `associate_with_predicted_intensities` — so the two-incompatible-flat-174-layouts trap can't
     reappear. A test pins the mapping (a `'b'` intensity lands on a `b` ion's m/z, `'y'` on a `y`).

   Consequence: **no new `fragment_mz` artifact** — the render generates spectra from the peptide ion
   (v1's approach), using the feature-space artifacts only for what mscore can't derive (the Prosit
   relative fragment intensities). The CE-metadata / fixed-CE-prior notes (§4) still stand.
2. **A frame-type-aware DIA sweep** (render module). Extend the ion model and emission:
   - `Ion` gains `precursor_mz` and `fragment_peaks: Vec<(tof, intensity)>` alongside the existing
     isotope `peaks`.
   - At each frame the sweep asks the layout: MS1 or MS2 (+ window group)?
     - **MS1 frame** → emit isotope `peaks` (today's path).
     - **MS2 frame** → for each scan in the precursor's mobility window, `t =
       transmission.apply_transmission(frame, scan, [precursor_mz])[0]`; if `t > ε`, emit
       `fragment_peaks` scaled by `abundance × elution_w(frame) × mobility_w(scan) × t`.
   - **No "one group" assumption.** Each MS2 frame is gated on its own; a precursor emits into *every*
     MS2 frame whose transmission clears `ε` at its scan — overlaps and soft edges mean that can be
     several. The total is the natural sum across frames. Likewise many precursors clear the gate at a
     given scan → their fragments co-accumulate into that MS2 spectrum (co-isolation, correct).
   - The diagonal still makes this **sparse**: a precursor clears the gate only in a narrow scan band
     of a few MS2 frames, so most `(precursor, scan)` pairs gate to ~0 and are skipped (subject to the
     `ε` cutoff — a logistic is never exactly 0, so `ε` is a required policy knob, not optional).
     Working-set and chunk-parallelism carry over from MS1.
3. **DIA schedule + frame layout — replay the reference's cycle, don't invent one.** Extract the
   reference's *complete frame-level cycle* (its exact MS1/MS2 pattern, its exact frame→window-group
   sequence including group ordering/repetitions/omissions, its MS1 cadence, and its frame
   timing/accumulation if represented) and **repeat that whole cycle** over our `n_frames`, truncating
   if we end mid-cycle. Do **not** synthesize "MS1 every N + round-robin groups" — a real schedule is
   not that regular, and DiaNN sees the difference. Copy `DiaFrameMsMsWindows` (the group→scan-range→
   isolation definitions) verbatim and construct `TimsTransmissionDIA` from it; emit our own
   `DiaFrameMsMsInfo` for our replayed frame→group.

   **Scan-grid portability is a precondition, not a given.** `ScanNumBegin/End` are only meaningful on
   the *same* scan grid / mobility ramp. Since we already take `num_scans` from the reference (the
   calibration copy requires it), the grid matches by construction — but the writer must **validate**
   it (every referenced scan range ⊂ `[0, num_scans)`, every MS2 frame maps to a defined group) and
   **reject** a reference whose grid or mobility-ramp orientation differs rather than silently
   distorting the diagonal. Scan direction is a known footgun — a flipped ramp inverts the diagonal.
4. **Writer extensions** (`TdfWriter`): write `DiaFrameMsMsInfo` (synthesised) and copy
   `DiaFrameMsMsWindows` (reference); set `MsMsType` per frame (0 vs 9) instead of a constant 0. The
   `Frames` copy-then-overwrite already carries `ScanMode=9`.
5. **The MS2 conservation oracle** (the correctness gate — see §5).

---

## 4. Architecture: one frame-type-aware sweep, not two passes

A **single sweep** over all frames, with the precursor active over its elution window (spanning MS1 and
MS2 frames alike), branching emission on frame type. Rationale: the precursor's elution/mobility
placement is identical for isotopes and fragments — only the peak set and the transmission gate differ.
One sweep keeps the active-set/working-set bound and the chunk-parallel decomposition intact (frame
ranges still partition cleanly; a chunk just needs the transmission table, which is read-only).

**Transmission model — staged, each stage isolated and testable.**

1. **Mono-m/z gating** (first): `t` = transmission of the precursor's monoisotopic m/z at
   `(frame, scan)`, one scalar, all fragments scaled by it. Fully exercises the diagonal geometry and
   is fine for broad windows — but it is **visibly wrong near diagonal boundaries**: precursor isotopes
   are `1/z` Th apart, so at a steep/narrow window edge different isotopes are selected differently.
   Worst for high charge and precursors near an edge — exactly the region diaPASEF exploits. Ship it
   **labelled as an approximation**, not as final.
2. **Isotope-aware parent transmission** (the first realism upgrade): scale the whole fragment series
   by the *effective* parent transmission
   `t_eff = Σ_i isotope_fraction_i × Q(parent_mz + i/z)` — the transmitted fraction of the actual
   isotope envelope, not a delta. This removes the boundary artifact.
3. **`PerFragment`** (later): models the **product-ion isotope envelopes after fragmentation** given
   the transmitted parent. It must **never** quadrupole-gate product m/z — fragments are below the
   window and are not re-selected; `PerFragment` only shapes fragment isotope structure.

**Collision energy is two separate effects, both currently constant** (call this out): CE changes the
*relative* b/y intensities *and* the *total* fragmentation yield. v1's fragment artifact is single-CE,
so v2 ships a **fixed-CE spectral prior** with constant yield across windows/scans. This is safe for
plumbing, file validity, and gating validation — but **not** a claim of realistic spectra or robust
DiaNN scoring. Preserve the copied per-window CE metadata now so CE-indexed fragments (a fragment-stage
extension) slot in later without a schema change.

---

## 5. The MS2 conservation oracle (gate before wiring)

Mirrors the MS1 oracle's independence discipline (an earlier Codex catch: a consistency check that
shares the emission's own math is tautological). Required, independent in decomposition and code:

- **Independent diagonal recomputation.** On a tiny fixture, recompute the transmission `Q(frame,
  scan, mz)` from the window definitions with a dead-simple logistic — *not* `TimsTransmissionDIA` —
  and compare every gated value. Recompute using the **exact edge convention `TimsTransmissionDIA`
  uses** (it may combine two logistic sides / interpret `width` a particular way, so "edge ≈ 0.5" is
  not assumed — read the convention off the code). Cover: window center, both soft edges, far outside,
  **scan-range boundaries, gaps, overlapping windows, group ordering, absent window settings, first/
  last frame of a cycle**, and 0- vs 1-based scan/frame indexing.
- **Transmitting-set, not "one group."** Independently compute the *set* of all `(frame, group)` with
  `Q > ε` for a precursor at its scan — including overlap and soft-edge contributions — and require the
  render to emit into exactly that set. (Replaces the wrong "exactly one group" invariant.)
- **Fragment mass conservation.** Total emitted fragment signal for a precursor equals `abundance ×
  (elution mass) × Σ_scan (mobility_w × t) × Σ_fragments intensity`, computed with no render code —
  **stated precisely**: whether fragment intensities are per-precursor normalized, and whether
  `abundance` is total-parent or monoisotopic. Account for `(scan, tof)` bin collisions in the check.
- **Fragment-not-regated negative test.** Move a fragment's m/z below/outside the isolation window →
  its emission must be **unchanged** (fragments are gated by the *parent*, never their own m/z).
- **Layout, not just intensity.** Assert the generated `DiaFrameMsMsInfo`, per-frame `Frames.MsMsType`
  (0 vs 9), and copied `DiaFrameMsMsWindows` match the intended replayed layout — MS1 frames have **no**
  `DiaFrameMsMsInfo` row and `is_precursor` returns true for them.
- **Metamorphic.** Duplicate a precursor → exactly 2× its fragment bins; a precursor outside every
  window at its scan → **zero** MS2 signal; frame-range partition still reproduces the whole (the MS1
  invariant, now with mixed frame types).
- **`is_precursor` consistency:** MS1 frames emit no fragments; MS2 frames emit no isotopes.

---

## 6. Data flow

```
peptides ─┬─ precursors (mz, charge, isotopes, abundance)     ┐
          ├─ peptide_rt (elution)                             ├─→ timsim-render (DIA)
          ├─ precursor mobility (m/z trend now; CCS later)    │      │ MS1 frames: isotopes
          └─ fragment_intensities + fragment_mz (NEW column)  ┘      │ MS2 frames: fragments × transmission
reference DIA .d ── DiaFrameMsMsWindows (copy) ─────────────────────┘   ↓
                                                                    MS1+MS2 .d  → openTIMS/DiaNN
```

---

## 7. Sequencing

1. **Spectrum building block (✅ DONE)** — `spectrum.rs`: `precursor_peaks` + `fragment_peaks` via
   mscore, attach-per-ion, unified precursor/fragment path. Next in this step: rewire the render to be
   **sequence-driven** (load annotated sequences, use `precursor_peaks` for MS1 in place of the
   `isotope_intensity` + `1.0033/z` placement) — the change that lands the unification.
2. **DIA schedule copy + frame layout** — read reference DIA windows, **replay the reference's full
   frame-level cycle** over `n_frames` (§3.3), build `TimsTransmissionDIA`; validate the layout
   (scan ranges ⊂ grid, every MS2 frame → a group) and unit-test `is_precursor`/`get_setting`.
3. **The MS2 oracle** (§5) — build it *before* the sweep emitter, as with MS1.
4. **Frame-type-aware sweep** — extend `Ion` + emission; precursor-mono-m/z gating first. Verify
   against the oracle every bin.
5. **Writer DIA tables** — `DiaFrameMsMsInfo` + `DiaFrameMsMsWindows` + per-frame `MsMsType`.
6. **Produce an MS1+MS2 `.d`, open in openTIMS/DiaNN.** The milestone: a DIA file a search engine can
   process — the point of the simulator.
7. **Refinements** (separately, each testable): isotope-aware gating → `PerFragment` → CE-indexed
   fragments.

---

## 8. Open questions

- **Reference DIA source.** G8602 is a real per-scan diagonal (20 groups × 918 scans); the coarser
  timsim-style layouts (e.g. 15 groups × few windows) are simpler. Copy which by default? (Lean: a
  real diagonal, so DiaNN sees realistic isolation.)
- **Cycle length** — derived from the reference's replayed cycle (not a CLI cadence knob; §3.3). A
  CLI override would mean re-timing the schedule, which changes the diagonal — out of scope for v1.
- **Mobility is still a placeholder m/z trend** until the CCS step. The diagonal *geometry* is correct
  regardless (transmission is a real function of the placed scan); only the scan *value* firms up with
  CCS. So MS2 does not block on CCS, and CCS slots in cleanly after.
- **Transmission cost at scale** — one logistic per `(active precursor, scan)` in MS2 frames. Sparse
  after gating, chunk-parallel; measure it in the throughput bench once wired (a `--dia` mode).

## 9. Implementation traps (Codex checklist — verify each in code, not by assumption)

- **MS1 frames have no `DiaFrameMsMsInfo` row**, and `is_precursor(frame)` returns true for them (it
  keys on *absence* from `frame_to_window_group`).
- **`ε` cutoff policy is mandatory** — the logistic is nonzero arbitrarily far from a window; without a
  cutoff every precursor emits a (tiny) tail into every MS2 frame. Define it once, apply it in both the
  render and the oracle.
- **Fragment charge ≤ precursor charge**, and the generated fragment m/z must match the exact
  convention the stored `frag_charge` assumes (proton mass, mono vs average).
- **Scan orientation.** Confirm the reference's scan→mobility direction; a flipped ramp silently
  inverts the diagonal. Validate, don't assume.
- **Pre-render layout validation.** Before emitting: every MS2 frame maps to a defined window group,
  and every referenced scan range lies within the target `num_scans`. Reject otherwise.
