# timsim v2 (necro) — capability roadmap

A **deliberate subset** of v1's surface, not a mechanical port. v1 (`imspy-simulation`) has ~100 config
knobs across simulator.py's `get_default_settings`; rebuilding all of them would resurrect the monolith we
just escaped. This doc tiers what to port, and — as importantly — what to **drop**.

Principle: port a capability when a real user/benchmark needs it, behind a metric ([[golden gate]]). "Would
a different MS tool want this?" → port. "Is it v1 plumbing / cosmetic / superseded by the DAG?" → drop.

## DONE — and already past v1's *validated* surface

The v1 **paper** (MCPRO-S-26-00597) validated **timsTOF only** (dia-PASEF, dda-PASEF, thunder-dda). We've
matched the timsTOF core AND gone beyond it:

- **Acquisition, lean + closed (render → search → score):** Bruker **DIA** (DiaNN) + **DDA** (Sage);
  **SCIEX** SWATH → mzML (DiaNN); **Thermo Orbitrap + Astral** `.raw` (DiaNN). Thermo/SCIEX were v1 *code*,
  never in the paper — so on acquisition we're ahead. The DIA render **replays whatever DIA scheme the
  reference `.d` carries**, so synchro/midia/Slice-PASEF come free with the right reference.
- **Sample:** HeLa (ID recall + FDP), scalable (5K validated, render ~10× faster + memory-gated).
- **Infra:** content-addressed DAG, GPU fragments, `--max-peptides`, the two-axis golden gate.

## P0 — MUST (the core realism + the top-two use cases)

1. **Realism: noise + spike-into-real** — the linchpin. Full design in **`REALISM_PLAN.md`**. Our clean
   FDPs (0.17–1.4%) reflect BOTH a noiseless render AND DiaNN 2.5 being well-calibrated. v1's FDR-inflation
   headline was **tool-specific** (Spectronaut on modifications; DiaNN 1.8's big discrepancies — neither
   tested here), so reproducing it needs noise **and** the tool axis (run the old versions). Noise is still
   P0: (a) **m/z-ppm scatter** (easy, all instruments); (b) **real-data-noise injection** from the
   reference `.d` (moderate, Bruker); (c) **spike-into-real-experiment** (`superimpose_on_reference`) —
   overlay synthetic ground-truth onto a *real* run, the strongest realism mode (elevated from P2 per
   David; only the spikes are labeled → needs a spike-recovery eval mode).
2. **HYE quant + fold-change eval.** The second-most-cited v1 axis, and cheap: the proteome is already
   configurable (`hye.toml`, HUMAN/YEAST/ECOLI). Needs (i) repoint FASTA paths (they broke when
   `SUBMISSION/zenodo` was cleared), (ii) a multi-condition design with the dilution ratios
   (HUMAN 0.65 / YEAST {0.15,0.30} / ECOLI {0.20,0.05}), (iii) a **quant-accuracy scorer** in `timsim_eval`
   (fold-change recovery per organism + cross-species leakage) — a genuinely new eval dimension beyond
   ID/FDP that generalizes to all instruments.

## P1 — SHOULD (distinct capability axes, clear demand)

3. **Phospho + FLR scoring.** PTM site-localization (a whole v1 benchmark). Needs: a phospho `mods` config
   (S/T/Y, ≥2 sites), a **phospho-capable fragment predictor** (v1 uses AlphaPeptDeep — Koina has it), and
   an **FLR scorer**. We can improve on v1 here: simulate both positional isomers in ONE run instead of
   v1's separate-runs-recovered-by-filename hack.
4. **HeLa complexity → 250k + the true-FDR curve.** We have 5K; scale the ramp and plot true-FDR vs
   density. Only meaningful **after noise (P0)** — the whole point is the FDR inflation, which needs noise.

## P2 — MAYBE (specialized; port on demand, not speculatively)

5. **HLA-I immunopeptidomics.** Peptide-*seeded* proteome (not FASTA-digest), non-tryptic, thunder-dda
   scheme, binomial 1+ charge. A distinct paper axis but narrow audience.
6. **MBR (match-between-runs).** Multi-run design + false-transfer eval (the PIP-ECHO split). Specialized
   DDA analysis; valuable but self-contained.
7. *(Plasma/PYE sample scenario — the spike-into-real *mechanism* is now P0 (`REALISM_PLAN.md`); Plasma/PYE
   is just a specific sample to apply it to.)*

## DROP — don't port (v1 plumbing, cosmetic, or superseded)

- **`from_existing`** — replay a prior `synthetic_data.db`. **Superseded** by the content-addressed DAG
  (re-request = cache hit; change a knob = targeted re-sim). No need.
- **`from_findings`** — drive the sim from real search results. Niche; not a simulation-fidelity concern.
- **Provenance / mzPROV** (Ed25519 signing embedded in the `.d`) — defer indefinitely; orthogonal to
  fidelity, adds a crypto dep.
- **Preview video generation** — pure cosmetic. Drop.
- **Waters SONAR** — fully-synthesized mzML (mirrors SCIEX, ~trivial) BUT very low demand. Build only if a
  collaborator actually needs Waters. Keep as a "known-cheap" note, not a task.
- **Explicit synchro/midia/Slice-PASEF configs** — come **free** via reference-`.d` scheme replay; no
  dedicated port.
- **Legacy intensity knobs** (`intensity_mean/min/max/value`) — already dead in v1 (`_LEGACY_IGNORED_KEYS`).
- **Binomial charge model** — we use site-specific; revisit only if a benchmark shows it matters.
- **v1's per-instrument builder sprawl** (waters/sciex/astral builders, `register_prediction_set`) — the
  lean render already covers the instruments we need; don't port the registry.

## Sequencing

P0 first (noise → then it's worth scaling HeLa and it makes HYE FDP real), then HYE quant (cheap, new eval
axis), then phospho. Each lands behind the golden gate. Revisit P2/DROP only when a concrete need appears.
