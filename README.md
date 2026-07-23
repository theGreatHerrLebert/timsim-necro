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
│                                yield/modify/frag-input/spectra/
│                                render-thermo (Thermo .raw) / render (Bruker .d, v2)
│ prediction (Python, LEAN):     timsim-ccs / timsim-rt / timsim-fragments
│                                   └─ timsim-predict → pepdl → mscorepy
│                                      (mscore + ms-chem pyo3 primitives; imspy-free)
│ search (external):             DiaNN
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
make predict-deps          # necroflow + timsim-predict → pepdl → mscorepy   (the lean prediction stack)
make rust-bins            # cargo install timsim-cli from git (crates.io deps only, no rustims)
python flow/timsim_flow.py --help         # drive the DAG
```

`requirements.txt` is the whole Python surface — two git dependencies, nothing from imspy.

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
- **Bruker `.d`** (`--bruker-reference <ref.d>`) — **lean v2**: the same feature-space chain plus CCS →
  `timsim-render`, a streaming imspy-free projector onto a reference `.d`'s DIA schedule. Verified: a
  60-protein run authors a valid 3000-frame DIA `.d` (177 MS1 + 2823 MS2, 15k windows). The default
  Bruker path (no flag) still uses the v1 `timsim` monolith (imspy) — it owns DDA and the DIA truth
  output v2 does not emit yet, so scored Bruker runs stay on v1 until `run_dia` writes an answer key.
- **SCIEX mzML** (`--sciex-config`) — still the v1 `timsim` build-from-`.wiff` (imspy); no native v2
  SCIEX renderer exists.

**No gaps left on the Thermo path** — structure → prediction → render → search → score ingests only small,
independently-versioned federated repos, zero imspy / zero rustims. The lean Bruker `.d` render is wired
and verified; Bruker-DIA scoring and SCIEX remain the two v1 touchpoints.

## Layout
```
flow/timsim_flow.py     the DAG (nodes, typed edges, command wiring)
flow/configs/           run configs (design*.toml, hela*.toml, hye.toml, sciex.toml, mods*.toml, fasta)
requirements.txt        the lean Python dependency surface (necroflow + timsim-predict)
Makefile                predict-deps / rust-bins / setup
```
