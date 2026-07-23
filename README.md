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
│                                render (Bruker .d, v2) / render-sciex (SWATH mzML, v2)
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

## Layout
```
flow/timsim_flow.py     the DAG (nodes, typed edges, command wiring)
flow/configs/           run configs (design*.toml, hela*.toml, hye.toml, sciex.toml, mods*.toml, fasta)
golden/                 the regression gate — a frozen 60-protein run scored two ways (sim realism +
                        tool benchmark) against its own answer key; `golden/run.sh`, see golden/README.md
requirements.txt        the lean Python dependency surface (necroflow + timsim-predict + timsim-eval)
Makefile                predict-deps / rust-bins / setup
```
