# golden — the regression gate

A small, frozen reference experiment that lands every change against a **number** instead of a vibe — the
spine of driving v2 toward v1 parity without regressing. The render's answer key is a *fixed oracle*;
which side you hold still decides what the number measures, so the gate has **two axes**:

| axis | fix | vary | catches |
|---|---|---|---|
| **sim** (default) | the tool (DiaNN) | the render / predictors / schema | simulation-realism regressions |
| **tool** (`--tool-axis`) | a frozen rendered dataset | the search engine / its version | tool regressions + cross-engine benchmarks |

They agree where they overlap: the tool axis run with DiaNN on the sim axis's own render reproduces the
sim-axis number **exactly** (Δ 0.00 pp) — a clean slice, and a correctness check on the freeze.

## Sim axis — is the simulation getting more realistic?

The same 60-protein experiment through **both** measurement paths, each closing structure → render →
search → score:

- **Thermo `.raw` DIA** — `render_thermo` → DiaNN (native via .NET) → `timsim_eval` score
- **Bruker `.d` DIA** — the lean `timsim-render` → DiaNN (native, no .NET) → the same score

The DAG is content-addressed, so within a run the two loops **share their whole feature space** (digest →
… → spectra) — computed once; only ccs (Bruker) + render + search + score differ. But content-addressing
keys on inputs/config/command, **not the binary**, so a rebuilt render binary at the same path would be a
cache hit. To read *current* code the sim axis **wipes the work dir each run** (use `--no-clean` only when
iterating and you don't need a true reading).

```bash
./run.sh                       # both loops, diff baseline.sim_axis, append history.jsonl
./run.sh --only bruker         # one instrument (e.g. tuning the Bruker render)
./run.sh --update-baseline     # rerun and REWRITE baseline.sim_axis — deliberately, and say so in the commit
```

## Tool axis — did the software regress / which tool wins?

Freeze a render you trust once, then run *only* search→score against it — cheap (~1–2 min, no
re-simulation), and it isolates the software from the simulation.

```bash
./run.sh --freeze                          # snapshot current renders as the frozen dataset (after a sim run)
./run.sh --tool-axis --tool diann          # DiaNN on the frozen dataset, diff baseline.tool_axis
./run.sh --tool-axis --tool diann --only bruker --update-baseline
```

Use it to (a) regression-test a DiaNN version bump (swap the binary, rerun), and (b) benchmark engines
head-to-head — `sage` / `fragpipe` slot into `gate.py`'s `TOOLS` registry via `timsim_eval`'s
`sage_executor` / `fragpipe_executor` (DiaNN is wired; the others are named slots).

Exit code is **1 on a regression** (detectable recall drops or FDP rises past 0.3 pp) on either axis, so
the gate drops straight into CI or a pre-commit check. Example:

```
  Thermo DIA  recall(det)  88.2% ->  88.2% [no change]         FDP 1.04% -> 1.04% [no change]
  Bruker DIA  recall(det)  26.6% ->  31.8% [improved +5.2pp]   FDP 0.87% -> 0.80% [improved -0.1pp]

  ✓ no regression vs baseline.
```

## Layout

```
config/                 the FROZEN small inputs — version-pinned so the baseline is stable
  tiny.fasta              60 proteins
  proteome.toml.tmpl      proteome spec template (gate.py fills in the FASTA path — portable)
  mods_basic.toml         fixed carbamidomethyl + light oxidation (no PTM combinatorics)
  tiny_design.toml        single-organism, 50 proteins expressed
gate.py                 the harness: sim axis + freeze + tool axis; parse metrics, diff baseline, log
run.sh                  env wrapper (venv + TIMSIM_BIN + DiaNN/.NET), then gate.py
baseline.json           pinned reference: {sim_axis:{thermo,bruker}, tool_axis:{inst:{tool}}}
dataset.json            manifest of the frozen tool-axis dataset (paths + content hashes)
history.jsonl           append-only run log (gitignored — per-machine telemetry)
```

The frozen dataset itself (`GOLDEN_DATASET`, ~0.9 GB of `.raw`/`.d` + truth) lives outside git; only its
`dataset.json` manifest is committed.

## The big templates are not committed

The Thermo template (`.raw`) and Bruker reference (`.d`) are hundreds of MB and live outside git. Point
at them with env vars (defaults suit the dev box):

| env var | what | default |
|---|---|---|
| `GOLDEN_THERMO_TEMPLATE` | a no-IMS Orbitrap/Astral `.raw` | `/scratch/timsim-demo/orbi_hela_dia.raw` |
| `GOLDEN_BRUKER_REF` | a reference DIA `.d` | `.../MIDIA-250K-NOISELESS-FRAG.d` |
| `GOLDEN_WORKDIR` | where the DAG materialises | `/scratch/timsim-demo/golden-work` |
| `NECRO` | venv with necroflow + timsim-predict + timsim-eval | `/scratch/timsim-demo/timsim-2/timsim-necro/NECRO` |

## What the numbers mean

- **recall(det)** — recall over the *strictest* denominator: precursors present, in a DIA window, with
  fragments, above the abundance floor. The honest "of what could be found, how much did DiaNN find."
- **FDP** — false-discovery proportion: DiaNN calls at q≤0.01 that don't match a truth precursor.
- **deciles** (in `history.jsonl`) — recall per abundance decile; the ladder should stay monotonic
  (high-abundance easier). A ladder that flattens or inverts is a realism regression the headline
  number can hide.

Absolute levels today: Thermo ~58% detectable recall, Bruker ~20% — the Bruker gap is a render-realism
target, and this gate is how you'll know when a change closes it.
