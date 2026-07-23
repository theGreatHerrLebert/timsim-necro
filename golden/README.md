# golden — the regression gate

A small, frozen reference experiment that both closed loops run end to end, so every change to the
pipeline lands against a **number** instead of a vibe. This is the spine of driving v2 toward v1 parity
without regressing: port a capability, run the gate, see whether detectable recall / FDP moved the way
you predicted — and whether anything else moved that shouldn't have.

## What it runs

The same 60-protein experiment through **both** measurement paths, each closing structure → render →
search → score:

- **Thermo `.raw` DIA** — `render_thermo` → DiaNN (native via .NET) → `timsim_eval` score
- **Bruker `.d` DIA** — the lean `timsim-render` → DiaNN (native, no .NET) → the same score

Because the DAG is content-addressed, the two loops **share their whole feature space** (digest → … →
spectra); it is computed once and only ccs (Bruker) + render + search + score differ. One gate run ≈ one
pipeline's worth of compute, not two.

## Use

```bash
./run.sh                    # run both loops, diff vs baseline.json, append to history.jsonl
./run.sh --only bruker      # one instrument (e.g. while tuning the Bruker render)
./run.sh --update-baseline  # rerun and REWRITE the baseline — do this deliberately, after a change
                            #   you've decided is the new reference (and say so in the commit)
```

Exit code is **1 on a regression** (detectable recall drops or FDP rises past 0.3 pp), so the gate drops
straight into CI or a pre-commit check. Example:

```
  Thermo DIA  recall(det)  58.3% ->  58.3% [no change]           FDP 1.63% -> 1.63% [no change]
  Bruker DIA  recall(det)  19.7% ->  24.9% [improved +5.2pp]     FDP 2.31% -> 2.10% [improved -0.2pp]

  ✓ no regression vs baseline.
```

## Layout

```
config/                 the FROZEN small inputs — version-pinned so the baseline is stable
  tiny.fasta              60 proteins
  proteome.toml.tmpl      proteome spec template (gate.py fills in the FASTA path — portable)
  mods_basic.toml         fixed carbamidomethyl + light oxidation (no PTM combinatorics)
  tiny_design.toml        single-organism, 50 proteins expressed
gate.py                 the harness: run both loops, parse metrics, diff baseline, log history
run.sh                  env wrapper (venv + TIMSIM_BIN + DiaNN/.NET), then gate.py
baseline.json           the pinned reference metrics (per instrument)
history.jsonl           append-only log of every run: {ts, commit, thermo{...}, bruker{...}}
```

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
