# timsim v2 — Snakemake port

A faithful [Snakemake](https://snakemake.github.io/) port of the necroflow DAG in
[`../flow/timsim_flow.py`](../flow/timsim_flow.py). Same tools, same command lines, same artifacts —
a side-by-side so you can compare the two orchestrators on the identical simulator.

Every resolved shell command is byte-identical to the one `timsim_flow.py` emits for the same
backend; this port only swaps *who assembles the DAG* (Snakemake instead of necroflow), not *what runs*.

## The idea it has to preserve

The necroflow flow's whole thesis is that the three axes are a **cache model, not a cross product**:

| axis | what it decides | cost |
| --- | --- | --- |
| **STRUCTURE** | which molecules exist | computed **once**, shared by every sample |
| **QUANTITY** | how much of each | the design; also **once** for the whole design |
| **MEASUREMENT** | how it is observed | **one job per sample** — the fan-out |

necroflow deduplicates the structure work with a content-addressed fingerprint. **Snakemake gets the
same dedup for free from the graph shape:** every rule above the measurement step has *no `{sample}`
wildcard*, so it has one output path and runs once no matter how many samples ask for it. Only the
`render_*` rule (and its downstream `search`/`score`) carries the `{sample}` wildcard.

```
proteome ─┬─ digest ─┬─ modify ── precursors ─┐
          │          │                        │
          │          └──── (ccs, rt) ─────────┼── frag_input ─ fragments ─ spectra ┐
          └─ design ── yield ─────────────────┘                                    │
                                 └───────────────────────── render_<backend>(sample) × N
```

You can see the collapse in the job counts from a dry run: request 2 samples on the `v1` backend and
you get **8 structure/quantity jobs (count 1 each) + 2 `simulate` jobs** — not 2 full pipelines.

## Invalidation: a spec file is a dependency, not a string

The flow is emphatic that hashing the *filename* `design.toml` (rather than its content) is the exact
failure it exists to kill — edit the fold change, and a filename-hash reports "up to date" and hands
you the previous experiment's run. necroflow solves it with per-node `hashes_file` invalidators.

Here the same guarantee is just an **honest input edge**: the spec TOMLs (`hye.toml`, `mods.toml`,
`design.toml`, the v1 config) are declared as rule `input:`s, so Snakemake restages exactly the rules
that read a file when its content changes. The Bruker reference gets the same treatment the flow gives
its hidden dependency:

- **Bruker v2** — the reference `.d`'s `analysis.tdf` (frame schedule + DIA windows + calibration) is
  declared as the render's input, mirroring `hashes_reference_d`.

(The lean SCIEX path takes no template file at all — its SWATH schedule is synthesised from command-line
params, so there is nothing hidden to hash.)

## Layout

```
Snakefile              config resolution, backend selection, `rule all`
config.yaml            every knob (mirrors the flow's argparse defaults one-for-one)
rules/structure.smk    proteome · digest · modify · precursors · ccs · rt · frag_input · fragments · spectra
rules/quant.smk        design · peptide_yield
rules/v1.smk           simulate           (backend: v1 — monolithic Bruker .d, owns DDA)
rules/thermo.smk       render_thermo · search · score
rules/bruker.smk       render · search_bruker · score_bruker   (lean v2 .d, imspy-free)
rules/sciex.smk        render_sciex · search_sciex · score_sciex   (lean ZenoTOF SWATH → open mzML)
```

## Backends

Pick one with `backend:` in `config.yaml` (or `--config backend=...`). Mirrors the flow's four `main()`
branches:

| backend | measurement | needs | phase 2 (search+score) |
| --- | --- | --- | --- |
| `v1` (default) | Bruker `.d` via the v1 `timsim` seam | `timsim_config` | — |
| `thermo` | Thermo `.raw` into a template | `template` | ✅ set `search_fasta` |
| `bruker` | lean v2 Bruker `.d` (imspy-free projector) | `reference_d` | ✅ set `search_fasta` |
| `sciex` | lean v2 SCIEX ZenoTOF SWATH → open mzML | — (synthesised SWATH) | ✅ set `search_fasta` |

Setting `search_fasta` on `thermo`/`bruker`/`sciex` appends DiaNN library-free search → score against
the render's co-emitted answer key — the number the whole simulate→search→score run exists to produce.
(A Bruker `.d` and an open mzML are DiaNN's native inputs on Linux; a Thermo `.raw` needs the .NET runtime.)

## Prerequisites

Same tool surface as the flow (see [`../README.md`](../README.md) / [`../Makefile`](../Makefile)):

- **Rust bins** (`timsim-proteome/digest/modify/precursors/design/yield/frag-input/spectra/
  render-thermo/render`) from `timsim-cli`. `config.yaml: timsim_bin` points at the dir holding them
  (default `target/release`).
- **Python prediction tools on `PATH`**: `timsim-ccs`, `timsim-rt`, `timsim-fragments`, the v1 `timsim`
  monolith, and `timsim_eval` (`make -C .. predict-deps`).
- **DiaNN** + the **.NET 8 runtime** for phase 2 (`config.yaml: diann`, `dotnet_root`). A Bruker `.d`
  is DiaNN's native input on Linux — no .NET; a Thermo `.raw` needs it.

## Run

```bash
# see the plan without running anything (this is what proves the structure dedup)
snakemake -n

# default v1 backend, 4 cores
snakemake --cores 4

# Thermo .raw + DiaNN search/score, override backend + external paths on the CLI
snakemake --cores 8 --config backend=thermo \
  template=/path/to/template.raw search_fasta=/path/to/db.fasta

# lean v2 Bruker .d, both samples
snakemake --cores 8 --config backend=bruker reference_d=/path/to/reference.d

# visualise the DAG
snakemake --dag --config backend=thermo template=/path/to/t.raw | dot -Tsvg > dag.svg
```

Per-rule `threads:` and `resources: mem_mb` are carried over from the flow's `threads=`/`ram=`
declarations (e.g. `simulate`/`render_*` reserve 8 GB, `search` 32 GB), so `--cores` / `--resources`
schedule the heavy nodes honestly — the flow learned that the hard way (a PTM-enriched design once
produced a 3.1 GB precursor table and the kernel OOM-killed two parallel renders).
```
