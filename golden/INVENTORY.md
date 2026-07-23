# Fixture inventory (dev-box scan, 2026-07-23)

Assets available to define a clean-slate **5K-protein** re-sim across three instruments. Paths are
dev-box-specific. Goal: pick ONE real anchor per instrument + the proteome, so the numbers mean something.

## A. Proteome sources (protein-level FASTA)

| proteins | what | path |
|---|---|---|
| 20,535 | **human UniProt + 172 contaminants** (canonical, real) | `hd02:data/fasta/human_20365_conts_172_..._validated.fasta` |
| 20,597 | human reference proteome UP000005640 | `hd02:data/tims2rescore/plasma_data/UP000005640_9606.fasta` |
| 2,500 | HeLa subset (what the DAG uses today) | `timsim-necro/flow/hela_subset.fasta` |
| 60 | current golden smoke slice | `golden/config/tiny.fasta` |

⚠️ `TIMSIM-HeLa5K-001.peptides.fasta` is **peptides** (`>pep_0`), a prior sim's output — NOT a proteome.

→ **5K proteins = a 5,000-entry slice of the real human UniProt FASTA.**

## B. Reference acquisitions (the anchor — real vs synthetic)

### Bruker timsTOF `.d`
| frames | grad | MS1/DIA | windows | real? | name |
|---|---|---|---|---|---|
| 21,473 | 41 min | 1023/20450 | 18,360 | **REAL** | `hd02:data/clusterQualityDev/G8027/raw.d` (dia-PASEF) |
| 11,217 | 21 min | 2244/8973 | 3,708 | **REAL** | `hd02:data/raw/synchro/synchro-hela.d` (synchro-PASEF HeLa) |
| 11,926 | 21 min | 746/11180 | 36 | **REAL** | `SUBMISSION/primitives/O240206_003_S1-B2_1_15479.d` |
| 34,155 | 60 min | 2010/32145 | 15,248 | SYNTH | `MIDIA-250K-NOISELESS-FRAG.d` (Slice-PASEF, current golden ref) |

### Thermo `.raw`
| size | what | path |
|---|---|---|
| 24 G | **REAL Astral, 60 min, DIA 3Da, 125ng CEA** | `timsim-rawbench-m0/templates/230724_AST_125ng_CEA_60mim_DIA3Da_7ms_Rep1.raw` |
| 21 G | **REAL Orbitrap HAP1, 60 min, DIA 2Th** | `timsim-rawbench-m0/templates/20231206_HAP1_1ug_60min_DIA_2Th_..._rep01.raw` |
| 3.5 G | REAL Astral nDIA 2Th, 15 min | `astral_ndia_2Th_15min.raw` |
| 664 M | current golden thermo template (smaller/trimmed) | `orbi_hela_dia.raw` |

### SCIEX ZenoTOF `.wiff`
| what | path |
|---|---|
| **REAL K562 ZenoTOF SWATH** (covid benchmark) | `PhantomBENCH/simulations/covid_sciex_native/*/…K562….wiff(2)` |

## C. Toolchain

| component | state |
|---|---|
| Rust render bins | `timsim-cli/target/release/`: `timsim-render` (Bruker v2), `timsim-render-thermo` (Thermo v2) |
| NECRO venv | necroflow 0.0.3, timsim-predict 0.1.0, pepdl 0.1.0, timsim-eval 0.1.0, koinapy 0.0.10 — **imspy-free** |
| DiaNN | 2.5.0 (Bruker `.d` native; Thermo `.raw` via .NET 8, present) |
| Sage | **not installed** (tool-axis sage slot can't run yet) |

## D. Decisions this surfaces

1. **Bruker anchor**: real `G8027/raw.d` (41 min dia-PASEF) or `synchro-hela.d` (21 min) — vs the synthetic
   MIDIA-250K we use now. A **real** anchor is the point of this exercise.
2. **Thermo anchor**: real `230724_AST…60min_DIA3Da` (Astral) or `HAP1 60min` (Orbitrap) — vs the trimmed
   `orbi_hela_dia.raw`. Note: 24 G/21 G templates are heavy to render against.
3. **Gradient**: pick ONE and match across the fixture where the anchors allow (~60 min if we want the big
   real templates; ~41 min if we anchor Bruker to G8027).
4. **SCIEX — writers DONE, lean wiring is the gap** (corrected): the native `.wiff` writer EXISTS and is
   validated — `sciexwiff` (github.com/theGreatHerrLebert/sciexwiff, `feature/wiffscan-codec`): pure-Rust
   `.wiff.scan` encode + GROW rebuild; native render = 241 prec / 217 prot = 84% of the mzML baseline
   (287/261); pwiz-validated; per-template profile for generalization. Plus `timsim-core/mzml.rs`
   (`render_db_to_mzml`) for open mzML. BUT both are wired to the **v1 world** — `render_db_to_mzml` reads a
   `synthetic_data.db`; `write_sciex_wiff` is a PyO3 method v1 `timsim` calls. Neither is driven by the lean
   v2 **parquet** feature space via a `timsim-cli` binary, so the DAG's `render_sciex` still uses v1 `timsim`
   (imspy). `sciexwiff` is **legal-held/private** (like `sciex-io`), so native `.wiff` must stay
   rustims-local; the open **mzML** path has no such constraint.
   → **DONE (mzML path):** `timsim-render-sciex` built + wired into the DAG (`--sciex`) — parquet feature
   space → synthesised SWATH schedule → open mzML via `timsim-core` (mzdata), no `.wiff`/`sciexwiff`.
   Verified fresh e2e through necro: 84.4% detectable recall, FDP 0.81% (on par with Thermo). Native
   `.wiff` output remains a rustims-local satellite reusing the validated `sciexwiff` writer (follow-on).
5. **Gradient render fix** (both Bruker + this re-sim): `timsim-render` currently defaults `--n-frames 3000`
   (5 min) instead of inheriting the reference's frame count — the bug behind the 26.6%. Fix before re-sim.
