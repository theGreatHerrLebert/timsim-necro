# v1 → v2 configuration mapping — Bruker dia-PASEF

> Companion to [`V1_V2_HEADTOHEAD.md`](V1_V2_HEADTOHEAD.md), phase 1 ("mapping table + independent
> sign-off, before any compute"). Written 2026-08-09 from the code, not from the names.
> **STATUS: DRAFT, awaiting independent sign-off.** Nothing here has been signed off; the (b) rows
> in particular are the ones a reviewer must check line by line.

## Scope

**Bruker dia-PASEF only** — the one mode v1's paper validated. Thermo, SCIEX, Waters and DDA are out
of scope and their settings are not mapped. Every claim below cites the file and line it came from.

The two surfaces compared are:

| | v1 | v2 |
|---|---|---|
| entry point | `timsim <config.toml>` (`imspy_simulation/timsim/simulator.py`) | `flow/timsim_flow.py:job` job TOML → `timsim-*` binaries |
| defaults | `get_default_settings()`, `simulator.py:286-575` | each binary's `clap` defaults + the flow's `SimpleNamespace`, `build_cfg`, `timsim_flow.py:1588-1620` |
| reference config | `/scratch/timsim-demo/v1-diapasef/config.toml` (+ its resolved `arguments-*.txt`) | `/scratch/timsim-demo/PhantomBENCH/simulations/v2/bruker_base.toml` |

**Which v2 flow is authoritative.** There are two copies of `timsim_flow.py` on this box and they are
NOT the same file. `/scratch/timsim-demo/timsim-necro/flow/timsim_flow.py` is 829 lines and still uses
the removed `necroflow.Rules` registry; `/scratch/timsim-demo/timsim-necro-repo/flow/timsim_flow.py`
is 1777 lines, is what `bruker_base.toml`'s `".pipeline"` key points at, and is the only one with the
Bruker-v2 render, the noise flags and the phospho/quant arms. **This document maps against
`timsim-necro-repo/flow/timsim_flow.py`.** Anyone benchmarking against the other file is benchmarking
a different program.

## Reading the classes

- **(a) exactly equivalent** — same semantics, same units, provably the same effect. Reserved for
  cases where the two implementations call the same primitive, or where the arithmetic is written
  out identically in both. If the argument is "the names match and it looks right", it is not (a).
- **(b) approximately equivalent** — comparable but differing, with the difference stated
  concretely (units, formula, magnitude). Every (b) row below names a number.
- **(c) not representable** — no counterpart in one direction or the other.

A row is (b) or (c) whenever equivalence could not be *shown*. Several v1 settings turned out to be
**dead** — present in the config surface and silently ignored by the code. Those are marked
`(c) — dead in v1` and matter mostly because the known-good v1 config sets them, which invites the
false impression that they are controlling something.

## Reference acquisition used for every number below

`/media/hd02/data/raw/dia/blanks/blanks-dia-PASEF/G241217_011_Slot2-2_1_16312.d`, read directly from
`analysis.tdf`:

| property | value |
|---|---|
| frames | 17,646 (1,103 MS1 + 16,543 MS2) → 16-frame cycle, 1.687 s |
| frame period | 1861.3 s / 17,646 = **0.10548 s/frame** |
| mobility scans | 927 |
| `1/K0` acquisition range | 0.600 – 1.600 (→ ≈ 0.00108 `1/K0` per scan, linearised) |
| m/z range | 99.99 – 1700.0, `DigitizerNumSamples` 396,224 |
| DIA windows | 36 (12 window groups × 3), isolation width 25 Th |
| **window-table collision energy** | **20.00 – 58.12 eV** |
| gradient | 1861.3 s (v1 config declares `gradient_length = 1860.0`) |

---

# The table

## 1. Inputs, geometry, acquisition

| v1 setting | v1 default/meaning | v2 counterpart | class | notes |
|---|---|---|---|---|
| `fasta_path` | one FASTA; organism recovered by substring-matching the header | `proteome_spec` → `[[source]] path/organism` (`spec.rs:24-37`) | **(a)** | Same protein universe when one source is declared. v2 additionally *declares* the organism instead of inferring it (`spec.rs:18-23`); irrelevant for a single-organism Bruker run. |
| `reference_path` + `use_reference_layout = true` | the real `.d` whose frame schedule, DIA windows and TOF/mobility calibration are replayed | `bruker_reference` → `timsim-render --reference-d` (`timsim_flow.py:679-685`); `DiaSchedule::from_reference` (`dia.rs:33-37`) copies `DiaFrameMsMsInfo` + `DiaFrameMsMsWindows` verbatim | **(a)** | v2 has no separate `use_reference_layout` switch — setting `bruker_reference` *is* the switch, and it selects the lean Bruker pipeline (`timsim_flow.py:1245-1251`). |
| `acquisition_type = "DIA"` | DIA vs DDA | `--dia`, hard-wired into `_RENDER_HEAD` (`timsim_flow.py:683`) | **(a)** | |
| `apply_fragmentation = true` | render MS2 | implicit — the DIA render always projects MS2 | **(a)** | |
| `gradient_length = 1860.0` | seconds; v1's RT model output is scaled onto this, then mapped to the reference's frame times | none. `--n-frames 0` inherits the reference frame count (`render.rs:391-398`); the RT index is linearly stretched onto `[0, n_frames-1]` (`render.rs:490`) | **(b)** | Same *effective* span here (1861.3 s reference vs 1860.0 declared, 0.07 % apart), but the mechanism differs: v1 anchors RT in **seconds** and can disagree with the reference; v2 anchors it in **frames** and cannot. Setting `gradient_length` to anything other than the reference's own length is a v1-only capability. |
| `reference_in_memory`, `num_threads`, `batch_size`, `use_gpu`, `gpu_memory_limit_gb`, `lazy_frame_assembly`, `frame_batch_size`, `silent_mode`, `log_level` | execution/resource knobs | necroflow `threads=`/`ram=` declarations (`timsim_flow.py:694`, `:601-605`), `--render-chunks`, `--no-parallel` | **(c)** | No model effect in either tool. Listed so the benchmark manifest can record them as *cost* variables, not accuracy variables. |
| `use_bruker_sdk = true` | write the `.d` via the Bruker SDK when available | none — v2 always writes through `ms-io`'s `TdfWriter` | **(c)** | Bounded: any writer-level difference in the produced `.d` is a v1-vs-v2 difference the benchmark cannot control for. |
| `emit_provenance = true`, `provenance_embed = true`, `provenance_key_path` | Ed25519-signed mzPROV self-disclosure, embedded in `analysis.tdf` | none — no `provenance`/`mzprov` symbol in `timsim-cli/src/bin/render.rs` | **(c)** | v1's default is ON, so the v1 arm's `.d` carries an extra table the v2 arm's does not. Confirm DiaNN ignores it (it should) or disable it for the benchmark. |

## 2. Digestion and the peptide space

| v1 setting | v1 default/meaning | v2 counterpart | class | notes |
|---|---|---|---|---|
| `cleave_at = "KR"`, `restrict = "P"` | trypsin, no cleavage before proline | `timsim-digest --enzyme` (default `trypsin`) → `Enzyme::new("KR", "P", true)` (`timsim-chem/src/enzyme.rs:122`) | **(a)** | Literally the same two strings. |
| `missed_cleavages = 2` | max missed cleavages enumerated | `--max-missed-cleavages 2` (`digest.rs:28-29`; flow `timsim_flow.py:1591`) | **(a)** | Same enumeration bound. What each tool *does* with the resulting peptides differs — see `digestion_efficiency` in §6. |
| `min_len = 7`, `max_len = 30` | peptide length filter | `--min-length 7 --max-length 30` (`digest.rs:30-33`; flow `:1592-1593`) | **(a)** | v2's own default is `max_length = 50`; the flow overrides it to 30, so the pair matches only because the flow says so. |
| `digest_proteins = true` | digest rather than read peptides | implicit | **(a)** | |
| `decoys = false` | append decoy proteins | none in `timsim-digest`/`timsim-proteome` | **(c)** | Off in the reference v1 config, so it costs the benchmark nothing — but v2 cannot reproduce a v1 decoy run. |
| `remove_degenerate_peptides = false` | drop peptides mapping to >1 protein | none; v2 keeps shared peptides as multiple *occurrences* by design (`digest.rs` report line, `timsim-cli/src/bin/digest.rs:206-207`) | **(c)** | Off in the reference config. v2's model is the stronger one here (occurrences are first-class), but it is not the same knob. |
| `n_proteins = 20000` | **sample** this many proteins from the FASTA; the rest do not exist (`simulate_proteins.py:228-238`) | `[design].n_proteins` (`spec.rs:60-65`) — a *presence* subset; excluded proteins stay in the structure with `amount_amol = 0` (`design.rs:918-934`) | **(b)** | Different objects. v1 removes proteins from the universe, so protein-level FDR is unanswerable for them; v2 keeps them at amount 0, so it is answerable. For a benchmark on a FASTA with ≤ `n_proteins` entries the two coincide (v1 clamps, `simulate_proteins.py:228-229`) — **use a FASTA no larger than `n_proteins` and this row collapses to (a)**. |
| `num_sample_peptides = 10000`, `sample_peptides = true`, `sample_seed = 41` | after abundance is assigned, `peptides.sample(n, random_state=sample_seed)` (`simulator.py:2051`) | `--max-peptides` (`digest.rs:40-49`; flow `--max-peptides`, `timsim_flow.py:1518`) + `[design].load_ng` | **(b)** | Both are a seeded uniform subsample of the digest, so as a *cap* they are comparable. The difference is what they mean: v2's `spec.rs:60-63` states the design position — "There is deliberately no peptide-count knob: peptides vanish by falling below the detection limit, not at random" — and `--max-peptides` is documented as a tractability cap, not a model parameter. v1's cap **is** the model: there is no detection limit, so `num_sample_peptides` is the only thing setting how many peptides exist. Consequence for the benchmark: matching the two by setting `--max-peptides = num_sample_peptides` compares v1's *whole sample model* against v2's *debug cap*, and silences v2's `load_ng` detection-limit behaviour entirely. **This is the single most distorting (b) row.** |
| `num_peptides_total = 250000` | ceiling on peptides carried out of the digest (`simulator.py:2013`) | none | **(c)** | Inert when it exceeds the digest size, which it does for the tiny/5 k fixtures. Verify per run. |
| `upscale_factor = 100000` | multiplicative scale on protein `events` (`simulate_proteins.py:176`) | none — v2's absolute scale comes from `[design].load_ng` mass balance (`design.rs:906, 923-936`) | **(b)** | v1's scale is an arbitrary integer multiplier; v2's is nanograms on column, converted through the molar-weighted mean MW of the present proteins. Not convertible without knowing the proteome; see §6. |
| `modifications = ""` | no modifications | `mods` spec is **mandatory** (`timsim-modify --mods`, `timsim_flow.py:449-452`) | **(b)** | v2 has no "no spec" mode. Closest match is a spec declaring only `Carbamidomethyl` at occupancy ~1.0 — but v1 applies carbamidomethyl as a **fixed** static mod on C (`simulate_proteins.py:193`, `static_mods = {"C": "[UNIMOD:4]"}`), i.e. occupancy exactly 1.0, whereas the PhantomBENCH spec uses `occupancy = 0.98` (`mods_covid.toml`). Set `occupancy = 1.0` for the benchmark and this becomes (a). |
| `phospho_mode = false` | phospho-specific handling | v2 mods spec + `--phospho` arm | **(c)** | Off; out of scope for this benchmark. |
| `proteome_mix = false`, `multi_fasta_dilution` | multi-FASTA mixing by dilution factor | `[[source]]` per organism + `[[condition]].mix` (`spec.rs:100-113`) | **(b)** | Off in the reference config. v2's model is richer (declared organism + mass fractions + `"rest"`); the two do not agree on what the number means (v1: a multiplicative dilution factor; v2: a mass fraction of `load_ng`). |

## 3. The ion layer: charge and isotopes

| v1 setting | v1 default/meaning | v2 counterpart | class | notes |
|---|---|---|---|---|
| `binomial_charge_model = true` + `p_charge` (default 0.8; **reference config uses 0.5**) | `Binomial(n = protonatable sites, p = p_charge).pmf(z)`; sites = N-term + every H/R/K (`mscore-0.5.0/src/algorithm/peptide.rs:382-391, 413-428`) | `timsim-precursors --charge-model binomial --charged-probability <p>` (`precursors.rs:64-68`, `timsim-chem/src/ionize.rs:113-114`) | **(a)** *for the PMF* | v2's Poisson-binomial convolution (`ionize.rs:256-302`) reduces exactly to the binomial when all site probabilities are equal (`ionize.rs:241-242`). **The flow does not pass `--charge-model`, so the default `site-specific` is what runs today** (`timsim_flow.py:1598`, `precursors.rs:64`). Requires an explicit override to reach (a) — see §7. |
| (as above, with the flow default) | — | `--charge-model site-specific`: N-term 0.93, R 0.97, K 0.95, H 0.80 (`ionize.rs:141-149`) | **(b)** | Same site *set* as v1, different per-residue probabilities. A clean tryptic peptide comes out 88.7 % 2+ under site-specific; under `p = 0.5` binomial with 3 sites the mode is 1+/2+. This changes the precursor m/z population, so it changes what a DIA window transmits. Not a small effect. |
| `max_charge = 4` | truncate at z ≤ 4 | `--max-charge 4` (`precursors.rs:69-70`) | **(a)** | Same bound. What happens to the discarded mass is *not* the same — next row. |
| `normalize_charge_states = true` | divide surviving z ≥ 1 states by their own partial sum (`imspy-predictors/.../ionization/predictors.py:121-126`), so the neutral fraction, the sub-threshold states and everything above `max_charge` are all redistributed into the survivors | v2 divides out z = 0 only (`ionize.rs:287`); z > `max_charge` is **deleted and the loss reported** (`ionize.rs:281-301`), so the fractions deliberately sum to < 1 | **(b)** | Concrete: for a peptide with 5 protonatable sites at `p = 0.5`, `P(z>4) = 0.031`. v1 inflates every surviving charge state by 1/(1 − 0.031) ≈ 1.032; v2 leaves the 3.1 % missing. Multiplicative and peptide-dependent, so it is an abundance bias correlated with peptide length. |
| `min_charge_contrib = 0.005` (**reference config uses 0.25**) | drop charge states below this *post*-normalisation share (`predictors.py:127`) | none found in `precursors.rs` / `ionize.rs` | **(c)** | At the reference config's 0.25 this is aggressive: it keeps only charge states carrying ≥ 25 % of the ionised mass, i.e. typically one or two per peptide. v2 keeps every state with fraction > 0. **This alone changes the number of precursors rendered.** |
| `charge_state_one_probability = 0.0` | additive bump to z = 1 | none | **(c) — inert in v1** | Only consumed on the deep-charge-model branch (`simulate_charge_states.py:56`), which `binomial_charge_model = true` does not take. No effect either way. |
| `isotope_k = 8`, `isotope_min_intensity = 1`, `isotope_centroid = true` | intended as the averagine envelope depth/floor | `--isotope-depth` (default 6, `precursors.rs:71-73`) | **(c) — dead in v1** | `simulate_precursor_spectra_averagine`, their only consumer (`simulate_precursor_spectra.py:40-60`), **has no caller**. The live path is `simulate_precursor_spectra_sequence` (`simulator.py:2319-2323`) → `mscore::algorithm::isotope::generate_precursor_spectrum` → `calculate_isotopic_spectrum(1e-3, 1e-9, 200, 1e-6)` (`mscore-0.5.0/src/algorithm/isotope.rs:582-589`). |
| *(the live v1 isotope envelope)* | exact-composition envelope, abundance threshold 1e-9, up to **200** peaks | exact-composition envelope, fixed depth **6** (`precursors.rs:72-73`, `:145`) | **(b)** | Both are exact from elemental composition — no averagine on either side. The truncation rule differs: v1 by abundance (1e-9) capped at 200 peaks, v2 by *count* (6). For a 1500 Da tryptic peptide the peaks beyond M+5 are ≪ 1 % and the difference is immaterial; for the long/heavy tail it is not. Raise `--isotope-depth` if the benchmark includes long peptides. |

## 4. Retention time, mobility and peak shape

| v1 setting | v1 default/meaning | v2 counterpart | class | notes |
|---|---|---|---|---|
| `rt_model = None` | imspy's own `DeepChromatographyApex` GRU, scaled onto `gradient_length` (`simulate_retention_time.py:60-70`) | `timsim-rt` → **Chronologer** by default (`timsim-predict/src/timsim_predict/rt.py:3`), emitting a portable RT *index* with `index_min`/`index_max` stamped in the artifact metadata | **(b)** | **Two different trained models.** Nothing in either codebase claims they agree. Elution *order* is broadly conserved between RT predictors, but per-peptide apexes are not, so co-elution structure — which is exactly what a DIA search must disentangle — differs. There is no flag that makes them the same model. |
| `sigma_lower_rt = None`, `sigma_upper_rt = None`, `sigma_alpha_rt = 4`, `sigma_beta_rt = 4` | EMG σ, sampled per peptide from a scaled Beta(4,4) over `[σ_mid·0.75, σ_mid·1.25]` with `σ_mid = gradient/3600·0.75 + 1.125` **seconds** (`simulate_frame_distributions_emg.py:11-29, 279-288`) | `timsim-render --sigma-frames` (default **30 frames**, `render.rs:82-83`), a symmetric Gaussian | **(b)** | Numbers: at `gradient_length = 1860`, v1 draws σ ∈ [1.134, 1.891] s, mean **1.513 s**. v2's 30 frames × 0.10548 s = **3.164 s**. **v2's chromatographic peaks are 2.09× wider in σ**, and symmetric where v1's are right-tailed. The flow does not pass `--sigma-frames`, so the default is what runs. |
| `k_lower_rt = 0`, `k_upper_rt = 10`, `k_alpha_rt = 1`, `k_beta_rt = 20` | EMG tail: `k ~ 0 + Beta(1,20)·10`, `λ = 1/(kσ)` (`simulate_frame_distributions_emg.py:279-288`) | `--peak-shape emg` (**the default**) with `--emg-k` defaulting to `V1_DEFAULT_EMG_K = 10/21` (`timsim-cli/src/bin/render.rs:98-106`, `render.rs:125`) | **(b)** | **Was (c) — closed 2026-08-09.** v2 now renders v1's EMG: the same `λ = 1/(kσ)` kernel, defaulting to v1's *mean* draw `E[k] = 10/21 = 0.476` (τ = kσ ≈ 0.72 s). v2's inversion is independent of v1's (closed-form CDF + survival-function bisection vs v1's 1000-step Riemann sum), and `tests/emg_v1_parity.rs` checks them against each other. **What remains (b): v1 SAMPLES `k` per peptide from `Beta(1,20)·10`; v2 applies the mean to every peptide.** So the population mean tail matches and the per-peptide spread does not — v1's `k` ranges over roughly `[0, 1.6]` at the 1st/99th percentile, v2's is a point mass at 0.476. The flow passes neither flag, so the default is what runs. |
| `target_p = 0.999`, `sampling_step_size = 0.0001`, `n_steps = 1000`, `remove_epsilon = 1e-4` | where the EMG profile is truncated and how finely it is integrated | `--n-sigma` (default 3.0, `render.rs:86-87`) | **(b)** | Different truncation rules: v1 truncates at a cumulative-probability target (0.999) with an ε floor; v2 truncates at a fixed ±3σ (which for a Gaussian is p = 0.9973). Comparable in intent, ~0.2 pp apart in tail mass, but v1's ε floor additionally deletes low-weight frames. |
| `min_rt_percent = 2.0`, `exclude_accumulated_gradient_start = true` | void-volume correction: peptides the RT model dumps into the first 2 % of the gradient are thinned to the median bin density, identity-keyed by sequence hash (`simulator.py:1944-1974`) | none | **(c)** | v1 default is ON and it **deletes peptides**. A v1 run and a v2 run therefore do not render the same peptide set even from an identical digest. Must be reported as a v1-only filter, or disabled in v1 (`exclude_accumulated_gradient_start = false`) for the matched arm. |
| `ccs_model = None` | imspy GRU CCS + per-ion CCS std | `timsim-ccs` → `pepdl.ccs.DeepPeptideIonMobilityApex` (`timsim-predict/src/timsim_predict/ccs.py:88`), CCS in Å² plus `ccs_std` | **(a)** | Same model, extracted: the DL-extraction work records CCS parity as **exact on 40,509 precursors** between `pepdl` and `imspy`. Both then convert CCS → `1/K0` by Mason–Schamp with the same constants (N₂, 28.013 u, 31.85 °C — `simulator.py:335-336`; `mobility_ce.rs:49-52`; `render.rs:1035`). |
| `use_inverse_mobility_std_mean = true`, `inverse_mobility_std_mean = 0.009` | keep the deep model's per-ion CCS-std *shape*, rescale the population **mean** to 0.009 `1/K0` (`simulator.py:2269-2280`); Gaussian in `1/K0` | `--sigma-scans` (default **4.0 scans**, `render.rs:84-85`), one width for **every** ion | **(b)** | Numbers: 0.00108 `1/K0` per scan on this reference → v1's mean σ ≈ **8.33 scans**; v2's is **4.0 scans** ≈ 0.00432 `1/K0`. **v2's mobility peaks are 2.08× narrower, and v2 has no per-ion width at all** — it throws away the `ccs_std` column `timsim-ccs` produces. (Scan↔`1/K0` linearised; the real map is Bruker ModelType-2 and mildly nonlinear.) |
| `projection_mode = "off"` | opt-in Rust re-projection of the frame/scan distributions | n/a | **(c)** | Off; the legacy kernels are what run. No v2 concept. |
| `re_scale_rt = false`, `rt_variation_std`, `ion_mobility_variation_std`, `intensity_variation_std` (all `None`) | per-replicate jitter on RT / mobility / intensity | `[variance] biological`, `biological_heterogeneity`, `technical` (`spec.rs:259-272`, `design.rs:1040-1073`) | **(b)** | Different axes. v1 jitters the *measurement* (RT, mobility, intensity) per run; v2 jitters the *amounts* per biological replicate, identity-keyed per protein, and gives technical replicates identical amounts by construction (`design.rs:1112-1114`). All off/absent in the reference v1 config, so the benchmark's single-replicate arm is unaffected — but a replicate-level comparison is **not** supported by this mapping. |

## 5. Fragmentation

| v1 setting | v1 default/meaning | v2 counterpart | class | notes |
|---|---|---|---|---|
| **collision energy (DIA)** | not a config key — v1 reads the **per-ion** CE off the reference `.d`: `(frame → window group, scan) → CollisionEnergy` from `DiaFrameMsMsWindows` (`timsim-core/src/handle.rs:1079-1109`, `ion_map_fn_dia` at `:1141-1183`). An ion transmitted in several window groups gets **several** CEs and several predicted spectra. On this reference that is **20.00 – 58.12 eV**. | `timsim-fragments --collision-energy` — **one scalar for the whole run**, `25.0` in the flow (`timsim_flow.py:601-605`, `:1539`) | **(b)** | Magnitude: the flow's flat 25.0 eV is **−33.1 to +5.0 eV** away from the CE v1 would have used for the same ion. And v2 still *writes the reference's own 20–58 eV window table into the output `.d`* (`dia.rs:26`, `:33-37`) — so the produced file declares a CE the fragment intensities were not predicted at. See the box below; this is the row that most constrains the benchmark. |
| *(same, with the opt-in per-precursor CE)* | as above | `timsim-frag-ce` (`timsim-cli/src/bin/frag_ce.rs`) → `timsim-fragments --collision-energies` (`imspy_simulation/timsim/jobs/fragments.py:404-410`) | **(b)** | Still (b), for two independent reasons — see the box. |
| `round_collision_energy = true`, `collision_energy_decimals = 0` | round the CE before prediction | none | **(c)** | v1 rounds to integer eV; v2 passes the float through. Sub-eV, but it means the two never feed the predictor bit-identical inputs even if the CE model is matched. |
| `intensity_model = None` / `fragment_intensity_model = None` | the bundled PROSPECT-fine-tuned timsTOF PyTorch model | `--frag-model ""` = local timsTOF (`timsim_flow.py:1538`, `:606-609`) | **(a)** *if the same `timsim-fragments` is on `PATH`* | Two packages ship a `timsim-fragments` console script: `imspy-simulation` (`pyproject.toml:42` → `imspy_simulation.timsim.jobs.fragments:main`) and `timsim-predict` (`pyproject.toml:23`). **Only the imspy one implements `--collision-energies`.** Which one the flow invokes depends on the active venv — pin it in the manifest. |
| `down_sample_factor = 0.5` | zero `int(nnz × 0.5)` of each spectrum's non-zero predicted fragments, sampled without replacement with probability ∝ `1/(value·10⁴)` — weakest fragments die first (`utility.py:481-538`) | `--floor` (default `1e-3`, `fragments.py`) — a deterministic relative-intensity cutoff | **(b)** | Different operations. v1 removes **half the peaks, stochastically, intensity-weighted**; v2 removes **all peaks below 0.1 % of the base peak, deterministically**. On a typical Prosit spectrum the 1e-3 floor removes far fewer than half the non-zero slots, so **v2's MS2 spectra carry more peaks than v1's**. Set `down_sample_factor = 0.0` in v1 for the matched arm and the row becomes near-(a). |
| `precursor_survival_min = 0.0`, `precursor_survival_max = 0.0` | fraction of precursor surviving the quad intact and bleeding into MS2 | `--precursor-survival-min/--precursor-survival-max` (`render.rs:109-116`) | **(a)** | Same knob, same default, identity-keyed per ion in v2. Off in the reference config. |
| `quad_isotope_transmission_mode = "none"`, `quad_transmission_min_probability`, `quad_transmission_max_isotopes` | quad-dependent isotope transmission | none in `timsim-cli/src/` | **(c)** | Set to `"none"` in the reference config, so no effect on this benchmark; a v1 run with `"per_fragment"` has no v2 counterpart. |

> ### Collision energy is **not** (a). Both reasons, with numbers.
>
> The head-to-head plan recorded CE as "now (a) — 0.0 max absolute difference over 11,164,886 rows /
> all 893 reachable scans". **I could not reproduce or locate that measurement**; nothing in
> `/scratch/timsim-demo` contains those figures. What does exist is
> `mobility_ce.rs:178-188`, which pins `collision_energy_at` against the closed form
> `CE = 54.1984 − 0.0345·scan` over scans 0..2000 with `assert_eq!(worst, 0.0)`. That is a real
> 0.0-difference result, but it is a *self*-consistency test of the ramp, not a v1↔v2 comparison.
>
> **Reason 1 — the flag is not wired.** `timsim-frag-ce` exists and `--collision-energies` exists,
> but no caller connects them. The flow's `fragments` node is
> `"timsim-fragments --precursors {…} --collision-energy {collision_energy} --model {frag_model} --out {…}"`
> (`timsim_flow.py:601-603`) — a scalar, and there is **no `frag_ce` node anywhere in
> `timsim_flow.py`** (grep for `frag_ce` over the flow, backend and GUI returns nothing). So today's
> v2 Bruker DIA render predicts every precursor at 25.0 eV.
>
> **Reason 2 — even wired, it is the wrong CE model for DIA.** `timsim-frag-ce` implements the
> **dda-PASEF ramp** (`pasef_policy(CE_BIAS = 54.1984, CE_SLOPE = −0.0345)`, `frag_ce.rs:66-71`),
> reusing v1's `dda_selection_scheme` constants. But v1's **DIA** path does not use that ramp at all —
> it reads the reference `.d`'s window-group CE table. Evaluating the ramp at the mid-scan of each of
> this reference's 36 windows gives **24.94 – 43.28 eV** against the table's **20.00 – 58.12 eV**:
> a per-window difference of **−14.85 to +9.72 eV** (mean +0.37 eV). Matching the DDA ramp to the DIA
> table is not a wiring job; it is a different derivation.
>
> A third, smaller mismatch: v1 emits one `(peptide, ion, sequence, charge, CE)` row **per distinct CE
> the ion was transmitted at** (`handle.rs:1171-1177`), so a precursor spanning two window groups gets
> two spectra. v2's model is one CE per `(sequence, charge)` by construction — that invariant is what
> the dedup rests on (`mobility_ce.rs:29-36`), and `timsim-fragments` *refuses the run* if it is
> violated (`fragments.py:411-413`). The two models are not merely differently parameterised; they
> have different cardinality.

## 6. Abundance

| v1 setting | v1 default/meaning | v2 counterpart | class | notes |
|---|---|---|---|---|
| `intensity_mean = 5`, `intensity_min = 4`, `intensity_max = 8`, `intensity_value = 6`, `sample_occurrences = true` | present in the reference config's `[peptide_intensity]` block | `[abundance]` in the design spec | **(c) — dead in v1** | All five are in `_LEGACY_IGNORED_KEYS` and are **deleted** during config translation (`simulator.py:232-236, 279-281`). Their only consumer, `simulate_peptide_occurrences` (`jobs/simulate_occurrences.py:5`), has no caller. They control nothing. |
| *(the live v1 protein abundance)* | uniform rank draw `np.random.uniform(1, 1e4)` indexed into a fixed hockey-stick curve `get_tenzer_hokey()`, × `upscale_factor` (`simulate_proteins.py:155-179`) — **draw-order-keyed**, assigned positionally at `:177` | `[abundance.<org>] source = "hockeystick"` → `hockey_stick(uniform_01(hash("rank", protein_id, seed)), decay = 0.06, tail = 1e-4)` (`design.rs:779-792`, `:111-113`; defaults `spec.rs:91-96`) | **(b)** | Same curve *family*, different keying and different normalisation. v1: a protein's abundance depends on its position in a bulk RNG draw, so adding a protein reshuffles every subsequent one (`simulator.py:1655-1658` says so). v2: keyed on the accession (`design.rs:412-414`), stable under any change to the protein set. Then v1 scales by `upscale_factor` while v2 renormalises to `load_ng` mass balance (`design.rs:997-1011`). The **marginal distributions are comparable; the per-protein assignments are not**, so a protein-level abundance-vs-recall comparison across the two tools is meaningless unless it is done on the marginal only. |
| *(flyability)* | `10 ** N(−2, 1)` truncated to `[1e-4, 1]`, drawn in bulk by row order and **multiplied into `events`, then cast to int32** (`simulate_peptides.py:47-60, 132-138, 156`) | `--flyability lognormal --flyability-median 1e-2 --flyability-sigma 1.0` → same distribution, keyed on `hash("flyability", sequence, attempt)`, kept as its **own column** (`ionize.rs:317-344, 371-406, 429-432`) | **(b)** | **Same distribution, different keying and different placement.** v1 conflates abundance × flyability into one truncated integer; v2 keeps `ion_amol = peptide_amol × modform_fraction × ionization_propensity × charge_fraction` factorised (`ionize.rs:421`). v1's int32 cast quantises away everything below 1 event — a hard, undeclared detection floor v2 does not have. |
| (no v1 setting) | v1 has no digestion-yield model; peptides are subsampled per protein by `np.random.choice` (`simulate_peptides.py:11-44`) and each inherits the protein's `events` **unchanged** (`:112`) | `timsim-yield --digestion-efficiency 0.9` (`yield.rs:48-53`, `:130`), an exact expectation `p_yield = b_i·b_j·Π(1 − p_eff(k))` (`timsim-chem/src/digest.rs:26-39, 489-491, 532-567`) | **(c)** | v2 models per-peptide yield; v1 does not. Every peptide of a protein has the same amount in v1. **This is a first-order difference in the intensity distribution and there is no v1 knob that reproduces it.** |
| `intensity_multiplier = 1.0` | multiplier on `events` — **only on the `--from-findings` path** (`simulator.py:1841`, `load_findings.py:504, 565`) | — | **(c) — inert on this path** | Not used on the FASTA path at all. |
| (no v1 setting) | — | `--intensity-scale` (`5.0e5` in the flow, `timsim_flow.py:1540`) — the detector-count quantum applied at quantisation (`render.rs:550, 689`) | **(c)** | v2-only. v1's absolute detector scale falls out of `events × relative_abundance`; there is no single knob. **Consequence: the two arms' absolute intensity scales are not commensurable**, and any threshold-crossing metric (`--min-peak-intensity`, DiaNN's own intensity handling) is affected differently. Calibrate `--intensity-scale` against a v1 run before the benchmark, and record the calibration. |
| `findings_reference_median`, `from_findings`, `findings_path`, `from_existing`, `existing_path` | alternative input seams | `--spike-into` is the nearest thing, and it is not the same | **(c)** | Off; out of scope. |

## 7. Noise

All v2 A1 noise reuses v1's own mscore parameterisation, so this block is the closest to (a) anywhere
in the mapping — but not all the way.

| v1 setting | v1 default/meaning | v2 counterpart | class | notes |
|---|---|---|---|---|
| `mz_noise_precursor` (bool) + `precursor_noise_ppm = 5.0` | `mz' ~ N(mz, (mz·ppm/1e6 / 3)²)` — **ppm is the 3σ envelope**, σ = 1.67 ppm at ppm = 5 (`mscore-0.5.0/src/data/spectrum.rs:315-322`) | `--noise-mz-ppm` → `m·(1 + N(0,1)·ppm·1e-6/3)` (`render.rs:1364-1379`); `0` = off | **(b)** | **The 3σ convention matches exactly** — same divisor, same primitive. Three differences remain: (i) v1 **redraws per scan contribution**, v2 draws **one offset per `(precursor_id, peak_index)` constant across the whole elution** (v2 documents this at `render.rs:1000-1004`) — marginal mass-error distribution identical, within-elution correlation opposite; (ii) v1 uses unseeded `rand::thread_rng()`, so **a v1 noise run is not reproducible**, while v2 is seeded (`--noise-seed`); (iii) v1 re-bins to 1e-6 m/z after noising (`spectrum.rs:324-337`), v2 goes straight to TOF. |
| `mz_noise_fragment` + `fragment_noise_ppm = 5.0` | same draw on MS2 peaks and on surviving-precursor bleed-through (`timsim-core/src/dia.rs:731-737, 776-782`) | `--noise-frag-ppm`, selected by `is_frag = (ms_level == 2)` (`render.rs:1364`, `:1562`, `:1682`) | **(b)** | Same three differences as the row above; the "precursor bleed-through uses the *fragment* ppm" convention matches. |
| `mz_noise_uniform = false` | when true, and on the Bruker DIA path where `right_drag` defaults to **True** (`builders/dia.py:166`): `U[mz·(1 − ppm/6e6), mz·(1 + ppm/2e6)]` — **asymmetric, mean +ppm/6** | `--noise-mz-uniform` → symmetric `U[−ppm, +ppm]·1e-6` (`render.rs:1370-1371`) | **(b)** | v2 implements v1's `right_drag = false` branch, which v1's DIA path never takes. Against v1's actual default the v2 uniform is **3× wider in total width and has no +0.83 ppm mean shift**. `render.rs:150` concedes it: "v1's asymmetric `right_drag` tailing variant is not ported." Both configs set this **false**, so it does not affect the benchmark — but it is not a substitute. |
| `add_real_data_noise = false` | sample real background peaks from the reference `.d` and add them to each rendered frame (`jobs/add_noise_from_real_data.py:101-123`) | `--noise-real-data` → `build_frame_noise` (`render.rs:768-853`, `:1341-1351`) | **(b)** | Same construction (per output frame, sample N reference frames of the matching class, filter, Bernoulli-downsample, deposit additively). Two divergences: (i) v1 **retries** when a sampled reference frame is empty, up to `(nf·8).max(16)` attempts (`ms-io-0.2.0/src/data/dia.rs:684-698`); v2 takes exactly `nf` picks with no retry, so v2 deposits slightly **less** background on a sparse reference; (ii) v1 merges background into signal in **m/z space at 1e-6** and re-derives TOF (`mscore-0.5.0/src/timstof/frame.rs:749-802` + `tdf.py:239`), v2 deposits the reference's **raw TOF index** (`render.rs:682-748`) — so a background peak can land ±1 TOF bin apart between the two. |
| `num_precursor_noise_frames = 5`, `num_fragment_noise_frames = 5` | reference frames sampled per output MS1 / MS2 frame | `--noise-precursor-frames`, `--noise-fragment-frames` (defaults 5/5, `render.rs:166-171`) | **(a)** | Same pools (MS1 by `ms_ms_type == 0`; MS2 by matching window group), same with-replacement draw, same defaults. |
| `precursor_sample_fraction = 0.2`, `fragment_sample_fraction = 0.2` | keep each sampled background peak with probability `p` (`mscore-0.5.0/src/timstof/frame.rs:584-605`) | `--noise-precursor-fraction`, `--noise-fragment-fraction` (0.2/0.2, `render.rs:176-181`) | **(a)** | Bernoulli(p) per peak in both; `<=` vs `<` on a continuous draw is immaterial. v2's draw is seeded. |
| `reference_noise_intensity_max = 9999999.9` (v1 code default 30) | keep reference peaks with intensity in `[1, max]`; **`scan_max` is hardcoded to 1000** (`ms-io-0.2.0/src/data/dia.rs:693-694`) | `--noise-intensity-max` (default 150000, `render.rs:172-175`); the scan bound is the run's actual `n_scans` (`render.rs:824-826, 838-840`) | **(a)** *on this reference* | The predicates are otherwise identical. v1's hardcoded `scan_max = 1000` would drop peaks on a ≥ 1000-scan run — this reference has **927** scans, so the two agree here. **On a different reference `.d` this row becomes (b).** Defaults differ by 5000× (30 vs 150000), so the value must be set explicitly in both. |
| `noise_frame_abundance = false`, `noise_scan_abundance = false` | per-frame / per-scan multiplicative jitter `signal·(1+U(0, nl))`, `nl ~ U(0,2)`, renormalised (`utility.py:57-61`) | `--ion-count-noise` — **Poisson counting statistics** per detector bin | **(b)**, deliberate divergence | **Was (c) — closed 2026-08-11, but NOT by porting v1.** v1's jitter is multiplicative, so a bright and a dim peak get the same RELATIVE wobble; real counting does the opposite. v2 draws `Poisson(expected counts)` instead — variance = mean, tested across `lambda` 0.5–5000. Applies to synthetic signal only (A2 background is real counts that already carry shot noise). Measured: peaks/scan 9.96→11.96 MS1, 5.34→6.92 MS2, intensity/peak 1078→898 and 80.4→57.6, TIC preserved to +0.5%. v1's variant is off in its own reference config, so nothing is lost by not matching it — but it is a divergence, not parity. |
| `rt_variation_std`, `ion_mobility_variation_std`, `intensity_variation_std` (all `None`) | per-run jitter on RT / mobility / intensity | `--run-rt-sd` (s), `--run-im-sd` (1/K0), `--run-intensity-cv` | **(b)** | **Was absent — added 2026-08-11.** Without it two technical replicates of one design were BIT-IDENTICAL. Keyed on `(sample, precursor)` so a run is reproducible while a different sample gets a different realisation — v2's identity-keying discipline applied to the measurement axis, where v1 draws from a global stream. Intensity is lognormal, not additive: injection variation scales with the amount and an additive term drives small amounts negative on a 17-order abundance axis. Verified unbiased (mean factor 1.000 at the requested CV). |
| `superimpose_on_reference = false` | overlay the synthetic run on the real reference `.d` (`add_noise_from_real_data.py:128-152`) | `--spike-into` (`render.rs:188-195`, `:862-920`) | **(b)** | Same arithmetic (`real + synthetic`), different failure mode: v1 silently passes through frames beyond the reference's frame count; v2 validates contiguity and MS-type agreement and **errors out** (`render.rs:401-406`). v1's `elif` precedence over `add_real_data_noise` vs v2's mutual-exclusivity error (`render.rs:370-374`) is also a behavioural difference. Off in both reference configs. |
| `down_sample_factor` (noise section) | see §5 — it is a fragment-prediction knob, not a renderer knob | — | — | Listed in `[noise_settings]` in the reference v1 config, which is misleading; it is consumed at `simulator.py:2441/2463`. |
| `seed = 41`, `sample_seed = 41` | seeds the **global** `np.random` stream — reproducible but order-dependent (`simulator.py:1649-1660`) | per-tool seeds: `timsim-digest --seed`, `timsim-precursors --seed` (default 42), `[design].seed` (default 42), `--noise-seed` | **(b)** | v1: one global stream, so any change upstream reshuffles every draw downstream. v2: identity-keyed hashes, so nothing reshuffles. The seeds are not interconvertible; matching them is not possible and not meaningful. **v2 has no equivalent of "run the same seed and get the same peptide *set*", because v2's set does not depend on draw order in the first place.** |

## 8. v2-only settings with no v1 counterpart (the reverse (c) rows)

These bound the comparison in the other direction: the v2 arm is being run with knobs the v1 arm
cannot express, so a difference in outcome may be attributable to them.

| v2 setting | meaning | class |
|---|---|---|
| `[design] load_ng` | total peptide mass on column, per run; sets the absolute amount scale by mass balance (`spec.rs:56-59`, `design.rs:906`) | **(c)** |
| `[abundance] source = "table"` | explicit per-protein amounts from a file — "the only source that gets protein IDENTITY right" (`spec.rs:87-89`) | **(c)** |
| `[[condition]].mix`, `.regulate`, `[variance]` | declared A/B design, per-protein fold changes, hierarchical replicate variance (`spec.rs:99-113`, `:259-272`) | **(c)** |
| `--digestion-efficiency` | per-site cleavage probability → peptide yield (`yield.rs:48-53`) | **(c)** |
| `--min-peak-intensity` | drop quantised `(scan, tof)` bins below this count (`render.rs:103-108`) | **(c)** |
| `--noise-only` | background-only control render, for FDP subtraction (`render.rs:182-187`) | **(c)** |
| `--render-chunks`, `--no-parallel`, `--no-memory-guard` | streaming/parallelism controls | **(c)**, no model effect |
| `timsim-yield --report` | missed-cleavage and truncation-loss accounting | **(c)**, diagnostic |

---

# Class counts

70 rows, of which 62 are v1→v2 and 8 are the reverse (§8).

| class | rows |
|---|---|
| **(a) exactly equivalent** | **16** |
| **(b) approximately equivalent** | **25** |
| **(c) not representable** | **29** (21 v1→v2, 8 v2→v1) |

Four of the 16 (a) rows carry a stated precondition and are **conditional**: the binomial charge PMF
(needs a flow edit), the fragment-intensity model (needs the `timsim-fragments` provider pinned),
`reference_noise_intensity_max` (holds only for references with < 1000 mobility scans), and
`n_proteins` (holds only when the FASTA is no larger than `n_proteins`). Read as 12 unconditional
(a) rows.

Four of the 21 v1→v2 (c) rows are **v1 settings that are dead in v1's own code** — the
`[peptide_intensity]` block plus `sample_occurrences`, the `isotope_k`/`isotope_min_intensity`/
`isotope_centroid` triple, `charge_state_one_probability`, and `intensity_multiplier` on this path.
They matter only because the known-good v1 config sets them, which invites the false impression that
they are controlling something.

The (c) rows that actually bite are: **the void-volume correction (`exclude_accumulated_gradient_start`,
default ON, deletes peptides), `min_charge_contrib`, v2's digestion-yield model, `--intensity-scale`,
EMG tailing, and mzPROV provenance.**

---

# Recommended paired configuration

The goal is **estimand B** of the head-to-head plan — a controlled comparison — reached by pushing as
many (b) rows toward (a) as the code allows, and *declaring* the rest. Every choice below is a
deliberate deviation from a default and is justified in one line.

## v1 config (deltas from `/scratch/timsim-demo/v1-diapasef/config.toml`)

```toml
[main_settings]
reference_path = "/media/hd02/data/raw/dia/blanks/blanks-dia-PASEF/G241217_011_Slot2-2_1_16312.d"
acquisition_type   = "DIA"
use_reference_layout = true
apply_fragmentation  = true
emit_provenance = false     # (c) mzPROV: v2 cannot emit it; keep the two .d files structurally comparable

[peptide_digestion]
num_sample_peptides = 10000   # matched to v2 --max-peptides; see the (b) note — this is the distorting knob
missed_cleavages = 2
min_len = 7
max_len = 30
cleave_at = "KR"
restrict  = "P"
n_proteins = 1000000          # >= |FASTA|, so v1's protein SAMPLING is disabled and the (b) row collapses to (a)
upscale_factor = 100000       # left at default; the absolute scale is calibrated on the v2 side instead

[distribution_settings]
gradient_length = 1861.0      # the reference .d's own span, so v1's seconds-anchored RT matches v2's frame-anchored one
target_p = 0.999
sampling_step_size = 0.0001
# sigma_*_rt / k_*_rt left at defaults -> sigma ~ 1.51 s, EMG tail tau ~ 0.72 s.
# NOT matchable to v2's 3.16 s symmetric Gaussian; declared, not tuned.

[charge_state_probabilities]
p_charge = 0.8                # matches v2 --charged-probability default; 0.5 in the old config was arbitrary
min_charge_contrib = 0.005    # v1's own default, NOT the old config's 0.25 -- 0.25 has no v2 counterpart at all
max_charge = 4
binomial_charge_model = true
normalize_charge_states = true

[noise_settings]
mz_noise_precursor = true
precursor_noise_ppm = 6.5     # 3-sigma envelope, the convention v2 shares exactly
mz_noise_fragment = true
fragment_noise_ppm = 6.5
mz_noise_uniform = false      # v1's DIA uniform branch is asymmetric and NOT ported to v2
add_real_data_noise = true
reference_noise_intensity_max = 150000.0   # match v2 explicitly; the two DEFAULTS differ 5000x
precursor_sample_fraction = 0.2
fragment_sample_fraction  = 0.2
num_precursor_noise_frames = 5
num_fragment_noise_frames  = 5
noise_frame_abundance = false
noise_scan_abundance  = false
down_sample_factor = 0.0      # v2 has no stochastic fragment cull; 0.0 removes the 50% peak deletion

exclude_accumulated_gradient_start = false   # v1-only peptide DELETION; off, so both tools render the same set
projection_mode = "off"
quad_isotope_transmission_mode = "none"
precursor_survival_min = 0.0
precursor_survival_max = 0.0
seed = 41
```

## v2 job TOML

```toml
".pipeline" = "/scratch/timsim-demo/timsim-necro-repo/flow/timsim_flow.py:job"
".requests" = ["raw"]

bruker_reference = "/media/hd02/data/raw/dia/blanks/blanks-dia-PASEF/G241217_011_Slot2-2_1_16312.d"
proteome_spec = "<one [[source]], the same FASTA v1 digests>"
mods          = "<Carbamidomethyl C, occupancy = 1.0>"   # v1 applies it as a FIXED mod; 0.98 would not match
design_spec   = "<one condition, one replicate, [abundance] source = hockeystick, load_ng calibrated>"
sample        = "A_R1"
max_peptides  = 10000        # == v1 num_sample_peptides; the matched cap, with the caveat above

collision_energy = 25.0      # DECLARED, not matched -- see the CE box; this is estimand A on the CE axis
intensity_scale  = 5.0e5     # calibrate against a v1 render before freezing; absolute scales are not commensurable

noise_mz_ppm    = 6.5        # same 3-sigma convention as v1
noise_frag_ppm  = 6.5
noise_real_data = true
noise_precursor_frames   = 5
noise_fragment_frames    = 5
noise_intensity_max      = 150000.0
noise_precursor_fraction = 0.2
noise_fragment_fraction  = 0.2
noise_seed = 0
```

### The three flow edits this configuration needs, and why

The job TOML above cannot express three of the matched settings, because the flow does not surface
them. All three are one-line changes to `timsim-necro-repo/flow/timsim_flow.py`; **none has been
made**, and until they are, the pairing above is not achievable.

1. **`--charge-model binomial --charged-probability 0.8`.** The `precursors` rule hard-codes
   `charge_model=cfg.charge_model` with `cfg.charge_model = "site-specific"`
   (`timsim_flow.py:1598`), and `--charged-probability` is not passed at all. Without this the v2
   arm runs site-specific charges (N-term 0.93 / R 0.97 / K 0.95 / H 0.80) against v1's binomial —
   a difference in the *precursor m/z population*, i.e. in what the DIA windows transmit.
2. **`--isotope-depth`.** Not passed (`timsim_flow.py:464-472`), so the envelope is 6 peaks against
   v1's threshold-based envelope. Raise it if the peptide set includes long peptides.
3. **`--sigma-frames` / `--sigma-scans`.** Not passed (`timsim_flow.py:679-685`), so the defaults
   (30 frames, 4 scans) run. The chromatographic width **cannot** be matched to v1's EMG regardless
   — but the *mobility* width can: `--sigma-scans 8.33` reproduces v1's 0.009 `1/K0` mean on this
   reference, and would move that row from (b) to near-(a) at the population level (v1's *per-ion*
   width still has no v2counterpart).

Freeze whichever choice is made, hash the flow file into the manifest, and report the flow's SHA
alongside the results.

---

# Claims bounded by this mapping

## What the benchmark **can** conclude

- **Digestion and the peptide space are directly comparable.** Enzyme, missed cleavages, length
  filter and (with `n_proteins` ≥ |FASTA|) the protein universe are (a). A difference in the number
  of *precursors* rendered is therefore attributable to the ion layer, not to the digest.
- **Precursor m/z is comparable.** Both compute exact-composition monoisotopic masses and both apply
  the identical mscore charge arithmetic when v2 is put in binomial mode. Prior work already
  measured precursor m/z as identical.
- **Mobility placement is comparable.** Same CCS model (parity measured exact on 40,509 precursors),
  same Mason–Schamp constants, same reference calibration. Where an ion *sits* in the mobility
  dimension is a matched quantity; how *wide* it is, is not.
- **A1 m/z noise is comparable in its marginal.** The ppm→σ convention is identical (ppm = 3σ), so a
  search engine's mass-error calibration sees the same distribution from both arms.
- **The A2 real-data background is comparable**: same pools, same per-peak Bernoulli, same intensity
  filter on this reference.
- **Cost is comparable** — wall time, peak RSS, bytes on disk, and the marginal cost of an additional
  sample — provided the workload is genuinely the same, which §1–§3 establish.
- Therefore: **"v1 and v2 produce comparable identification performance on Bruker dia-PASEF under
  their respective recommended configurations, with the differences enumerated in the mapping"** —
  the head-to-head plan's **estimand A**.

## What the benchmark **cannot** conclude

- **It cannot conclude "as accurate" in the sense of estimand B.** Estimand B requires a matched CE
  model. There is none: v2's flow uses a flat 25 eV where v1 uses the reference `.d`'s 20–58 eV
  window table (−33.1 to +5.0 eV per ion), and the opt-in `timsim-frag-ce` implements a *different*
  CE model (dda-PASEF ramp, −14.85 to +9.72 eV against the same table) with a different cardinality
  (one CE per `(sequence, charge)` vs one per transmitted `(ion, window group)`). Estimand B on the
  CE axis is **blocked, not merely unwired**.
- **It cannot attribute an MS2-scoring difference to anything in particular.** Three things move
  together on the MS2 axis: the CE model (above), the peak-culling rule (v1 stochastic 50 % vs v2 a
  0.1 % floor — mitigable by `down_sample_factor = 0.0`), and the fragment-model provenance (two
  packages ship `timsim-fragments`). Report MS2 metrics as *jointly confounded* unless the arms are
  varied one at a time.
- **It cannot conclude anything about chromatographic peak shape or extracted-ion-chromatogram
  quality.** v1: EMG, σ ≈ 1.51 s, right tail τ ≈ 0.72 s. v2: symmetric Gaussian, σ = 3.16 s. **2.09×
  wider and untailed**, with no knob that reconciles them. Any metric sensitive to peak shape — DiaNN
  co-elution scores, apex-based quantification, peak-width features — is measuring this difference,
  not the simulators' relative fidelity.
- **It cannot conclude anything about mobility peak width or per-ion mobility resolution.** v2 is
  **2.08× narrower** in σ and discards `ccs_std` entirely, so it has no per-ion width where v1 has
  one. Matching the population mean via `--sigma-scans 8.33` fixes the first half, not the second.
- **It cannot compare absolute intensities, dynamic range, or anything with an intensity threshold in
  it.** v1's scale is `hockey-stick × upscale_factor × flyability`, cast to int32 (an undeclared
  1-event floor); v2's is `load_ng` mass balance × `--intensity-scale`. There is no conversion.
  Report intensity results as **rank or ratio only**, never as counts.
- **It cannot make any per-protein or per-peptide abundance claim.** v1's abundance is draw-order
  keyed, v2's is identity-keyed. The marginal distributions are the same family; the assignment of a
  value to a *named protein* is unrelated between the two. A protein-level scatter plot across arms
  would be noise.
- **It cannot claim the two rendered the same peptide set** unless
  `exclude_accumulated_gradient_start = false` and `min_charge_contrib = 0.005` are set in v1. At the
  reference config's values, v1 deletes low-RT peptides that v2 keeps, and keeps only charge states
  above 25 % that v2 keeps in full.
- **It cannot claim anything about replicates, differential expression, or fold-change recovery.**
  v1 has no design axis; v2's `[[condition]].regulate`, `[variance]` and `load_ng` are (c). Any
  multi-sample result is a v2 capability statement, not a comparison.
- **It cannot claim anything about DDA, Thermo, SCIEX or Waters.** Out of scope by construction, and
  v2's A1 m/z noise is DIA-only anyway (`run_dda`'s projector applies no ppm scatter,
  `render.rs:1122-1123`).
- **Reproducibility is asymmetric and must be stated.** v1's m/z noise uses unseeded
  `rand::thread_rng()`, so **a v1 noise arm cannot be reproduced bit-for-bit**, and the
  self-consistency check ("each tool run twice and shown to agree with itself") will fail for v1 by
  construction on the noise arms. Pre-register a tolerance, or run the self-consistency check on the
  noiseless arm.

---

# Flagged: counterparts that exist but are untested at benchmark scale

Each of these has a v2 implementation, and each has evidence only well below the scale a benchmark
ladder would reach.

1. **`--charge-model binomial`** (`precursors.rs:112-115`). Claimed v1-exact at `ionize.rs:113-114`
   and reducible to the binomial in principle (`ionize.rs:241-242`), but **the flow has never run
   it** — every rendered artifact on disk used `site-specific`. Unexercised at any scale.
2. **`--collision-energies` / `timsim-frag-ce`** (`frag_ce.rs`, `fragments.py:404-413`). Unit-tested
   (`test_fragment_mobility_ce.py`, `mobility_ce.rs:171-232`) and never wired into a pipeline. The
   claimed 11,164,886-row / 893-scan agreement measurement is **not present anywhere in this tree**;
   treat the capability as untested end-to-end.
3. **`--noise-real-data` at full-proteome scale.** The realism port is recorded as complete and
   flow-wired, but the A2 background sampler holds a per-frame working set; its cost at 9 M
   precursors is not measured. The `fragments` node's RAM declaration carries an explicit warning
   that an earlier "streamed to ~1 GB" estimate was wrong at scale (`timsim_flow.py:586-600`).
4. **`--spike-into`** (`render.rs:862-920`). Validated on small runs; the exact 1:1 frame mapping it
   requires (`render.rs:399-406`) has not been exercised against a 17,646-frame reference.
5. **`--max-peptides` at the benchmark cap.** `design_covid_full.toml` records that a small cap
   makes planted effects unobservable and that "a DE cohort needs the cap off (or ≥ a few 100 k)".
   The interaction between the cap and recall at the 10 k setting proposed above is not characterised.
6. **`--sigma-scans` / `--sigma-frames` as anything but their defaults.** Both are plumbed to the
   binary and neither has ever been set by the flow, so no rendered artifact exists at a non-default
   peak width.

---

# What I could not determine from the code

- **The claimed CE parity measurement.** "0.0 max absolute difference over 11,164,886 rows / all 893
  reachable scans" appears in no file, log or test under `/scratch/timsim-demo`. The nearest real
  result is `mobility_ce.rs:178-188`, a self-consistency assert against the closed-form ramp. If the
  measurement exists elsewhere, it should be attached to the manifest; if it does not, the CE row
  stands at (b) as written.
- **Which `timsim-fragments` the flow actually invokes.** Two packages install that console script
  and only one implements `--collision-energies`. No venv on this box currently has either on
  `PATH`, so I could not resolve it. **This must be pinned in the manifest before any run.**
- **Whether v1's window-group CE genuinely varies per *rendered* ion** on this reference, as opposed
  to varying across the window table. The table spans 20–58 eV and `ion_map_fn_dia` keys CE by
  `(frame, scan)`, so it should — but I did not execute v1 to confirm the realised per-ion CE
  distribution. The −33.1/+5.0 eV figure is derived from the table, not from a v1 run.
- **Whether v2's `--sigma-scans 8.33` really reproduces v1's mobility width.** The conversion
  linearises a Bruker ModelType-2 calibration. The nonlinearity is mild but unquantified here.
- **Whether the `int32` cast in v1's `events` (`simulate_peptides.py:156`) removes a material number
  of peptides** at the benchmark's `upscale_factor`. It is a hard floor at 1 event; how many
  peptides it deletes depends on the flyability draw and was not computed.
- **v2 seed consistency across tools.** `timsim-design`'s seed lives in the TOML (default 42) and
  `timsim-precursors`'s is a flag (default 42). I found nothing enforcing that they agree; a
  mismatched pair would be silently accepted.
- **`num_peptides_total = 250000`** — whether it binds before or after `num_sample_peptides` for a
  given FASTA. It is passed to `simulate_peptides` (`simulator.py:2013`) and inert when it exceeds
  the digest size, which is the case for the small fixtures, but it should be checked per run.
