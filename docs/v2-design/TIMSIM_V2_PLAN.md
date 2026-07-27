# timsim v2 — implementation plan

Companion to `TIMSIM_V2_SPEC.md`. **Revision 4** — S0 is built; this records what shipped, what
it cost, and what building taught us that planning did not.

**Status: S0 complete.** `timsim-chem`, `timsim-schema`, `timsim-cli` (four binaries). Workspace
green: 26 suites, 240 tests. **v1 untouched and still passing.**

---

## 0. Implementation philosophy

**Go slow.** The tool took years to get right; the failure mode we are avoiding is a half-finished
rewrite that leaves two broken things instead of one working one.

- **Small, self-contained, well-tested command-line tools.** Each does one thing and ships with its
  own tests before it is wired into anything.
- **Rust first, for everything that can be Rust.** The early protocol is pure chemistry and touches
  no model.
- **Explicit, permanent exception: the deep predictors stay Python.** RT, CCS, charge, and
  fragment-intensity models are PyTorch. **ONNX is not a goal.** The CLI/schema contract makes the
  implementation language a private detail, so this costs nothing architecturally.
- **Two frontends, one core.** Every Rust crate exposes a `clap` CLI (for necroflow, which shells
  out) *and* PyO3 bindings (so `imspy` remains a Python library). The CLI is not a replacement for
  the Python API — it is a second door onto the same room.
- **No tool merges without an independent oracle** (§6.1). Going slow is what *buys* the ability to
  write independent reference calculators; you cannot do that on a deadline.
- **v1 keeps working throughout.** Every strangler step is a legitimate stopping point.

---

## 1. The case, from the evidence

Not a speculative redesign. The benchmark suite accompanying the TimSim paper
(`/scratch/timsim-demo/SUBMISSION/timsim-bench`, ~6.2k LOC + 21 notebooks) **already implements, in
pandas, most of what v2 provides** — badly, once per benchmark, unreusably. **The refactor is
repatriation, not new capability.**

| What the paper needed | How it was actually done | v2 status |
|---|---|---|
| A/B fold change (HYE) | two simulations; condition recovered from a **filename substring**: `df["sample"] = df.filename.apply(lambda s: "B" if s.find("001.d") == -1 else "A")` (`quant_utils.py:11`). Expected ratios (1.0/0.667/3.0) **typed into the plot cell** — the simulator never knew them. | ✅ `timsim-design`: the mixture is the spec; `true_log2fc` is **derived** |
| Quantification truth | `intensity = events * relative_abundance`, re-derived in ≥3 places | ✅ `amount_amol`, first-class |
| Protein-level FDR | **abandoned, in writing** (`IdentificationBenchmark` cell 18): *"proteins were chosen from a subset … to increase simulation speed … cannot be reported since calculation needs to take degenerate peptides into account"* | 🟡 `peptide_occurrences` gives exact degeneracy (**138,173** repeat/shared occurrences a list column would destroy). Still needs a `protein_groups` artifact + scoring policy |
| Phospho localization | **10 simulations** (phospho on site A, then site B), site recovered from the filename, truth **fabricated**: `SIM["site_probability"] = 1.0` | ✅ `timsim-modify` places the mod (truth by construction); the adapter annotates each modform (`[UNIMOD:id]`) so it **fragments as modified** — verified 80% of fragmented ions modified, 770 backbones rendered as distinct positional isomers; `timsim-localization` emits the site-determining fragments to join against them, separating *unresolvable* from *engine missed*. Remaining: wire the join into `groundtruth_eval` |
| MBR peak matching | apex/FWHM **not emitted**; ~200 lines re-derive them from stringified arrays in SQLite, with a hand-rolled cache | 🔴 needs `run_features.parquet` (realised per-run shifts + cross-run correspondence) |
| Replicates / variance | separate whole simulation runs; the MBR notebook opens three DBs *just to check the replicates differ* | ✅ declared biological/technical variance, kept separable |
| Reproducibility | **zero `.toml` configs in the entire git history**, despite the README telling reviewers to run `timsim config.toml`. **The suite cannot regenerate its own inputs.** The experiment grid lives in *directory names*, re-parsed with regex in ≥6 places | ✅ TOML specs; necroflow `__grid` next |
| Column names | `retention_time_gru_predictor` hardcoded in **raw SQL in the benchmark repo** (`MBRBenchmark.py:845`) | ✅ `timsim-schema` |
| Blank subtraction | the `.d` is built on a **real blank**, so real peptides leak in and are scored as false positives — 5 copies of a subtraction workaround | 🔴 open (§8.6) |

**Absent entirely from v1:** TMT/isobaric, cross-instrument, entrapment FASTA, dilution/LOD series,
any sweep mechanism.

**Scale.** `TIMSIM-HELA-001` = 588 MB for 50k peptides. ~40 datasets → **~300 GB across 6 Zenodo
records**, split *"due to Zenodo's 50 GB limit"*. Compute is not theoretical — **it already cost a
benchmark.**

---

## 2. Architecture

A **cache model, not a cross product** (SPEC §1). Structure → Quantity, then
Structure × Quantity × Measurement → Observation. Known leaks (DDA selection, suppression changing
charge fractions, run-order carryover) are **stated in the spec, not papered over**.

---

## 3. Non-negotiables

- v1 keeps working. No big-bang cutover.
- A **scientific-contract oracle with independent reference calculators** before the physics changes.
- Every step ships something independently useful.
- **The stage CLI + schema contract is the deliverable, not necroflow.** If necroflow doesn't work
  out, the decomposition stands. This de-risks betting on a v0.0.3 single-author framework.

---

## 4. Strangle, don't rewrite

### ✅ S0 — schema + the protocol tools *(DONE)*

```
timsim-proteome → timsim-digest → timsim-design → timsim-yield
    FASTA          STRUCTURE        QUANTITY        QUANTITY
                   (once)         (the mixture)    (per sample)
```

Human proteome, end to end: 20,590 proteins → 2.5M peptides → 2.7M occurrences → 6 samples /
12 runs. Mass balance exact to **1.4e-12 ng**. Structure 96 ms once; quantity 29 ms per sample.

### ⏭ S0.5 — the candidate-universe measurement *(NEXT — blocking)*

Enumerate the **unpruned** modform × charge universe for a real proteome and *count it*. 2.5M
peptides is the base; the multiplier is unknown. If GPU prediction over it is intractable, the
candidate bound must be tightened **chemically** — never by detectability (SPEC §1). Needs
`timsim-modify`. **A number, not an argument, and the last thing that could force a design change.**

### S1 — cut sagepy from the digest
The cleavage rules are already independently confirmed against Sage (§6.2). Removing the *dependency*
also removes a **benchmarking confound**: using sagepy's digestion to generate data that benchmarks
Sage is a shared assumption between simulator and tool under test.

### S2a — export CCS alongside 1/K₀ *(cheap; do early)*
The predictor **already computes CCS** (`ccs/predictors.py:572`) and converts on the way out.
Mason–Schamp **already exists in Rust, already parameterised by gas mass and temperature**
(`imspy-core/chemistry/mobility.py:24`). Nobody varies the defaults. Persist CCS alongside legacy
1/K₀, round-trip compatible. v1 unaffected.

> **Correction from an earlier revision.** This was called "nearly free". The *physics function* is
> free; the *migration* is not — it needs the predictor API to expose CCS **and its uncertainty**,
> a schema migration touching every legacy reader, and a defined uncertainty transform. And
> **"CCS-calibration benchmarking falls out free" was simply wrong**: benchmarking a calibration
> requires *modelling* one. That is S2b.

### S3 — the strangler adapter
Write the new Parquet artifacts into v1's legacy SQLite tables, so v1 runs unchanged on v2 inputs.
**This is the step that makes everything reversible.**

### S4 — Parquet in Rust, render layer, vendor registry
The biggest single chunk (§8.1).

### S2b — measurement-owned mobility calibration + scan projection
Sequenced **here, not with S2a** — it needs the measurement abstraction to exist. Only after this
does CCS-calibration benchmarking become real, and it needs a **held-out reference panel**.

### S5 — necroflow DAG, grids, the payoff

---

> **Measurement-axis ordering (corrected).** Peptides separate in solution and only then reach the
> ESI needle, so **LC precedes ionisation** — and ionisation is *coupled* to it (eluent composition
> at elution; who co-elutes). The chain is `load → LC → ionise → ion optics → fragment → acquire →
> render`. Consequence: the structure axis carries m/z, isotopes, and *latent propensities* only;
> the realised charge distribution and response are measurements. The tool to build is
> `timsim-precursors` (structural), not a structure-axis `timsim-ionize`.

---

## 5. What we explicitly do NOT do

- **No ONNX.** Permanent exception for the deep predictors.
- **No cluster/cloud executor.** Necroflow is local-only. SLURM scripts existed in the bench repo and
  were gitignored — if grids need a cluster, that's a separate conversation.
- **No TMT in v2.0.** The mapping model must not *preclude* it; that is all we claim.
- **No GUI rewrite** until the config schema settles (it holds a second, duplicated config model).
- **No new science features** during migration.

---

## 6. The oracle

We change the physics **on purpose**, so "v2 == v1" is **not** the acceptance criterion.

**The correction that inverted revision 1:** treating v1's validation suite as a hard constraint
would make the *legacy schema* the de facto v2 schema — `validate/comparison.py` defines truth as
`events * relative_abundance`, `retention_time_gru_predictor`, `inv_mobility_gru_predictor`. The
suite must be **ported to the new ground truth, not preserved as a constraint on it.**

Aggregate search-level comparison (ID counts, FDR, CVs) **cannot establish equivalence**: identical
aggregates can conceal wrong precursor-level physics, and legitimate changes move all three anyway.

### 6.1 The anti-self-consistency rule

**Expected values must come from an implementation other than the one under test** — an independent
calculator, a held-out measured panel, or a hand-auditable fixture. Never the same stage code that
emits the ground truth. *This is what "go slow" buys.*

### 6.2 Oracles built so far

| Claim | Verified by | Result |
|---|---|---|
| We cut in the right places | **Sage**, `timsim-chem/xcheck/` | identical peptide sets, 5 bounds, real proteome, **zero divergence** |
| The yield maths is right | **Monte Carlo**, 400k simulated digests, non-uniform blocking | agrees within 3σ binomial at every occurrence |
| Mass balance holds | proven analytically; **reproduced by counting** in the simulation | exact |
| The RNG is Gaussian | mean, σ, skew, **kurtosis**, 1/2/3σ coverage | excess kurtosis +0.005 (was +0.091) |

### 6.3 Still required, before the physics changes

- Reference panels: known peptide RT order; published CCS values **per charge**; empirical
  charge-state distributions.
- Chromatography: peak-width and co-elution distributions against the real technical-replicate
  measurements already made in `notebooks/data_modeling/intensity_30min_hela.ipynb`.
- **Counterfactual panels** — small, locked, precursor-level:
  - **DDA Top-N threshold:** two co-eluting precursors straddling the cutoff; changing *one amount*
    must change selection, co-isolation, and exclusion history exactly as specified.
  - **Suppression is local, not global:** two analytes far apart in RT — raising one must *not*
    suppress the other.
  - **Design-conditional detectability:** a peptide observable only in condition B stays in the
    shared candidate universe and is **not lost to A-driven pruning**.
  - **Positional phospho isomers:** equal mass, controlled ratios, both coexisting, with
    **site-determining fragments** differing.
  - **CCS panels** across charge and gas, checked against an **independent** implementation.
### 6.4 The list of INTENDED physics changes

We change the physics on purpose. Every entry below is a *deliberate* departure from v1 — **any diff
not on this list is a bug.** This list is the acceptance criterion; `v2 == v1` is not.

| # | change | v1 | v2 | why |
|---|---|---|---|---|
| 1 | **Charge-state model** | binomial — every basic site at p=0.8, so **histidine is treated exactly like arginine** | site-specific Poisson-binomial (N-term 0.93, R 0.97, K 0.95, **H 0.80**) | a binomial cannot express an ordering of basicity, so it is necessarily wrong about both ends of it: it calls **33% of clean tryptic peptides singly-charged**, which does not happen. v2 gives them 11.3% / 88.7% at 2+. The two models are tuned to **agree on the aggregate** of an observable tryptic population (~79% 2+, ~18% 3+) and to **disagree on which peptides** carry the charge — which is the thing a search-engine benchmark actually tests. `--charge-model binomial` reproduces v1 exactly. |
| 2 | **Modifications** | *variable mods*: "oxidation on M, phospho on STY, **max 3 per peptide**" → generate every combination | **occupancy**: the fraction of *this site*, across all molecules, that carries the mod. The modform distribution follows from it | "max N variable mods" is a **database search parameter** — the space a search engine looks in, not the population in a tube. Simulating from it makes the simulator and the search engine agree by construction, and a benchmark built on that agreement cannot detect the failure it most needs to. Occupancy is the chemist's number and it is measurable (phospho ~0.02, Met-ox ~0.05, carbamidomethyl ~0.98 = the alkylation efficiency). Measured on HYE: **18.2M modforms over 3.5M peptides**, mean truncation 0.61% at floor 1e-3, and `Σ abundance_fraction = 1 − loss` holds to 5.9e-5. |
| 2b | **Cleavage-blocking mods** | none — a modified lysine is cleaved exactly like a bare one | `p_eff(k) = p · (1 − blocking_occupancy(k))`, read from `modifications.parquet` | acetyl-K, GG-K, trimethyl-K and TMT-K **physically stop trypsin**, and the missed cleavage that forces at the modified residue *is the evidence that localises the site*. A diGly simulation whose protease ignores the GG is not a simulation of a diGly experiment. Measured (GG on K at 5% occupancy, 926,382 sites blocked): missed cleavages **0→90.3% / 1→8.9% / 2→0.9%** becomes **0→88.0% / 1→10.7% / 2→1.3%**. The missed cleavage is **not a parameter** — it emerges. |
| 3 | **CCS, not 1/K₀** | stored `1/K₀` (`inv_mobility_gru_predictor`) — implicitly N₂ at 305 K, the gas hardcoded and never varied | store **CCS** on the ion (structure, `timsim-ccs`, deep model); derive `1/K₀` per run via Mason–Schamp with explicit `drift_gas_mass`/`drift_temp` (measurement) | `1/K₀` is not a property of the ion — it is what a *particular* drift tube measures, and it moves with the gas though no molecule changed. CCS is the ion. Splitting them is B14, and it is what makes **cross-instrument** an experiment rather than a rebuild: same `precursor_ccs`, run with a different gas, is instrument B (measured Ar/N₂ mean factor 1.19). CCS is recovered losslessly by inverting the model's own default-gas `1/K₀` (round-trip 2e-16), so **at the default gas the numbers are bit-identical to v1's** — the gas is the only knob that moves them. |
| 3b | **Mobility peak width** | deep model's per-ion CCS std, rescaled so the population mean hits a target (`use_target_mean_std`) — the deep std is informative in *shape* but too wide in *scale* | same calibration, kept: **`ccs_std_model = "predicted"`** (deep per-ion shape, mean rescaled to `inverse_mobility_std_mean`), with **`"proportional"`** (a fixed CV of 1/K₀) as the model-free alternative | getting the CCS spread right was hard-won in v1; the deep std carries real per-ion information (some peptides predicted more confidently) but its absolute scale is miscalibrated. Rescale (not shift) preserves the ratios and fixes the scale. Verified: rescaled mean 0.0091 → 0.0090, per-ion range 0.0041–0.0146 preserved. **Caveat:** the rescale is a population operation (each ion's std depends on the batch mean) — matches v1 exactly, mild B8 non-locality. |
| 4 | **RT: index, not seconds** | GRU RT, linearly stretched so each SAMPLE's min..max fills [0, gradient] — a peptide's RT moved with whatever else was in the tube | **Chronologer** index per peptide (structure, `timsim-rt`); map index → seconds per run with a **fixed reference range** carried in the artifact (measurement) | seconds are not a property of the peptide; hydrophobicity is. Storing the index and mapping per run (B14) makes runs comparable across gradients — the RT analog of cross-instrument mobility. And it is *more physical* than v1: a fixed reference range (index span over the whole peptide space) means a peptide lands at the same gradient **fraction** regardless of sample. Verified end to end: index [-0.02, 22.73] → [0, 600]s, 4.4 MB .d; portability + sample-independence pinned by tests. Chronologer weights live locally; the upstream `chronologer` pkg is a dependency. |

**Everything else so far is a refactor with identical numbers**, verified: the digest reproduces
Sage's peptide sets exactly, the yield maths is Monte-Carlo confirmed, and the fold changes are
exact. Change #1 is the **first** deliberate divergence — so a v2-vs-v1 diff in charge states is
expected, and a diff in anything else is not.

- **An explicit, reviewed list of intended physics changes** — any diff not on the list is a bug.

---

## 7. What building taught us that planning did not

### 7.1 B13 — never re-enter or infer a fact the artifact already knows

The single most productive invariant, and it is **not in any earlier revision of the spec**. Four
bugs were the same bug:

- `timsim-yield` **inferred** protein length as `max(end)` over occurrences → the length filter
  discards C-terminal peptides → inferred length fell below some cleavage sites → boundary lattice
  non-monotonic → integer underflow.
- `timsim-yield` **re-took** `--max-missed-cleavages` → digest at 4, yield at a defaulted 2 → every
  3- and 4-missed-cleavage occurrence silently got **zero yield**. A wrong number, no error.
- A required column arriving **nullable** was accepted → read with `.value(i)` → **garbage rather
  than a failure**.
- The v1 adapter **redrew flyability** in Python while `precursors.parquet` already carried
  `ionization_propensity` for every peptide. Same distribution, different key (blake2b on the
  sequence vs. the Rust side's FNV/SplitMix identity key) — so the two were **independent draws of
  one physical quantity**. Measured: per-peptide correlation **r = +0.008**, and **43.6% of peptides
  differed by more than 10×** between the propensity that built their intensity and the propensity
  recorded as their ground truth. The declared chain
  `ion_amol = peptide_amol × ionization_propensity × charge_fraction` silently stopped inverting.
  Caught by Codex, not by me — and note **what** hid it: both draws share a marginal distribution,
  so every histogram, median and summary statistic looked perfect. **A quantity can be right in
  aggregate and meaningless per entity**, and only a per-entity check sees it. Fixed by reading the
  artifact; the chain now inverts to 0.0 relative error out of the `.d`.

Fix in every case: the fact travels **with the artifact**, and the flag that let two stages disagree
**no longer exists**. The schema crate is not just about column *names* — anything a consumer must
retype by hand is a place two stages can silently diverge. The fourth instance adds a corollary:
**this applies to random draws too.** An identity-keyed RNG makes a value reproducible *within* one
implementation; it does nothing to make two implementations agree. If a drawn value is ground truth,
exactly one component may draw it.

### 7.2 Four instruments, four failure classes — none substitutes for another

| Instrument | Caught |
|---|---|
| Unit tests (45) | boundary conditions, invariants, model semantics |
| Independent oracles (3) | that the *maths* is right |
| **`codex review --uncommitted`** (6 bugs) | **integer overflow, input validation, cross-stage desync** — none of which the oracles could see, because they live *outside* the model |
| **End-to-end run on real data** (1 bug) | the length-inference underflow. **No unit test would have found it** — it only appears when the length filter discards a C-terminal peptide |

`codex review --uncommitted` is fast (~2 min), diff-native, and clearly earns its keep. Run it after
every tool.

### 7.3 Tests that cannot fail against their own bug are theatre

Written twice, caught twice, both times only because the principle had just been stated out loud:

- The first `u8`-wrap regression test used `max_missed = 2` — but the loop breaks at 3 and never
  reaches 256, so it would have **passed against the buggy code**.
- The `gauss` kurtosis threshold was 0.10; the bug measured 0.091 — it would have **squeaked
  through**. Re-measured, recalibrated to 0.03.

**Calibrate every regression test against the implementation it was written to catch.** If you can't
show it failing, it isn't a test.

### 7.4 Claims must be measured, not asserted

The RNG defect was real (excess kurtosis +0.091, σ 1.0076) but I first billed it as "silently
corrupts every abundance draw". It doesn't — it mildly fattens the tails. Measured, corrected, moved
on. **A wobble reported as a catastrophe costs credibility that the next real catastrophe needs.**

### 7.5 A parameter asserted from a true fact is still asserted — `p(H) = 0.20`

The site-specific charge model shipped with `histidine: 0.20`, reasoned from a fact that is *true*
("histidine is a weaker gas-phase base than arginine") but that does not contain a number. The
number was invented. The relevant chemistry — His pKa ≈ 6.0 against an ESI solvent at pH 2–3, so it
enters the droplet **fully protonated** — puts it *below* K/R, not at a fifth of them. Correct value
~0.80.

It mattered far more than a parameter tweak should, and the reason is structural: **a fully-tryptic
peptide ends in K/R and has no internal K/R, so its only third basic site is almost always a
histidine.** The 3+ fraction of a tryptic digest is therefore very nearly a direct readout of this
one number. At 0.20 the simulator produced **5.7% 3+** against a real 20–30%; at 0.80, **18.3%**.
Every individual peptide's charge distribution looked perfectly plausible the whole time.

Three lessons, in increasing order of how much they cost:

1. **The unit test encoded the belief, not the physics.** `histidine_contributes_far_less_charge_
   than_lysine` asserted `Δ_K > 3·Δ_H` — it passed for exactly as long as the bug existed, and it
   failed the moment the bug was fixed. A test written from the same wrong intuition as the code
   defends the bug. The replacement asserts the *ordering* (R > K > H, which is the real physics)
   plus an **aggregate** the bug destroys but no single peptide reveals.
2. **The justification in the doc comment was a rigged comparison.** It claimed the binomial's mean
   charge of 2.81 was wrong against "~2.2–2.4 for real tryptic peptides" — comparing the *unweighted
   enumerated peptide space* (dominated by long missed-cleavage K/R-rich peptides nobody observes)
   against an *abundance-weighted, m/z-filtered, observed* literature figure. Different populations.
   That comparison would have confirmed whichever model I already preferred, which is precisely why
   it felt like evidence. **Validate on the observable set, never on the enumeration.**
3. **I nearly deleted the right model to protect a wrong constant.** The first diff looked bad
   (v2 94% 2+ vs v1's more-realistic 78.8%) and the tempting read was "the site-specific model is
   broken, go back to the binomial". It wasn't the model; it was one number *inside* the model — and
   the only reason the error was findable at all is that the site-specific model **has a slot for
   histidine**. The binomial cannot be wrong about histidine because it cannot represent it. **A
   model that can be wrong in a specific way is worth more than one that can only be wrong
   diffusely** — do not trade a diagnosable model for an undiagnosable one because the diagnosable
   one just told you something you didn't want to hear.

---

## 8. Sharp edges (remaining)

### ~~8.-1 v1 cannot reproduce its own ground truth~~ — FIXED

**Measured, not inferred.** The same config, the same seed, run twice, on the v2 path:

```
  peptides run A / run B      : 1,675 / 1,703
  in both                     : 1,425          Jaccard 73.0%
  of the SHARED peptides:
    identical events (abundance) : 100.0%
    identical RT                 : 100.0%
```

`np.random.seed` is never called anywhere on the simulation path — the only call in the package is in
a GUI plot. Every `np.random.*` draw comes from numpy's unseeded global RNG, seeded from OS entropy
at import. So **a timsim run is not reproducible**, and for a tool whose entire product is a ground
truth against which search engines are scored, that means you cannot re-derive the answer key of a
published benchmark. You can only keep the `.d`.

Note precisely what does and does not move, because it is the whole argument for **B8**: every
quantity v2 owns — abundance, flyability, charge, m/z — returns **bit-identical**, because it is
drawn from `hash(seed, entity_id)` rather than from draw order. The entire 27% churn is v1's, and it
enters through unseeded bulk draws (`simulate_peptides.py:16,32,51,200`, `digest_fasta.py:106`,
`simulate_frame_distributions_emg.py:298`) that decide *which peptides survive*.

This was nearly missed. A 0.6-point drift in the charge distribution between two runs looked like it
might be a regression from the modform change, and the only reason it was correctly attributed is
that the same config was run twice. **A nondeterministic simulator makes every A/B comparison
unfalsifiable** — any diff can always be the RNG.

**Fixed.** A `seed` config key (default 41) now seeds numpy, `random` and torch before any job runs.
Verified end to end:

```
  seed 41 vs seed 41 (rerun) : Jaccard 100.0%   ground truth byte-identical: YES
  seed 41 vs seed 99         : Jaccard  74.9%   ground truth byte-identical: NO
```

and on v1's own FASTA path, `simulate_proteins`' abundance draws are identical under one seed and
differ under another.

**This is a real behaviour change, and it is the point.** v1 produced *replicates* by re-running the
simulator, and the only thing that made two replicates differ was the unseeded RNG. Replicates must
now vary `seed` explicitly. Variation you did not ask for is not technical variance — it is an
uncontrolled variable that happens to look like it.

**What this does NOT buy.** A global seed makes a run reproducible but leaves it **order-dependent**:
add one peptide and every draw after it reshuffles. The stronger property — which the v2 structure
axis already has — is identity-keyed randomness, `hash(seed, entity_id)`, so a peptide's value depends
on the peptide and not on its neighbours. That remains the eventual fix for these jobs. This is the
one that stops the bleeding.

### 8.0 The venv runs a *copy* of imspy-simulation, not the source tree

`venvA` has `imspy_simulation` installed into `site-packages` as a **copy**, not as an editable
install. The v2 adapter therefore only runs with

```bash
PYTHONPATH=<repo>/packages/imspy-simulation/src
```

and without it the simulator **silently falls back to the FASTA path**: `SimulationConfig.__getattr__`
returns any key present in the TOML, so `config.from_v2` reads back `True` from a config the running
code has never heard of, the `if config.from_v2:` branch does not exist in the installed copy, and
the `elif` runs a perfectly normal v1 digest instead. Nothing errors. The `arguments-*.txt` dump even
records `from_v2 = True`, because it dumps the config dict rather than the code path taken.

This cost a real debugging detour — the run *looked* configured and produced a plausible `.d`. Two
defences, neither yet built: **(a)** `SimulationConfig` should reject unknown keys rather than
absorbing them, which turns a silent misconfiguration into a startup error; **(b)** the v2 branch
should log the resolved `simulator.py` path, so "which code am I running" is in the artifact.

### 8.1 Rust reads SQLite by column name — everywhere
`rustdf/src/sim/handle.rs` (1944 L), `dia.rs` (1430), `dda.rs` (1086), `lazy_builder.rs` (774) each
open their own `rusqlite` connection and read columns by string. **Largest single work item; most
likely to be underestimated.**

### 8.2 Constructors with filesystem side effects
`acquisition.py:_setup()` writes SQLite tables as a *constructor side effect*; `TDFWriter::__init__`
creates the `.d` before any simulation happens. **You cannot dry-run today.**

### 8.3 The `acquisition_builder` god object
Jobs reach *through* it (`acquisition_builder.tdf_writer.helper_handle.im_lower`). Cannot pass
through a typed file-product edge.

### 8.4 Mixture coupling / ion suppression — the biggest realism gap
v1 models it **not at all**. Genuinely hard; may be v2.1.

### 8.5 Two config models, 131 flat keys
`get_default_settings()` (131 keys, no schema, `__getattr__` — a typo raises `AttributeError` an hour
into a run) **and** `gui/state.py` (8 dataclasses, independently duplicated defaults).

### 8.6 The reference-blank dependency
Real peptides from the blank leak into "synthetic" data and are scored as false positives — hence
five copies of a blank-subtraction workaround. **A scientific-validity issue, not ergonomics.**

### 8.7 Necroflow's own gaps
- **`Node.fingerprint` is an uncached property recursing over all ancestors** (`nodes.py:77`);
  `Node.key` calls it, and `key` is a dict key in every hot loop. Exponential re-hashing on diamond
  DAGs — **exactly our shape**. `_accumulated_config` (`dag.py:69`) has the same defect. **Blocker
  for large grids. Cheap to fix. Raise upstream now.**
- `_accumulated_config` flattens ancestor configs into one namespace and "assumes config key names
  are unique across the pipeline" — with 131 keys including `path` and `name`, collisions silently
  corrupt the provenance record.
- **No retries, no timeouts.** A GPU OOM is a hard failure. `Constraints(ram=...)` is admission
  bookkeeping; nothing cgroups the process.
- Local `ThreadPoolExecutor` only; one instance per outdir.

### 8.8 Bugs found in v1's *submitted* code — fix independently of this plan
- `table_concat.py:287-291` — appends `maxquant_df` where `spectronaut_df` is intended, and vice versa.
- `ParameterOptimizationExample` searches `...-001.d` but loads ground truth from `...-002/`.
- `MBRBenchmark.py:851` — `'monoisotopic-mass'` in single quotes is a SQL string literal, not a column.

---

## 9. Decisions still needed

1. **Candidate-universe size** — blocking; S0.5.
2. `amount_amol` — "in sample" or "on column"?
3. Contaminants — modelled, or injected via a contaminant FASTA?
4. How far do we take the LC scope declaration (SPEC B1)?
5. Does source response ship in v2.0, or v2.1?
6. Does the reference-blank approach survive v2? (§8.6)
7. **Are we allowed to change the numbers?** If published benchmarks must stay byte-reproducible,
   this plan needs a compatibility mode.
