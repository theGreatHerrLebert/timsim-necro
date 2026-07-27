# timsim v2 — design spec

**Revision 4.** Incorporates two rounds of Codex review on the plan, two on the code, and —
most importantly — what we learned by actually building the first four tools. Where the code
and an earlier revision of this document disagree, **the code won**, and this revision says
so explicitly.

Implemented so far: `timsim-chem`, `timsim-schema`, `timsim-cli` (`timsim-proteome`,
`timsim-digest`, `timsim-design`, `timsim-yield`). v1 untouched.

---

## 0. What we are designing against

> "You push one button down and three stages downstream the program explodes because a
> hard-coded column name is no longer GRU."

Root cause: **the schema is an implicit contract, enforced nowhere** until something reads a
column that isn't there. `retention_time_gru_predictor` — named after a model replaced years
ago — is hard-coded in the v1 Rust SQL reader (`rustdf/src/sim/handle.rs:210`) in 66 places
across two languages, Python searches for it *by substring* (`simulator.py:1540`), and it has
now spread into raw SQL in the published benchmark repo (`MBRBenchmark.py:845`).

Three layers now stand against it, in order of how early they catch you:

| Layer | Catches | When |
|---|---|---|
| Necroflow `NodeType` | wrong artifact into wrong stage | graph construction, before execution |
| `timsim-schema` validate-on-read | wrong/missing/renamed/**retyped/nullable** column | stage entry, before compute |
| Parquet typing | wrong dtype | file read |

### The generalisation we learned by building (B13)

Building the tools taught a lesson the original spec did not contain, and it turned out to be
the **most productive invariant in the whole document**:

> **B13 — Never re-enter, or infer, a fact the artifact already knows.**

Three separate bugs were the same bug:

- `timsim-yield` *inferred* protein length as `max(end)` over occurrences. But the length
  filter discards C-terminal peptides, so the inferred length fell **below** some cleavage
  sites, the boundary lattice stopped being monotonic, and `end - start` underflowed.
- `timsim-yield` *re-took* `--max-missed-cleavages`. Digesting at 4 and yielding at a
  defaulted 2 gave every 3- and 4-missed-cleavage occurrence a **zero yield** — a wrong
  number, silently, with no error at all.
- A required column arriving *nullable* was accepted, and read with `.value(i)` — which on a
  null returns **garbage rather than failing**.

The fix in every case: the fact travels **with the artifact** (protein length as a column;
enumeration bounds as Arrow metadata; nullability as part of the conformance check), and the
flag that let two stages disagree **no longer exists**. Anything a consumer must retype by
hand is a place two stages can silently diverge.

---

## 1. Architecture: a cache model, not a cross product

The true shape is a **partial order with one fan-out** — not three orthogonal axes:

```
   STRUCTURE  ──────────────►  QUANTITY
   (molecules, properties)     (per-sample amounts)
   GPU-expensive, once         cheap linear algebra
        │                           │
        └──────────┬────────────────┘
                   ▼
        STRUCTURE × QUANTITY × MEASUREMENT  ──►  OBSERVATION
                                                 (one run / one file)
```

Quantity **depends on** structure. The value is that **structure is reusable candidate/property
data**, computed once and shared. It buys **caching**, not orthogonality.

### Where it genuinely leaks — stated, not papered over

- **DDA couples quantity into measurement, hard.** Abundance changes Top-N ranking, threshold
  crossings, co-isolation, and the entire dynamic-exclusion history.
- **Source suppression changes observed charge-state *fractions*, not merely total flux** — so
  the charge partition is mixture- and retention-dependent, not a fixed per-method matrix.
- **Observed fragment intensities are not quantity-independent** (isolation transmission,
  co-isolates, AGC, detector saturation). Only the *intrinsic pattern* caches (§6.7).
- **Run order couples runs to each other** via carryover, contamination, and batch drift.

### Why the structure/quantity split is *required*

If `amount_amol` lives in the structure, every sample has a different fingerprint and
**nothing is shared** — a 20-sample A/B design pays 20× GPU prediction. The benchmark suite
proves this is not hypothetical: every multi-run question was answered by running the
simulator N times into N directories.

### Detectability is **not** part of the structure

- **Candidate universe** (structure) — unpruned, bounded by *chemistry*, never by detectability.
- **`design_eligibility`** — a separately versioned artifact derived from
  structure × quantity × instrument.

> ⚠ **Open, and must be answered with a number:** without early filtering, we GPU-predict the
> full candidate universe. `timsim-digest` gives us the first half — **2.5M peptides** for the
> human proteome at 2 missed cleavages. The modform × charge multiplier is the missing factor,
> and it is the **S0.5** measurement. If it explodes, the candidate bound must be tightened
> *chemically* — never by detectability.

---

## 2. Design invariants

- **B1 — Method-invariant, with declared scope.** The structure axis is **not "device-free"**.
  It is invariant under a *declared equivalence class of methods*: a stated LC chemistry family
  and a source-neutral reference state. Recorded as `timsim.scope` in every artifact.
- **B2 — Parametric, not discretised.** Vectors appear only when an instrument discretises onto
  its own grid.
- **B3 — Structure carries no amounts.** *(Enforced by a test:
  `the_structure_axis_carries_no_quantities`.)*
- **B4 — Structure carries no detectability.**
- **B5 — `amount_amol` is a sample-domain quantity and stops at the source.** (§5)
- **B6 — Generative, not enumerative.** We model what a sample *contains*, not what a search
  engine would *consider*.
- **B7 — Content-derived IDs, with collision detection.** `peptide_id = hash(sequence)`, never a
  row index. A 64-bit hash is **not** a uniqueness guarantee: `timsim-digest` runs a
  `CollisionCheck` and fails loudly rather than merging two peptides into one row. Hash recorded
  as `fnv1a-splitmix64/v1`.
- **B8 — Identity-keyed randomness.** Seeds derive from `(global_seed, entity_id)`, never from
  draw order. (§9)
- **B9 — Sample ≠ Run.** Many-to-many.
- **B10 — Two vocabularies, one physics.**
- **B11 — Mixture coupling is real.** Absolute response is not computable on the structure axis.
- **B12 — One name, one random variable.** Never merge distinct variances into one `_sigma`.
- **B13 — Never re-enter or infer a fact the artifact already knows.** (§0)
- **B14 — Place a stage where the *physics* lives, not where today's *model* needs it.**
  A cheap stage in the right slot costs nothing; a stage in the wrong slot costs a
  migration when the model improves. (§8.2)

---

## 3. `timsim-schema` — the contract *(implemented)*

A Rust crate owning the Arrow schema for every artifact. **The only place a column name is
spelled.** Stage code writes `peptide_quantities::AMOUNT_AMOL`, never `"amount_amol"`.

Every artifact carries Arrow key-value metadata:

```
timsim.schema_version = "2.0"
timsim.table          = "peptide_occurrences"
timsim.axis           = "structure"          # structure | quantity | design | measurement
timsim.producer       = "timsim-digest/0.1.0"
timsim.scope          = "lc:c18_water_acn_fa0.1; source:reference"     # B1
timsim.bounds.max_missed_cleavages = "2"     # B13 — travels WITH the artifact
timsim.bounds.min_length = "7"
timsim.bounds.max_length = "50"
```

Validation on read checks table, version, required columns, **dtypes, and nullability**. Extra
columns pass through (a stage may annotate) but are surfaced when something else is wrong,
because an unexpected column is usually a renamed one:

```
peptide_quantities.parquet: table "peptide_quantities" does not conform to schema 2.0.
  - missing column "amount_amol" (Float64)
  ? unexpected columns present: amount_amol_gru_predictor — is one of these a renamed column?
  This is a schema mismatch at a stage boundary, not a bug in this stage.
```

---

## 4. The design axis — samples, runs, mixtures *(implemented: `timsim-design`)*

**Sample ≠ Run.** LFQ is 1→1; fractionation is 1→N; **TMT is N samples → 1 run**. Conflate them
and TMT is unrepresentable and fractionation is a hack — which is exactly how A/B support ends
up bolted on.

> **Honest scope note.** A `channel` column does **not** make TMT supported. Isobaric tags alter
> the chemical entity and need labelling efficiency, isotopic impurity, and reporter cross-talk.
> The claim is only that the mapping model does not *preclude* it. Out of scope for v2.0.

### The mixture is the spec; the fold change is **derived**

At the bench you specify a **mixture**, not a fold change. v1's benchmark suite has the expected
HYE ratios (1.0 / 0.667 / 3.0) **typed into a plotting cell**, because the simulator never knew
them.

```toml
[design]
reference = "A"
load_ng   = 200

[abundance]
HUMAN = { source = "table",     path = "human_paxdb.tsv" }
YEAST = { source = "lognormal", sigma = 2.0 }              # σ≈2 → ~7 orders of dynamic range
ECOLI = { source = "lognormal", sigma = 2.0 }

[[condition]]
name = "A"
mix  = { HUMAN = 0.65, YEAST = 0.30, ECOLI = 0.05 }        # MASS fractions
replicates = 3
technical_replicates = 2

[[condition]]
name = "B"
mix  = { HUMAN = 0.65, YEAST = 0.15, ECOLI = "rest" }      # → 0.20
replicates = 3
technical_replicates = 2

[variance]
biological = 0.15
technical  = 0.05
```

`"rest"` fills to 1.0, so a mixture cannot be mistyped into not summing. Yields
`true_log2fc = 0 / −1 / +2` — computed, not typed.

### Mass → molar, and what cancels

```
amount_i = total_amol(org) · f_i ,   total_amol(org) = mass_ng(org) · 1e9 / ⟨M⟩_org
```

where `⟨M⟩_org = Σ f_i · MW_i` is the molar-weighted mean **average** MW (not monoisotopic —
`--load-ng` is a *bulk* mass).

- **⟨M⟩ cancels from the fold change**, because `f_i` is structural and identical across
  conditions: `ratio = mass_org(B)/mass_org(A)`, exactly.
- **⟨M⟩ does *not* cancel from absolute amounts.** E. coli's smaller proteins mean a given mass
  contributes *more molecules*, which shifts what sits above the detection limit. Mean MW
  therefore affects **sensitivity**, not **fold change** — and only a model carrying the
  sequences can get that right.

### Compositional bias is recorded, not pretended away

Regulation is applied, then the load is restored to `--load-ng` — because that is what happens at
the bench. So a large spike-in **compositionally dilutes the background**, and because
`true_log2fc` is computed from **final** amounts, the answer key records it. A simulator that
wrote the *nominal* fold change into the truth column would hand search engines an answer key
that is quietly wrong.

### Mass balance

```
Σ_proteins  amount_amol × MW_average × 1e-9  ==  load_ng      exactly, per sample
```

Measured error on the human proteome: **1.4e-12 ng**.

### Replicates, and why Sample ≠ Run

- **Biological** replicate — different material. Amounts genuinely differ; variance applied here.
- **Technical** replicate — *the same tube injected twice*. Amounts are **identical**; all the
  variation is on the measurement axis. Two `Run`s, one `Sample`.

---

## 5. `amount_amol` — and where it stops

Absolute attomoles. Composable (digestion splits it, enrichment reweights it, loading normalises
it), correctly scaled (200 ng across ~7 orders → low end at single-digit amol, where real
detection limits sit), legible to both audiences (`copies = amount_amol × 6.022e5`).

Replaces `num_sample_peptides` — an informatician knob with no physical meaning — with
`--load-ng`, a chemist knob with one.

### The quantity ladder (B5)

```
sample domain    amount_amol           digest → enrich → loss → load
    ↓ SOURCE — competitive, mixture-coupled (§7.3)
ion domain       ion flux / current
    ↓ acquisition
counts domain    ion counts
    ↓ detector
signal domain    ADC values            gain, saturation, dead time, baseline, noise
```

Each arrow is a **transfer function, not a multiplication**. Nothing past the source is called
`amol`.

---

## 6. Structure-axis tables *(implemented)*

### 6.1 `proteome.parquet`
`protein_id`, `sequence`, **`length`**, `description?`, `organism?`, `is_contaminant`.

> **`length` is carried even though it is derivable from `sequence`.** Downstream stages need it
> to rebuild the cleavage-boundary lattice, and inferring it from the occurrence table is
> **wrong** — the length filter discards C-terminal peptides, so `max(end)` under-reports. This
> is B13, and it was a real underflow bug.

No amounts. Organism is a **declared column**, not a substring match on the FASTA header — v1
recovers it with `protein_id.str.contains("HUMAN")` and silently *drops* peptides shared between
organisms.

### 6.2 `peptides.parquet`
`peptide_id`, `sequence`, `length`, `mass_monoisotopic`.

> `mass_monoisotopic`, never just `mass`. Monoisotopic is for m/z; **average** is for bulk
> quantities. Conflating them is a real error, so neither gets the bare name.

### 6.3 `peptide_occurrences.parquet` — the digest operator **and** the inference answer key
`peptide_id`, `protein_id`, `start`, `end`, `n_missed_cleavages`.

**There is no `p_yield` column.** Earlier revisions of this spec put one here. That was wrong:
yield depends on digestion efficiency *and* on cleavage-blocking modification occupancy — and
**acetylation is regulated**, so cleavage probabilities differ between conditions. A yield here
would make the structure condition-dependent and silently destroy sharing.

Three jobs, which is why a `protein_ids` *list column* is insufficient:
1. It **is** the sparse digest operator.
2. It preserves **per-occurrence contributions** — how much each protein contributed to a shared
   peptide, not merely which ones did.
3. The same sequence can occur **more than once within one protein**. Measured on the human
   proteome: **138,173 more occurrences than unique (protein, peptide) pairs.** Sage discards
   these (`enzyme.rs:316`), because a search engine does not care. We do.

> **Necessary but not sufficient for protein-level FDR.** An answer key also needs
> indistinguishable-protein classes, group definitions, database-present-but-absent proteins,
> contaminants, entrapments, and a stated scoring policy. Separate artifact, separate task.

### 6.4 `cleavage_sites.parquet`
`protein_id`, `position`. Where the enzyme *can* cut — needed to recompute yields per condition
without re-running the enumeration. Termini are **not** listed: they are boundaries, not cleavage
events.

### 6.5 `modforms.parquet` *(not yet implemented — `timsim-modify`)*
`modform_id`, `peptide_id`, `mods: list<struct{pos, unimod_id, mass_delta}>`, `mass_monoisotopic`.

**Positional isomers are distinct `modform_id`s** and coexist in one run. v1 cannot do this at
all — the benchmark suite ran **10 separate simulations** (phospho on site A, then site B) and
recovered the site *from the filename*, then fabricated the ground truth as
`SIM["site_probability"] = 1.0`.

> **Necessary but not sufficient for a localisation benchmark.** Also required: site coordinates
> in **protein-residue space** (peptide-relative `mods.pos` breaks for shared peptides and
> repeats — join through §6.3); **condition-specific occupancy**; modified-fragment prediction;
> and **fragment annotations identifying which fragments are site-determining**. Without the last,
> we have an identification benchmark in a localisation costume.

### 6.6 `structure.parquet` *(not yet implemented — `timsim-precursors`)*
Precursors plus predicted properties: **m/z** and the **isotope envelope** (pure chemistry),
`rt_index` ∈ [0,1] *(valid only within the declared LC chemistry family)*, `ccs` in Å² *(per charge
state — conformers differ with charge)*, `charge_propensity`, `ionization_propensity`.

> Note these are **propensities**, not realised values. A peptide's *actual* charge distribution and
> *actual* ionisation response depend on the eluent composition when it elutes and on what co-elutes
> with it — both LC outputs. They are therefore measurements (§8.2), not structure.

**Deliberately absent:** seconds, 1/K₀, peak widths, frame indices, intensities, detectability.

**B12 — three different variances, three columns (or none):** conformer heterogeneity is
*molecular*; predictor uncertainty is *epistemic* (a property of the model, not the ion); device
spread is *measurement*. A single `ccs_sigma` silently sums them and makes the truth
uninterpretable.

### 6.7 `intrinsic_fragment_pattern.parquet` *(not yet implemented)*
Keyed by `(modform_id, charge, activation, collision_energy)`. **Quantity-independent**, so shared
across every sample at a given CE.

> **The name is load-bearing.** *Observed* fragment intensity depends on isolation transmission,
> co-isolates, AGC, and saturation — all mixture-dependent. Reserve "fragment intensity" for
> measured spectra. Same discipline as `amount_amol` stopping at the source.

---

## 7. Modifications — where they live, and why

**"Variable mods" is search-engine vocabulary.** A search engine has variable mods with a
`max_variable_mods` budget because it enumerates candidates. **A sample has site occupancies.**

### 7.0 There are no "fixed" and "variable" mods. There is only occupancy.

Carbamidomethyl-C is not a "fixed mod" at 100% — alkylation is ~95–99% efficient, which is *why
free-cysteine peptides appear in real data*. Set `occupancy = 0.98` and that falls out with no
special case. "Fixed" is occupancy ≈ 1; "variable" is occupancy < 1. **The chemist's parameter is
the alkylation efficiency.**

### 7.1 Only cleavage-**blocking** mods reach back into the digest *(implemented)*

Acetyl-K, ubiquitin-GG, trimethyl-K and TMT-K **physically block trypsin**. They are states of the
protein *before* the protease sees it, so they live in protein-residue coordinates and must be
applied in the digest — not as post-hoc decoration.

```
p_eff(k) = p_cleavage · (1 − blocking_occupancy(k))
```

**No protein-modform enumeration** (which would be 2ⁿ). A blocking mod simply lowers the firing
probability at one site.

**Why this matters scientifically:** a missed cleavage at site *k* now has two causes — the enzyme
failed (`(1−occ)(1−p)`), or the site was modified (`occ`). Same peptide *sequence*, different
*modform*. That coupling is the signature diGly and acetylome experiments rely on: **the missed
cleavage at the modified lysine is how you find it.** A model that decorates peptides after
digestion cannot reproduce this at all.

### 7.2 Everything else is downstream, ordered by *when it happens at the bench*

- **Pre-digestion, protein-level, blocking** → the digest (§7.1).
- **Pre-digestion, protein-level, non-blocking** → phospho-STY, carbamidomethyl-C. Protein-residue
  coordinates; do not change cleavage.
- **Post-digestion, peptide-level** → pyroglutamate (forms on the *peptide* N-terminus, which does
  not exist until after digestion), Met oxidation from handling.

### 7.3 Truncate modforms on **probability mass**, not a count budget

Never enumerate *protein* modforms. Enumerate *peptide* modforms, bounded by sites-per-peptide
(typically ≤4). And do not use `max_variable_mods = 3` — that is the search-engine budget again:
arbitrary, and blind to how much material it discards.

Instead: a **minimum modform abundance floor**, with the omitted mass **measured**. A peptide with
four S/T/Y at 2% occupancy is ~92% unmodified, ~7.8% singly-modified, ~0.24% doubly — a
probability floor captures that automatically; a count budget does not.

Exactly the same discipline as `max_missed_cleavages`: chemist sets occupancy, informatician sets
a floor, tool reports the omitted mass.

### 7.4 No protein sub-population — and this is *provable*

Two models: (a) independent sites at 50% occupancy each; (b) 50% of protein copies "modified" and
carrying *both* sites. Digest, and measure: **peptide 1 is 50% modified and peptide 2 is 50%
modified, in both models.** Identical observables.

Bottom-up proteolysis **destroys the linkage**. Two peptides from one protein copy become two
independent molecules in the tube, and no MS experiment recovers their joint distribution. The
two models are *empirically indistinguishable*. (Top-down would see it. We are not doing top-down.)

**Correlation only survives digestion when the sites land on the same peptide** — where the joint
*is* observable (a doubly-phosphorylated peptide has a distinct mass and fragment pattern), and
where it is biologically real (priming phosphorylation, heptad repeats). Independence
**under-produces multiply-modified peptides**, which are exactly the hard cases for localisation —
so a within-peptide correlation parameter directly tunes benchmark *difficulty*.

**The one case where a protein sub-population is real** — modified copies degraded faster, say — is
a statement about *how much*, not *what exists*. It belongs on the **quantity axis**, as
condition-dependent occupancy. Which is machinery we need anyway.

---

## 8. Measurement-axis stages *(not yet implemented)*

> **Ordering matters, and an earlier revision had it wrong.** Peptides separate **in solution** on
> the column and only *then* reach the ESI needle. So **LC comes before ionisation**, and
> ionisation is *coupled to the LC* — which is precisely why charge and response cannot be
> structural (see below). The chain is:
>
> ```
> load → LC → ionise (source) → ion optics → fragmentation → acquisition → render
> ```

**8.1 LC method — RT and peak shape. First, because everything downstream depends on it.**
`rt_index × gradient → seconds`, valid **only within the declared chemistry family**. Retention
*order* changes with stationary phase, temperature, pH, and ion-pairing. The LC method also
generates the EMG (width, tailing) — peak shape is dominated by column chemistry, particle size,
flow, gradient slope, and dwell volume, **not by the molecule**.

Its outputs are what make the source computable at all: **when** each precursor elutes (hence the
eluent composition it meets), and **who co-elutes with it**.

**8.2 Source (ionisation) — charge is not device-free, and not LC-free either.** The ESI charge
envelope and the ionisation response both depend on:

- **the eluent composition at the moment of elution** — a peptide coming off at 5% ACN ionises
  differently from one coming off at 40%. That is a *function of when it elutes*, i.e. an LC output.
- solvent, pH, additives, flow, emitter geometry, spray voltage, source temperature;
- **what else is in the spray at that instant** (§8.3).

So the structure axis may carry only the **latent propensities** — gas-phase basicity (from basic
residues) and surface activity. The *actual* charge distribution and the *actual* response are
**measurements**, and they are computed here, after the LC.

> This is why the tool once planned as a structure-axis `timsim-ionize` is misconceived. What is
> structural is `timsim-precursors`: enumerate modform × charge, compute **m/z** and the **isotope
> envelope** (pure chemistry, no LC needed), and attach the latent propensities. *Ionisation
> proper* is a measurement stage that consumes the LC output.

#### The placement principle (B14)

> **B14 — Place a stage where the *physics* lives, not where today's *model* needs it.**

v1 uses a simple **binomial charge-state model** that does not need the LC at all. It would
therefore be tempting to leave charge on the structure axis, where it is cheap and cacheable. That
is a trap: ionisation *physically* happens after the column, so if the model is ever coupled —
eluent composition at elution, competition from co-eluting ions — the column would have to **move
between axes**, breaking every artifact and every cached fingerprint.

So the source stage lives on the measurement axis from day one, and the *model* is a switch:

| `charge_model` | Uses | Cost | Status |
|---|---|---|---|
| `binomial` *(default)* | structural propensity only | ~free; output identical across runs | **v2.0 — parity with v1** |
| `coupled` | + eluent composition at elution, + local co-eluting load | needs the LC output | v2.1 (§8.3) |

Turning coupling on later is then a **config change, not a schema migration**. A cheap stage
sitting in an expensive slot costs nothing; a stage in the *wrong* slot costs a migration.

**The run must record which model ran.** A benchmark that reports "12% CV" is uninterpretable
unless it knows whether suppression was in play. This is a provenance field on the run, not an
afterthought.

#### Ionisation efficiency ("flyability") — and the model that would silently break it

v1 **random-samples** ionisation efficiency. Replacing that with a predictor is a real realism
upgrade — but **not every model called "flyability" predicts the same quantity**, and picking the
wrong one corrupts the quantification benchmark in a way that looks fine in every plot.

**The question to ask of any candidate model: was abundance controlled in its training data?**

| Model kind | What it actually predicts | Where it belongs |
|---|---|---|
| **Response factor** — trained on **equimolar** synthetic peptides (ProteomeTools-style) | ionisation efficiency, cleanly separated from abundance | ✅ structural `ionization_propensity` |
| **Observability / detectability** — trained on P(observed) across public repositories | a *marginalised observation*: digestion + LC + source + suppression + instrument + search — **and abundance** | measurement axis, as an empirical **replacement** for the suppression-and-LOD chain |

**Why this matters.** Observability correlates strongly with protein abundance, because abundant
proteins get observed. Use such a model as an intrinsic propensity and you multiply the
**explicitly modelled** `amount_amol` by a factor that **already encodes abundance** — squaring the
abundance effect, inflating the effective dynamic range, and wrecking the LOD behaviour. Composing
it *with* a suppression model double-counts suppression on top of that.

An observability model and an explicit suppression/LOD chain are **alternatives, not composable
factors**. Pick one.

#### Decision for v2.0: keep the random draw as the default, and support both

**Default: `lognormal`** — v1's `generate_normal_efficiency` (median 1e-2, ~1 order of spread,
clipped to [1e-4, 1]), but **identity-keyed** so it is reproducible and A/B-consistent.

The reasoning is not "it's simpler". It is that **a random draw is wrong but *unbiased*, while the
fine-tuned model is plausible but *biased*** — and a bias does not average out, correlates with the
very quantity being measured, and is inherited invisibly by every downstream benchmark. For a
**ground-truth generator**, wrong-and-honest beats plausible-and-tilted.

**Be clear-eyed about what random costs, though.** It does not *add* a bias; it *removes* a real
effect. Real detectability is strongly sequence-biased, so some proteins are systematically hard
(few flying peptides). Randomise it and every protein gets an even shot at clearing a two-peptide
rule — making **protein-level sensitivity unrealistically uniform**, and protein inference **easier
than reality**, which flatters every tool benchmarked against it.

| source | status | rationale |
|---|---|---|
| `lognormal` *(default)* | v1 parity | robust, unbiased per peptide, reproducible, no external dependency |
| `model = "pfly_base"` | the realistic option | equimolar ⇒ a genuine response factor; use the **graded expected response** `Σ P(class)·level`, never `1 − P(non_flyer)` |
| `model = "pfly_finetuned"` | **refused as a structural propensity** | ρ = 0.76 with abundance — would square the abundance effect |

**And run both.** If tool *rankings* move between random and base-PFly flyability, that is a genuine
result — **search-engine comparisons are sensitive to the detectability model** — and nobody has
measured it. If they do not move, the cheap default is vindicated. Either answer is worth having, and
it falls out of the switch for free.

The run **records which source was used**, so a benchmark reporting CVs or sensitivities is
interpretable.

#### PFly: which variant, and why it matters (measured, not speculated)

The candidate model is **PFly** (Wilhelm lab, JPR 2025). It exists in two variants, and **they are
not interchangeable**:

| variant | trained on | target label | abundance controlled? | belongs on |
|---|---|---|---|---|
| **base** | ProteomeTools — >1M **intended-equimolar** synthetic peptides | **MS1-intensity tertiles** (weak/intermediate/strong; non-flyer = synthesised, never seen) | **yes, by construction** | ✅ structural `ionization_propensity` |
| **`pfly_2024_fine_tuned`** | Sinitcyn et al. (PXD024364), 6 human cell lines | *"the number of cell lines from which the peptide was identified"* — i.e. **observability** | **no** | measurement axis only |

**The fine-tuned model's abundance confound is documented by its own authors:**

> *"Peptides originating from low-abundance proteins may be misclassified not due to their ability
> to fly, but rather because of their low concentration in the sample."*

> Spearman **ρ = 0.76** between protein abundance (mean iBAQ) and per-protein precision, *p* < 0.001.

ρ = 0.76 is not a subtle effect. Multiply that prediction by our explicitly-modelled `amount_amol`
and **abundance enters the signal twice** — inflating the dynamic range and moving the LOD cut, while
every plot still looks entirely plausible.

**Koina serves only the fine-tuned variant** (checked against the live registry:
`pfly_2024_fine_tuned`, `pfly_2024_fine_tuned_core`, `pfly_2024_preprocessing` — no base model).
PFly is in **DLOmix** with source freely available, so the ProteomeTools-only weights are the ones to
obtain.

**And use the graded output.** v1 computes `efficiency = 1 - P(non_flyer)`
(`ionization/predictors.py:513`), collapsing four ordered classes to a binary — a strong flyer and a
weak flyer both score ≈ 1.0, discarding the entire dynamic range the model was trained to produce.
The base model's classes **are** MS1-intensity tertiles, so the correct quantity is the **expected
response**: `Σ_k P(class_k) · level_k`.

> v1 has Koina flyability **off by default**, so no published result is affected. But it is a live
> footgun: enabling `use_koina_model` for flyability would silently square the abundance effect.

**The status quo is not neutral either.** Identity-keyed random flyability preserves A/B
consistency (§10), but it randomises *which* peptides are detectable. Real detection is strongly
**sequence-biased**, and that bias governs protein coverage, how many proteins clear a two-peptide
rule, and therefore how hard protein inference really is. A random detectable subset makes protein
inference **easier than reality — flattering every tool benchmarked against it**.

> Note the `plausible_charges` helper used for the S0.5 count is *not* a model — it enumerates
> `1..n` from the basic-residue count purely to bound the universe. The real baseline is v1's
> `BinomialChargeStateDistributionModel` (`ionization/predictors.py:75`), which is what
> `timsim-precursors` should carry forward as the structural propensity.

**8.3 Source response — mixture coupling. The biggest realism gap in v1, which models it not at
all.** Real ESI is competitive: response depends on the *local co-eluting ion load*. This is the
first point where we know **who co-elutes with whom** — which is only knowable *after* the LC, and
is exactly why it cannot live on the structure axis.

**8.4 Ion optics — CCS → 1/K₀.** Mason–Schamp with gas mass and temperature. **CCS depends on ion,
gas, and temperature — not pressure**, because 1/K₀ is a *reduced* mobility, already normalised to
standard number density. Pressure enters the **device model** (TIMS ramp, resolution, scan
calibration), which lives here.

> The transfer function **already exists and is already correctly parameterised**:
> `imspy-core/chemistry/mobility.py:24` → `ccs_to_one_over_k0(ccs, mz, charge, mass_gas=28.013,
> temp=31.85, ...)`, implemented in Rust with a parallel variant. Nobody ever varies the defaults.
> And the predictor **already computes CCS internally** (`ccs/predictors.py:572`) and converts on
> the way out. The physics is right; it is called in the wrong *place*.

**8.5 Fragmentation — not CE-only.** Activation type, collision gas, isolation width, **co-isolation
/ chimeric spectra**, scan timing.

**8.6 Acquisition — DDA needs competition.** Precursor competition and dynamic exclusion make DDA
selection a *stochastic, abundance-dependent* sampling. Multi-run designs are what make DDA
benchmarking meaningful — above all **match-between-runs**, which cannot be evaluated with one file.

**8.7 `run_features.parquet` — the per-run *realised* ground truth**

`rt_apex_s`, `rt_start_s`, `rt_end_s`, `fwhm_s`, `im_apex`, `im_fwhm`, `observed`,
`missingness_cause` (`below_lod` | `not_selected` | `suppressed` | `excluded` | `absent`), and
`feature_correspondence_id` — **the MBR answer key**.

> Apex and FWHM are **not** closed-form functions of EMG `(μ, σ, λ)`; they need numerical
> root-finding. And MBR truth needs *realised* per-run shifts and deformation plus an explicit
> cross-run correspondence artifact. Emitting only nominal LC parameters **relocates** v1's
> reconstruction problem rather than solving it.

---

## 9. The digest *(implemented)*

### 9.1 Two stages, not one

- **`Enumerator`** → *which peptides exist*. Structure. Depends on the enzyme and the truncation
  bound, and nothing else. Shared by every sample.
- **`YieldModel`** → *how much of each*. Quantity. Depends on digestion efficiency and blocking
  occupancy, both of which are condition-dependent.

Measured on the human proteome: **structure 96 ms once; quantity 29 ms per sample.**

### 9.2 The model — analytic expectation

```
p_yield = b_i · b_j · Π_{k internal to (i,j)} (1 − p_eff(k))
   where b_x = 1  if x is a protein terminus   ← NOT a cleavage event
             p_eff(x)  otherwise
```

Amounts sum over **occurrences**, not proteins. Deterministic, **seed-free**, identical at any
thread count.

**What it is not:** an idealised independent-site model. It ignores site-specific kinetics,
protease:substrate ratio and time, enzyme autolysis, sequential cleavage of intermediates,
missed-site correlation, structural accessibility, sequence context beyond the proline rule, and
competition in double digests. **A useful baseline — never "correct molar amounts".**

### 9.3 Mass balance — the invariant that makes it self-checking

Every realisation *partitions* the protein, so every residue lands in exactly one peptide with
probability 1. Therefore, over the complete enumeration:

```
Σ  p_yield(o) · len(o)  =  L        exactly
```

Independence is the **only** assumption, so it survives arbitrary per-site `p_eff` (hence
blocking), and holds at p=0 and p=1. **Verified by Monte Carlo** — 400,000 actually-simulated
digests with non-uniform blocking — which also reproduces mass balance by *counting*.

### 9.4 `--max-missed-cleavages` is a *measured* error bound

```
$ timsim-yield --digestion-efficiency 0.90 ...
  missed cleavages     :  0 → 90.2%  1 → 8.9%  2 → 0.9%   (molar, observed — NOT a parameter)
  truncated at n=2     : omits 0.3502% of residue mass   (measured, not a formula)
    (bounds read from the structure artifact, not re-entered)
  length filter        : omits 23.82% of enumerated mass
```

> Earlier revisions claimed `(1−p)^(n+1)`. That is the *count* probability under an **infinite**
> interior-site model; finite proteins, terminal peptides, and length filters all break it, and
> discarded **mass** ≠ discarded **count**. We enumerate anyway — so measure.

Truncation loss and filter loss are reported **separately**: one is a numerical bound, the other a
modelling choice, and collapsing them would hide both.

### 9.5 Cleavage rules — independently confirmed

`timsim-chem/xcheck/` digests the same FASTA with **Sage** and with us. Human proteome, 20,590
proteins: **identical peptide sets at 0, 1, 2, 3, and 5 missed cleavages. Zero divergence, both
directions.** Cleavage rules and FASTA parsing are derived from Sage (MIT, Michael Lazear);
Sage's `Digest`, decoys, and missed-cleavage budget are deliberately **not** carried over.

---

## 10. Determinism

```rust
seed_for(entity) = hash(global_seed, entity_id)     // identity-keyed, NOT draw-order
```

- Parallel execution is **byte-identical regardless of `--threads`**.
- Adding a protein does not reshuffle other proteins' peptides.
- Raising `n_proteins` 5k → 10k yields a **superset**, not a different sample.

**And it gives A/B consistency for free.** Flyability and charge propensity are *structural*
molecular properties — the same in every sample. They are not *stored and reused*; they are
**rederived and identical**, because they hash the same entity id. Nothing can drift.

> Box–Muller must be fed by a hash with **full avalanche**. Raw FNV under Box–Muller measures
> excess kurtosis +0.091 and σ 1.0076 (real heavy tails); with a SplitMix64 finaliser, +0.005 and
> 0.9995. Guarded by `gauss_is_actually_gaussian`, whose thresholds are calibrated to **fail** the
> unfinalised version. Mean and σ alone would not have caught it — **kurtosis is the load-bearing
> assertion**.

> `imspy_simulation/dispatch.py` already contains a `derive_seed` implementation — 309 lines,
> wired to nothing. It was the right idea, arrived at too early.

---

## 11. Crate layout *(implemented)*

```
timsim-chem      # fasta, enzymes, digest (structure/quantity split), masses, design, ids
timsim-schema    # THE Arrow schemas + validate-on-read. No column name lives elsewhere.
timsim-cli       # timsim-proteome | timsim-digest | timsim-design | timsim-yield
timsim-chem/xcheck/   # independent oracle: our cleavage rules vs Sage's
```

Separate binaries, not `timsim <verb>`: necroflow's cache identity can hash the binary, so one
binary per tool means changing the renderer does not invalidate your digest.

Static binaries, no environment — which also sidesteps necroflow's lack of conda/container support
for these stages, and lets `NodeType.invalidator` hash the binary itself, so cache identity is
exact and automatic rather than a hand-maintained `recipe_identity` string.

---

## 12. Open questions

1. **Candidate-universe size (S0.5)** — the blocking measurement. 2.5M peptides is the base;
   the modform × charge multiplier is unknown. Needs `timsim-modify`.
2. **`amount_amol` — "in sample" or "on column"?** Currently: *in sample at that point in the
   protocol*, with the load stage rescaling so `Σ(nᵢ·Mᵢ)` matches `--load-ng`.
3. **Contaminants** — modelled, or injected via a contaminant FASTA? (Currently: FASTA, tagged.)
4. **How far do we take the LC scope declaration (B1)?** One index per chemistry family, or a
   calibrated retention model?
5. **Does source response (§8.3) ship in v2.0, or v2.1?** Biggest realism gap; genuinely hard.
6. **Does the reference-blank approach survive v2?** Building the `.d` on a real blank means real
   peptides leak in and are scored as false positives — hence five copies of a blank-subtraction
   workaround in the benchmark repo. **A scientific-validity issue, not ergonomics.**
