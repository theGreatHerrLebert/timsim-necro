# P1.3 — Phospho + FLR (site-localization benchmark)

Can a search engine put the phosphate on the **right residue**? A distinct eval axis: not "is the peptide
right" (recall) or "is the ID real" (FDP), but "is the SITE right" — scored by the **false-localization rate
(FLR)**. And we improve on v1: simulate both **positional isomers in ONE run** (co-eluting), instead of v1's
separate-runs-recovered-by-filename hack.

## State — the sim can already do this (the good news)

- **Positional isomers are generated.** `enumerate_modforms` emits one modform per site combination, so
  phospho@S3 and phospho@T5 are two distinct modforms of the same peptide (confirmed: 1652 peptides get ≥2
  phospho modforms on a 2000-peptide HYE digest with `mods.toml`). Each modform carries `mod_positions` +
  `mod_names` — the **ground-truth site**.
- **Fragment m/z is SITE-SPECIFIC.** `annotate()` inserts `[UNIMOD:21]` at the phospho residue, so mscore's
  `calculate_product_ion_series` shifts every b-ion ≥ the site and every y-ion covering it by +79.966. Two
  isomers therefore have **distinguishable fragments** — localization is possible in the current pipeline.
- **Intensity is backbone-only (mod-blind):** the predictor sees the plain sequence (`frag_input` =
  `[precursor_id, sequence, charge]`), so phospho changes to fragmentation *propensity* are not modelled.
  For LOCALIZATION this is acceptable (the evidence is the site-specific *m/z*, which is exact); it only
  costs intensity realism. Upgrade path: Koina AlphaPeptDeep (phospho-aware) for the intensity axis — noted,
  not required for the first cut.
- **Ground-truth site** for a precursor: `precursor_id → modform_id` (precursors) →
  `mod_positions` where `mod_names == "Phospho"` (modforms).

## The work

### 1. A phospho-focused config
`mods_phospho.toml`: Phospho (UNIMOD 21, STY) at an enriched occupancy, **drop GG/Oxidation** so peptides
get clean single-/few-phospho isomers (mods.toml's diGly+phospho combos muddy the isomer set). Keep
Carbamidomethyl (fixed). Consider restricting the benchmark to peptides with **≥2 STY sites** (localization
is trivial with one site).

### 2. DiaNN phospho search (the one real unknown)
Search with variable Phospho on STY + localization output. DiaNN reports the localized `Modified.Sequence`
and a per-site confidence; **verify the exact columns on a real phospho report** (`PTM.Site.Confidence` /
`PTM.Q.Value` / localized `Modified.Sequence` positions) — this is the only thing not yet confirmed and
gates the scorer's parse. Likely flags: `--var-mod UniMod:21,79.966331,STY --monitor-mod UniMod:21`
(+ `--relaxed-prot-inf`); confirm against the installed DiaNN 2.5.

### 3. FLR scorer — `timsim_eval/v2_flr_eval.py` (NEW)
- **Truth:** phospho precursors → true site(s) (modform `mod_positions`), keyed on the backbone sequence.
- **Observed:** DiaNN's localized modified sequence → localized site(s) + confidence.
- **Metrics:**
  - **FLR** = among confidently-localized single-phospho PSMs (confidence ≥ τ), the fraction whose site ≠
    the true site. Report an **FLR-vs-confidence curve** (the honest form) + FLR at DiaNN's reported cutoff.
  - **Localization recall** = fraction of true phosphosites confidently + correctly localized.
  - **Isomer discrimination** (the v1 improvement): for peptides with ≥2 co-eluting isomers, does the engine
    recover each site, or collapse them / pick the abundant one? Report per-isomer recovery + a confusion
    over sites.
- **Eligibility/coverage tables** (as in quant): identified / localizable (≥2 STY) / confidently-localized.

### 4. Flow
`mods_phospho.toml` + a phospho `search` variant (var-mods phospho + localization flags) + a `score_flr`
node → `v2_flr_eval`. Gate on a `--phospho` flag. Bruker DIA first (mirrors the other closed pipelines).

## Review resolutions (codex)

- **Two DISTINCT tasks, don't conflate rows with molecules.**
  1. **PRIMARY = isolated single-isomer empirical FLR** (the calibration metric). Restrict to features where
     exactly ONE truth isomer occupies the precursor/RT region. `FLR(τ) = #{accepted calls, wrong site} /
     #{accepted calls}` where "accepted" = **localization** confidence ≥ τ (NOT the peptide/precursor
     q-value). It is **empirical FLR vs simulator truth**, not the engine's self-estimated FLR.
  2. **SECONDARY = co-eluting-isomer component recovery** (a separate mixture-resolution benchmark, NOT
     folded into FLR). Truth is at the **component level** (each isomer: site, abundance, RT profile).
     One-to-one match each truth isomer → a reported localized feature (site + backbone + charge + RT
     tolerance); missing = FN, wrong-site = FP; report a site confusion. If the engine emits a single call
     for an unresolved mixture, score vs the **dominant-abundance** isomer and LABEL it "dominant-isomer
     classification," never general FLR. **Do NOT score a row against "the isomer it was quantified as"** —
     that's circular unless provenance is retained independently.
- **Always pair FLR with coverage.** Report **FLR-vs-threshold AND correct-localization-recall-vs-threshold**
  curves + one pre-registered operating point (accepted calls at an estimated 1% localization FLR). Without
  recall, an engine scores FLR 0 by declining to localize.
- **Noiseless is a SANITY UPPER-BOUND, not the benchmark.** Exact m/z makes localization near-trivial except
  where site-determining ions don't bracket the competing sites / are thresholded out. Run a **factorial**:
  {isolated, co-eluting} × {noiseless, A1/A2 noise}. Hardness comes from missing/again-thresholded
  site-determining ions, isomer overlap (adjacent sites → fewer discriminating cleavages), co-elution
  degree + relative abundance, resolution/tolerance, interference. Specify what A1/A2 perturb (m/z error /
  intensity / background / fragment interference).
- **Mod-blind intensity is a real caveat, not just "intensity realism."** It changes which site-determining
  ions cross detection thresholds → can bias FLR either way. Call this a **controlled m/z-only localization
  benchmark**; Koina AlphaPeptDeep is the upgrade (validate its site-specific fragment encoding first).
- **Eligibility = ≥2 CANDIDATE phospho sites under the SEARCH RULES** (not just ≥2 STY): account for
  blocked/fixed-modified residues and termini. Exact **site-set equality** for the single-phospho primary;
  keep **multi-phospho separate** (site-level precision/recall + peptide-level exact-site-set accuracy).
- **Separate ID from localization.** A wrong backbone/charge is an ID error (FDP), NOT a localization error —
  exclude it from the FLR denominator. Handle unlocalized/ambiguous/duplicate/decoy/missing-confidence/ties/
  no-site-determining-ions explicitly.
- **Truth key = component, not backbone.** Include charge, precursor/isolation context, run, RT feature, and
  the simulator's precursor/modform id — the same backbone recurs across distinct features.
- **DiaNN localization confidence is UNVERIFIED — the concrete blocker.** Do not define FLR around columns
  that may not exist. Run one small phospho search first; if DiaNN 2.5 exposes a valid per-site confidence,
  use it; else treat it as a tool-specific adapter (rank/assignment) problem.

## Original open questions (resolved above)
1. **Is noiseless site-specific m/z a FAIR FLR benchmark, or too easy?** The sim's fragment m/z is exact, so
   with no interference a search should localize perfectly (FLR≈0) except where isomer fragments genuinely
   overlap (shared ions). Does a meaningful FLR *require* A1/A2 noise + co-elution (ties to the realism
   track), or is the isomer m/z overlap alone enough to stress localization? Likely: run the FLR curve both
   noiseless and with noise, and expect noise to be where FLR becomes non-trivial.
2. **FLR definition granularity:** per-PSM site-level vs per-site aggregated; single-site only vs multi-site
   peptides; how to score a peptide with 2 true sites where the engine gets 1 right, 1 wrong.
3. **Co-eluting isomers:** when both isomers are present, DiaNN reports one localized form per precursor
   feature — is the right estimand "did the reported form match the more abundant isomer" or a mixture
   deconvolution? Start simple: score the reported localization against the isomer it was quantified as.
4. **Mod-blind intensity:** does backbone-only intensity make localization artificially easy or hard? The
   site evidence is m/z (exact), so intensity mostly affects which *ions are observed above threshold*.
   Note as a known limitation; the Koina AlphaPeptDeep upgrade is the fix if it matters.
5. **DiaNN localization columns** — the concrete blocker; resolve by running one small phospho search.

## Validation
- `mods_phospho.toml` produces clean phospho isomers (peptides with ≥2 single-phospho modforms).
- A rendered phospho run searched by DiaNN yields localized sites; the FLR scorer recovers a low FLR on a
  noiseless render (sanity) and a realistic FLR under A1/A2 noise.
- Deterministic given the same render + report; unit-test the FLR math on synthetic (truth site, localized
  site, confidence) triples.
