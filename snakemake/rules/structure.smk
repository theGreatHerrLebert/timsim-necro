# ── STRUCTURE axis + shared spectral prediction ──────────────────────────────
# Which molecules exist, and their instrument-independent properties. NONE of these rules carry a
# `{sample}` wildcard, so each runs exactly once and is shared by every sample and every backend —
# the necroflow "structure deduplicated by fingerprint" line, but here it is just the graph shape.

# proteome: FASTAs → proteins. The spec is an INPUT (necroflow: `invalidator = hashes_file("spec")`),
# so editing hye.toml restages the proteome and everything below it.
rule proteome:
    input:
        spec = PROTEOME_SPEC,
    output:
        f"{STRUCT}/proteome.parquet",
    shell:
        "{config[timsim_bin]}/timsim-proteome --spec {input.spec} --out {output}"


# digest: proteins → peptides. Three co-outputs of ONE call (one computation), exactly as the
# necroflow rule returns a 3-tuple.
rule digest:
    input:
        proteome = f"{STRUCT}/proteome.parquet",
    output:
        peptides       = f"{STRUCT}/peptides.parquet",
        occurrences    = f"{STRUCT}/peptide_occurrences.parquet",
        cleavage_sites = f"{STRUCT}/cleavage_sites.parquet",
    shell:
        "{config[timsim_bin]}/timsim-digest --proteome {input.proteome} "
        "--out-peptides {output.peptides} --out-occurrences {output.occurrences} "
        "--out-cleavage-sites {output.cleavage_sites} "
        "--max-missed-cleavages {config[max_missed_cleavages]} "
        "--min-length {config[min_length]} --max-length {config[max_length]}"


# modify: peptides → modforms + the modification spec as an artifact. The mods TOML is an INPUT
# (necroflow: `hashes_file("mods")`). modifications.parquet is read by BOTH modify and yield, so the
# two can never disagree about an occupancy.
rule modify:
    input:
        peptides = f"{STRUCT}/peptides.parquet",
        mods     = MODS_SPEC,
    output:
        modforms      = f"{STRUCT}/modforms.parquet",
        modifications = f"{STRUCT}/modifications.parquet",
    shell:
        "{config[timsim_bin]}/timsim-modify --peptides {input.peptides} --mods {input.mods} "
        "--out-modforms {output.modforms} --out-modifications {output.modifications} "
        "--floor {config[floor]}"


# precursors: modforms → ions (m/z, isotope envelopes, charge fractions). STRUCTURE.
rule precursors:
    input:
        peptides = f"{STRUCT}/peptides.parquet",
        modforms = f"{STRUCT}/modforms.parquet",
    output:
        f"{STRUCT}/precursors.parquet",
    shell:
        "{config[timsim_bin]}/timsim-precursors --peptides {input.peptides} "
        "--modforms {input.modforms} --out {output} "
        "--charge-model {config[charge_model]} --seed {config[seed]}"


# ccs: precursors → collision cross section. STRUCTURE (instrument-independent). Python tool from PATH.
# Consumed by the v1 and Bruker (ion-mobility) backends; not by Thermo/SCIEX (no ion mobility).
rule ccs:
    input:
        precursors = f"{STRUCT}/precursors.parquet",
        peptides   = f"{STRUCT}/peptides.parquet",
    output:
        f"{STRUCT}/precursor_ccs.parquet",
    shell:
        "timsim-ccs --precursors {input.precursors} --peptides {input.peptides} --out {output}"


# rt: peptides → retention-time index. STRUCTURE; deep model (Chronologer by default). Python tool.
rule rt:
    input:
        peptides = f"{STRUCT}/peptides.parquet",
    output:
        f"{STRUCT}/peptide_rt.parquet",
    shell:
        "timsim-rt --peptides {input.peptides} --out {output}"


# ── shared spectral prediction (Thermo + Bruker measurement branches) ─────────
# frag_input is STRUCTURE (no CE/model); fragments is MEASUREMENT but instrument-shared — one
# `frag_model` per run, so both are wildcard-free and predicted once for every sample.

# frag_input: precursors + modforms → the frozen (precursor_id, [UNIMOD]-sequence, charge) input.
rule frag_input:
    input:
        precursors    = f"{STRUCT}/precursors.parquet",
        peptides      = f"{STRUCT}/peptides.parquet",
        modforms      = f"{STRUCT}/modforms.parquet",
        modifications = f"{STRUCT}/modifications.parquet",
    output:
        f"{STRUCT}/fragment_prediction_input.parquet",
    shell:
        "{config[timsim_bin]}/timsim-frag-input --precursors {input.precursors} "
        "--peptides {input.peptides} --modforms {input.modforms} "
        "--modifications {input.modifications} --out {output}"


# fragments: annotated input → predicted fragment intensities. `frag_model` is quoted so the default
# empty string (local timsTOF) passes as a valid empty argument, matching the necroflow command.
rule fragments:
    input:
        frag_input = f"{STRUCT}/fragment_prediction_input.parquet",
    output:
        directory(f"{STRUCT}/fragment_intensities"),
    shell:
        "timsim-fragments --precursors {input.frag_input} "
        "--collision-energy {config[collision_energy]} --model '{config[frag_model]}' --out {output}"


# spectra: precursors + fragments → instrument-independent MS1 isotope + MS2 fragment spectra.
rule spectra:
    input:
        precursors           = f"{STRUCT}/precursors.parquet",
        peptides             = f"{STRUCT}/peptides.parquet",
        modforms             = f"{STRUCT}/modforms.parquet",
        modifications        = f"{STRUCT}/modifications.parquet",
        fragment_intensities = f"{STRUCT}/fragment_intensities",
    output:
        directory(f"{STRUCT}/ion_spectra"),
    shell:
        "{config[timsim_bin]}/timsim-spectra --precursors {input.precursors} "
        "--peptides {input.peptides} --modforms {input.modforms} "
        "--modifications {input.modifications} "
        "--fragment-intensities {input.fragment_intensities} --out {output}"
