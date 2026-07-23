# ── QUANTITY axis ─────────────────────────────────────────────────────────────
# The design (the mixture) and the peptide amounts it implies. In the necroflow flow BOTH of these
# run once for the WHOLE design — the per-sample column lives inside protein/peptide_quantities, and
# the fan-out only happens later when the measurement node selects a sample with `--v2-sample`. So
# these rules are wildcard-free too: one design → one set of quantity tables shared by every render.

# design: proteome × design spec → the whole-experiment mixture. Four co-outputs of one call. The
# design TOML is an INPUT (necroflow: `hashes_file("spec")`), and note the seed lives INSIDE the
# spec, not on the command line — it is a property of the experiment.
rule design:
    input:
        proteome = f"{STRUCT}/proteome.parquet",
        spec     = DESIGN_SPEC,
    output:
        samples            = f"{QUANT}/samples.parquet",
        runs               = f"{QUANT}/runs.parquet",
        sample_run_map     = f"{QUANT}/sample_run_map.parquet",
        protein_quantities = f"{QUANT}/protein_quantities.parquet",
    shell:
        "{config[timsim_bin]}/timsim-design --proteome {input.proteome} --spec {input.spec} "
        "--out-samples {output.samples} --out-runs {output.runs} "
        "--out-sample-run-map {output.sample_run_map} "
        "--out-protein-quantities {output.protein_quantities}"


# peptide_yield: structure × mixture → peptide amounts (amol), per sample, in one table. Takes
# `modifications` so a blocking mod (acetyl-K, GG-K, TMT-K) actually stops the protease.
rule peptide_yield:
    input:
        proteome           = f"{STRUCT}/proteome.parquet",
        occurrences        = f"{STRUCT}/peptide_occurrences.parquet",
        cleavage_sites     = f"{STRUCT}/cleavage_sites.parquet",
        protein_quantities = f"{QUANT}/protein_quantities.parquet",
        modifications      = f"{STRUCT}/modifications.parquet",
    output:
        f"{QUANT}/peptide_quantities.parquet",
    shell:
        "{config[timsim_bin]}/timsim-yield --proteome {input.proteome} "
        "--occurrences {input.occurrences} --cleavage-sites {input.cleavage_sites} "
        "--protein-quantities {input.protein_quantities} --modifications {input.modifications} "
        "--digestion-efficiency {config[digestion_efficiency]} --out {output}"
