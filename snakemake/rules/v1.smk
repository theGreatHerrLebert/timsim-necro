# ── MEASUREMENT: Bruker .d via the monolithic v1 `timsim` seam (the default) ──
# The strangler seam: v1's LC / ion mobility / fragmentation / acquisition, driven from v2 artifacts
# via the `--v2-*` flags. This backend owns DDA. The v1 config TOML is an INPUT (necroflow: RawData
# `invalidator = hashes_file("config")`), so editing v1.toml restages the run.
#
# NOTE: the config file is passed as `{input.v1_config}`, NOT `{config[...]}` — a bare `{config}` in a
# Snakemake shell block would interpolate the whole Snakemake config dict, not our v1 TOML path.
rule simulate:
    input:
        v1_config          = V1_CONFIG,
        proteome           = f"{STRUCT}/proteome.parquet",
        peptides           = f"{STRUCT}/peptides.parquet",
        occurrences        = f"{STRUCT}/peptide_occurrences.parquet",
        peptide_quantities = f"{QUANT}/peptide_quantities.parquet",
        precursors         = f"{STRUCT}/precursors.parquet",
        precursor_ccs      = f"{STRUCT}/precursor_ccs.parquet",
        peptide_rt         = f"{STRUCT}/peptide_rt.parquet",
    output:
        raw = directory(f"{OUT}/v1/{{sample}}/raw"),
    threads: 2
    resources:
        mem_mb = 8192,
    shell:
        "mkdir -p {output.raw} && timsim {input.v1_config} --save-path {output.raw} "
        "--v2-proteome {input.proteome} --v2-peptides {input.peptides} "
        "--v2-occurrences {input.occurrences} "
        "--v2-peptide-quantities {input.peptide_quantities} "
        "--v2-precursors {input.precursors} --v2-ccs {input.precursor_ccs} "
        "--v2-rt {input.peptide_rt} --v2-sample {wildcards.sample} --seed {config[seed]}"
