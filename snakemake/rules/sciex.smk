# ── MEASUREMENT: SCIEX ZenoTOF SWATH → open mzML (lean v2 projector) + phase 2 ─
# The imspy-free counterpart of the old v1 `.wiff` build: it reuses the SAME frag_input→fragments→
# spectra feature-space chain as the Thermo/Bruker paths, then projects the instrument-independent
# `ion_spectra` onto a SYNTHESISED fixed-width SWATH schedule and writes open mzML — no `.wiff`, no
# template file. The SWATH params live in the command (config.yaml), so there is no config TOML and no
# hidden template dependency to hash. No CCS (SCIEX has no ion mobility). Co-emits the per-precursor
# answer key so a DiaNN search closes search→score exactly like Thermo/Bruker.

rule render_sciex:
    input:
        precursors         = f"{STRUCT}/precursors.parquet",
        peptide_rt         = f"{STRUCT}/peptide_rt.parquet",
        ion_spectra        = f"{STRUCT}/ion_spectra",
        peptide_quantities = f"{QUANT}/peptide_quantities.parquet",
    output:
        mzml  = f"{OUT}/sciex/{{sample}}/sciex.mzML",   # a single mzML FILE (not the v1 directory)
        truth = f"{OUT}/sciex/{{sample}}/truth.parquet",
    threads: 2
    resources:
        mem_mb = 8192,
    shell:
        "{config[timsim_bin]}/timsim-render-sciex --precursors {input.precursors} "
        "--peptide-rt {input.peptide_rt} --ion-spectra {input.ion_spectra} "
        "--peptide-quantities {input.peptide_quantities} --sample {wildcards.sample} "
        "--gradient-length-s {config[sciex_gradient_s]} --cycle-time-s {config[sciex_cycle_s]} "
        "--mz-min {config[sciex_mz_min]} --mz-max {config[sciex_mz_max]} "
        "--window-width {config[sciex_window_width]} "
        "--collision-energy {config[collision_energy]} --intensity-scale {config[intensity_scale]} "
        "--frag-model '{config[frag_model]}' --out {output.mzml} --truth {output.truth}"


# search_sciex: DiaNN library-free over the open mzML — DiaNN's native open input, no .NET, no vendor SDK.
rule search_sciex:
    input:
        mzml  = f"{OUT}/sciex/{{sample}}/sciex.mzML",
        fasta = SEARCH_FASTA if SEARCH_FASTA else [],
    output:
        diann = directory(f"{OUT}/sciex/{{sample}}/diann"),
    threads: 16
    resources:
        mem_mb = 32768,
    shell:
        "mkdir -p {output.diann} && {config[diann]} "
        "--f {input.mzml} --fasta {input.fasta} --out {output.diann}/report.parquet "
        "--fasta-search --predictor --gen-spec-lib --qvalue {config[qvalue]} "
        "--threads {config[search_threads]} --met-excision --cut 'K*,R*' "
        "--missed-cleavages {config[max_missed_cleavages]} "
        "--min-pep-len {config[min_length]} --max-pep-len {config[max_length]} "
        "--var-mods 1 --unimod35"


# score_sciex: the same instrument-agnostic timsim_eval scorer, keyed on the SWATH answer key.
rule score_sciex:
    input:
        diann    = f"{OUT}/sciex/{{sample}}/diann",
        truth    = f"{OUT}/sciex/{{sample}}/truth.parquet",
        peptides = f"{STRUCT}/peptides.parquet",
    output:
        f"{OUT}/sciex/{{sample}}/metrics.json",
    shell:
        "python -m timsim_eval.v2_thermo_eval "
        "--report {input.diann}/report.parquet --truth {input.truth} "
        "--peptides {input.peptides} --fdr {config[qvalue]} --out {output}"
