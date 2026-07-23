# ── MEASUREMENT: Thermo .raw (no ion mobility) + phase 2 ─────────────────────
# The fan-out. `{sample}` selects a column from the shared peptide_quantities table (`--sample`), so
# this is the ONLY thing that differs between samples — everything upstream is computed once.

# render_thermo: author the feature space into a real Thermo .raw template. Three co-outputs of one
# command — the .raw, its per-precursor answer key, and a durable run manifest — so they can never
# drift from the data. Restages when the template file changes (declared as an input).
rule render_thermo:
    input:
        precursors         = f"{STRUCT}/precursors.parquet",
        peptide_rt         = f"{STRUCT}/peptide_rt.parquet",
        ion_spectra        = f"{STRUCT}/ion_spectra",
        peptide_quantities = f"{QUANT}/peptide_quantities.parquet",
        template           = config["template"] or [],
    output:
        data_raw = f"{OUT}/thermo/{{sample}}/data.raw",
        truth    = f"{OUT}/thermo/{{sample}}/truth.parquet",
        manifest = f"{OUT}/thermo/{{sample}}/manifest.json",
    threads: 2
    resources:
        mem_mb = 8192,
    shell:
        "{config[timsim_bin]}/timsim-render-thermo --precursors {input.precursors} "
        "--peptide-rt {input.peptide_rt} --ion-spectra {input.ion_spectra} "
        "--peptide-quantities {input.peptide_quantities} --sample {wildcards.sample} "
        "--template {input.template} --intensity-scale {config[intensity_scale]} "
        "--frag-model '{config[frag_model]}' --method {config[method]} "
        "--expected-ce {config[collision_energy]} "
        "--out {output.data_raw} --thermo-truth {output.truth} --manifest {output.manifest}"


# search: DiaNN library-free over the .raw. Reads .raw natively via the .NET 8 runtime (DOTNET_ROOT +
# on PATH). The FASTA is a content-hashed dependency (declared as an input).
rule search:
    input:
        data_raw = f"{OUT}/thermo/{{sample}}/data.raw",
        fasta    = SEARCH_FASTA if SEARCH_FASTA else [],
    output:
        diann = directory(f"{OUT}/thermo/{{sample}}/diann"),
    threads: 16
    resources:
        mem_mb = 32768,
    shell:
        "mkdir -p {output.diann} && "
        "DOTNET_ROOT={config[dotnet_root]} PATH={config[dotnet_root]}:$PATH {config[diann]} "
        "--f {input.data_raw} --fasta {input.fasta} --out {output.diann}/report.parquet "
        "--fasta-search --predictor --gen-spec-lib --qvalue {config[qvalue]} "
        "--threads {config[search_threads]} --met-excision --cut 'K*,R*' "
        "--missed-cleavages {config[max_missed_cleavages]} "
        "--min-pep-len {config[min_length]} --max-pep-len {config[max_length]} "
        "--var-mods 1 --unimod35 --reanalyse --relaxed-prot-inf"


# score: the DiaNN report against the render's answer key. Content-addressed to the exact .raw/truth/
# DB that produced it — the number the whole simulate→search→score run exists to produce.
rule score:
    input:
        diann    = f"{OUT}/thermo/{{sample}}/diann",
        truth    = f"{OUT}/thermo/{{sample}}/truth.parquet",
        peptides = f"{STRUCT}/peptides.parquet",
    output:
        f"{OUT}/thermo/{{sample}}/metrics.json",
    shell:
        "python -m timsim_eval.v2_thermo_eval "
        "--report {input.diann}/report.parquet --truth {input.truth} "
        "--peptides {input.peptides} --fdr {config[qvalue]} --out {output}"
