# ── MEASUREMENT: lean v2 Bruker .d (WITH ion mobility) + phase 2 ─────────────
# The imspy-free projector (`timsim-render`): it reuses the SAME frag_input→fragments→spectra chain as
# the Thermo path, then PROJECTS the instrument-independent spectra onto a reference .d's DIA
# schedule. Bruker has ion mobility, so — unlike Thermo — it consumes precursor_ccs (CCS → 1/K0).

# The reference .d is a directory; its cheap analysis.tdf metadata DB carries the frame schedule / DIA
# windows / calibration the render reads. Declaring it as an input restages the render when the
# reference changes (necroflow: `hashes_reference_d` hashes exactly this file).
REFERENCE_TDF = os.path.join(config["reference_d"], "analysis.tdf") if config.get("reference_d") else []

rule render:
    input:
        precursors         = f"{STRUCT}/precursors.parquet",
        peptide_rt         = f"{STRUCT}/peptide_rt.parquet",
        ion_spectra        = f"{STRUCT}/ion_spectra",
        precursor_ccs      = f"{STRUCT}/precursor_ccs.parquet",
        peptide_quantities = f"{QUANT}/peptide_quantities.parquet",
        reference_tdf      = REFERENCE_TDF,
    output:
        raw   = directory(f"{OUT}/bruker/{{sample}}/data.d"),
        truth = f"{OUT}/bruker/{{sample}}/truth.parquet",
    threads: 2
    resources:
        mem_mb = 8192,
    shell:
        "{config[timsim_bin]}/timsim-render --precursors {input.precursors} "
        "--peptide-rt {input.peptide_rt} --ion-spectra {input.ion_spectra} "
        "--precursor-ccs {input.precursor_ccs} "
        "--peptide-quantities {input.peptide_quantities} --sample {wildcards.sample} "
        "--reference-d {config[reference_d]} --dia --intensity-scale {config[intensity_scale]} "
        "--out {output.raw} --truth {output.truth}"


# search_bruker: DiaNN library-free over the .d. Unlike the Thermo .raw, a Bruker .d is DiaNN's NATIVE
# input on Linux — no .NET runtime needed.
rule search_bruker:
    input:
        data_d = f"{OUT}/bruker/{{sample}}/data.d",
        fasta  = SEARCH_FASTA if SEARCH_FASTA else [],
    output:
        diann = directory(f"{OUT}/bruker/{{sample}}/diann"),
    threads: 16
    resources:
        mem_mb = 32768,
    shell:
        "mkdir -p {output.diann} && {config[diann]} "
        "--f {input.data_d} --fasta {input.fasta} --out {output.diann}/report.parquet "
        "--fasta-search --predictor --gen-spec-lib --qvalue {config[qvalue]} "
        "--threads {config[search_threads]} --met-excision --cut 'K*,R*' "
        "--missed-cleavages {config[max_missed_cleavages]} "
        "--min-pep-len {config[min_length]} --max-pep-len {config[max_length]} "
        "--var-mods 1 --unimod35"


# score_bruker: identical scorer to `score`, keyed on the Bruker answer key. `timsim_eval` is
# instrument-agnostic — the truth schema is the same 8 columns — so the .d closes search→score
# exactly like the Thermo path.
rule score_bruker:
    input:
        diann    = f"{OUT}/bruker/{{sample}}/diann",
        truth    = f"{OUT}/bruker/{{sample}}/truth.parquet",
        peptides = f"{STRUCT}/peptides.parquet",
    output:
        f"{OUT}/bruker/{{sample}}/metrics.json",
    shell:
        "python -m timsim_eval.v2_thermo_eval "
        "--report {input.diann}/report.parquet --truth {input.truth} "
        "--peptides {input.peptides} --fdr {config[qvalue]} --out {output}"
