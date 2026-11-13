"""
CHIP example on gray cohort VEP csvs (scripts copied from 00-pmbb/somatic).

    snakemake -s scripts/chip/chip.smk --configfile config/chip.yaml -j 2

Per caller -> {outdir}/{caller}/:
  chip_prep        chip_example_prep.py: calls.csv.gz + sites.csv + sample_qc.csv (all PASS)
  chip_seq_filter  chip_seq_filtering.py: CH genes, MANE, orientation, DP / alt depth / VAF
  chip_whitelist   chip_whitelist.py: wl_* flags
  chip_recurrence  chip_recurrence.py: carriers per variant (whitelist + manual review only)
  chip_artifact    chip_artifact.R: review table (no batch / age test on gray)
"""

OUT = config["outdir"]
S = "scripts/chip"


def opt(c, k):
    return config["callers"][c].get(k, config[k])


localrules: chip_all, chip_prep, chip_seq_filter, chip_whitelist, chip_recurrence, chip_artifact

rule chip_all:
    input:
        expand(f"{OUT}/{{caller}}/chip_artifact.csv", caller=config["callers"]),


rule chip_prep:
    input:
        csv=lambda w: config["callers"][w.caller]["input"],
        script=f"{S}/chip_example_prep.py",
    output:
        calls=f"{OUT}/{{caller}}/calls.csv.gz",
        sites=f"{OUT}/{{caller}}/sites.csv",
        qc=f"{OUT}/{{caller}}/sample_qc.csv",
    log:
        f"{OUT}/{{caller}}/logs/chip_prep.log",
    shell:
        "python {input.script} -i {input.csv} -o {OUT}/{wildcards.caller} > {log} 2>&1"


rule chip_seq_filter:
    input:
        calls=rules.chip_prep.output.calls,
        sites=rules.chip_prep.output.sites,
        script=f"{S}/chip_seq_filtering.py",
    output:
        csv=f"{OUT}/{{caller}}/chip_seq.csv",
        stats=f"{OUT}/{{caller}}/chip_seq_stats.tsv",
    log:
        f"{OUT}/{{caller}}/logs/chip_seq_filter.log",
    params:
        orient=lambda w: "" if opt(w.caller, "require_orientation") else "--no-orientation",
    shell:
        r"""
        python {input.script} --calls {input.calls} --sites {input.sites} --chrom all \
            --whitelist-dir {config[whitelist_dir]} --min-dp {config[min_dp]} \
            --min-alt-depth {config[min_alt_depth]} --min-vaf {config[min_vaf]} {params.orient} \
            -o {output.csv} --stats {output.stats} > {log} 2>&1
        """


rule chip_whitelist:
    input:
        csv=rules.chip_seq_filter.output.csv,
        script=f"{S}/chip_whitelist.py",
    output:
        f"{OUT}/{{caller}}/chip_calls.csv",
    log:
        f"{OUT}/{{caller}}/logs/chip_whitelist.log",
    shell:
        "python {input.script} -i {input.csv} -o {output} --whitelist-dir {config[whitelist_dir]} > {log} 2>&1"


rule chip_recurrence:
    input:
        calls=rules.chip_whitelist.output[0],
        qc=rules.chip_prep.output.qc,
        script=f"{S}/chip_recurrence.py",
    output:
        f"{OUT}/{{caller}}/chip_recurrence.csv",
    log:
        f"{OUT}/{{caller}}/logs/chip_recurrence.log",
    shell:
        r"""
        python {input.script} -i {input.calls} --sample-qc {input.qc} \
            --min-carriers {config[min_carriers]} -o {output} > {log} 2>&1
        """


rule chip_artifact:
    input:
        calls=rules.chip_whitelist.output[0],
        rec=rules.chip_recurrence.output[0],
        qc=rules.chip_prep.output.qc,
        script=f"{S}/chip_artifact.R",
    output:
        variants=f"{OUT}/{{caller}}/chip_artifact.csv",
        calls=f"{OUT}/{{caller}}/chip_calls.artifact.csv",
    log:
        f"{OUT}/{{caller}}/logs/chip_artifact.log",
    params:
        exempt=",".join(config.get("exempt_genes") or []),
    shell:
        r"""
        Rscript {input.script} --calls {input.calls} --recurrence {input.rec} --sample-qc {input.qc} \
            --p-age {config[p_age]} --p-batch {config[p_batch]} --exempt "{params.exempt}" \
            --out-variants {output.variants} --out-calls {output.calls} > {log} 2>&1
        """
