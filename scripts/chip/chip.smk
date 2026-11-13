"""
CHIP example on gray cohort VEP csvs (scripts copied from 00-pmbb/somatic).

    snakemake -s scripts/chip/chip.smk --configfile config/chip.yaml -j 2

Per caller -> {outdir}/{caller}/:
  chip_seq_filter  chip_seq_filtering.py: CH genes, MANE, DP / alt depth / VAF; + sample_qc.csv (all PASS)
                   input = basic_filtered.csv (or raw vep csv), samples_from = raw csv for the sample N
  chip_whitelist   chip_whitelist.py: wl_* flags
  chip_recurrence  chip_recurrence.py: carriers per variant (whitelist + manual review only)
  chip_artifact    chip_artifact.R: batch / age / germline flags -> review table
"""

OUT = config["outdir"]
S = "scripts/chip"


localrules: chip_all, chip_seq_filter, chip_whitelist, chip_recurrence, chip_artifact

rule chip_all:
    input:
        expand(f"{OUT}/{{caller}}/chip_artifact.csv", caller=config["callers"]),


rule chip_seq_filter:
    input:
        csv=lambda w: config["callers"][w.caller]["input"],
        samples=lambda w: config["callers"][w.caller].get("samples_from") or [],
        script=f"{S}/chip_seq_filtering.py",
    output:
        csv=f"{OUT}/{{caller}}/chip_seq.csv",
        stats=f"{OUT}/{{caller}}/chip_seq_stats.tsv",
        qc=f"{OUT}/{{caller}}/sample_qc.csv",
    log:
        f"{OUT}/{{caller}}/logs/chip_seq_filter.log",
    params:
        samples=lambda w, input: f"--samples-from {input.samples}" if input.samples else "",
        drop=",".join(config.get("drop_cols") or []),
    shell:
        r"""
        python {input.script} -i {input.csv} {params.samples} --sample-qc-out {output.qc} \
            --drop-cols "{params.drop}" \
            --whitelist-dir {config[whitelist_dir]} --min-dp {config[min_dp]} \
            --min-alt-depth {config[min_alt_depth]} --min-vaf {config[min_vaf]} \
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
        qc=rules.chip_seq_filter.output.qc,
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
        qc=rules.chip_seq_filter.output.qc,
        script=f"{S}/chip_artifact.R",
    output:
        variants=f"{OUT}/{{caller}}/chip_artifact.csv",
        calls=f"{OUT}/{{caller}}/chip_calls.artifact.csv",
    log:
        f"{OUT}/{{caller}}/logs/chip_artifact.log",
    params:
        exempt=",".join(config.get("exempt_genes") or []),
        batch_test=config.get("batch_test", True),
        age_test=config.get("age_test", True),
        germline_test=config.get("germline_test", True),
    shell:
        r"""
        Rscript {input.script} --calls {input.calls} --recurrence {input.rec} --sample-qc {input.qc} \
            --p-age {config[p_age]} --p-batch {config[p_batch]} \
            --batch-test {params.batch_test} --age-test {params.age_test} \
            --germline-test {params.germline_test} --p-binom {config[p_binom]} --gnomad-max-af {config[gnomad_max_af]} \
            --exempt "{params.exempt}" \
            --out-variants {output.variants} --out-calls {output.calls} > {log} 2>&1
        """
