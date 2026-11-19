"""
Somatic filtering on the cohort VEP csvs (rules/annotate.smk merge_samples).

    snakemake -s rules/filter.smk --configfile config/config.yaml -j 2

Per caller -> {outdir}/{caller}/:
  sequence  sequence.py: gene / on target / depth / indel length / alt orientation; pass_minad + pass_strand flags
All callers -> {outdir}/:
  merge     merge.py: one row per call, per-caller columns; kept if called by both callers (both_callers) and >= 1 caller passes filter, minad and strand
  whitelist whitelist.py: CHIP whitelist (wl_*) + tumour hotspots (wl_hotspot) -> protected; protected rows also in whitelist_protected.csv
Stats -> {reports}/: sequence_stats.{caller}.tsv, merged_stats.tsv, whitelist_stats.tsv
Report (not part of the default target; needs pandoc >= 2.8, e.g. export RSTUDIO_PANDOC=~/.conda/envs/pandoc/bin):
    snakemake -s rules/filter.smk --configfile config/config.yaml -j 1 sequence_report merge_report
  sequence_report  sequence_report.Rmd -> {reports}/sequence_report.html
  merge_report  merge_report.Rmd -> {reports}/merge_report.html + {outdir}/merged_{mutect2,vardict}_only.csv
"""

F = config["filter"]
OUT = F["outdir"]
REP = F["reports"]
TARGETS_BED = config["resources"]["targets_bed"]


localrules: filter_all, sequence, merge, whitelist, merge_report, sequence_report

rule filter_all:
    input:
        f"{OUT}/whitelist.csv",


rule sequence:
    input:
        csv=lambda w: F["callers"][w.caller],
        script="scripts/sequence.py",
    output:
        csv=f"{OUT}/{{caller}}/sequence.csv",
        removed=f"{OUT}/{{caller}}/sequence_removed.csv",
        stats=f"{REP}/sequence_stats.{{caller}}.tsv",
    log:
        f"{OUT}/{{caller}}/logs/sequence.log",
    shell:
        r"""
        python {input.script} -i {input.csv} -o {output.csv} --removed {output.removed} --stats {output.stats} \
            --genes "{F[genes]}" --targets-bed "{TARGETS_BED}" \
            --min-dp {F[min_dp]} --max-len {F[max_len]} \
            --min-ad-snv {F[min_ad_snv]} --min-ad-indel {F[min_ad_indel]} --strand-min {F[strand_min]} > {log} 2>&1
        """


rule merge:
    input:
        csv=expand(f"{OUT}/{{caller}}/sequence.csv", caller=F["merge"].values()),
        script="scripts/merge.py",
    output:
        csv=f"{OUT}/merged.csv",
        removed=f"{OUT}/merged_removed.csv",
        stats=f"{REP}/merged_stats.tsv",
    log:
        f"{OUT}/logs/merge.log",
    params:
        callers=" ".join(f"--{k} {OUT}/{v}/sequence.csv" for k, v in F["merge"].items()),
        ignore=",".join(F["mutect2_ignore_tags"]),
        both="--both" if F.get("both_callers") else "",
    shell:
        r"""
        python {input.script} {params.callers} --mutect2-ignore "{params.ignore}" {params.both} \
            -o {output.csv} --removed {output.removed} --stats {output.stats} > {log} 2>&1
        """


rule whitelist:
    input:
        csv=rules.merge.output.csv,
        hotspots=F["hotspots"],
        script="scripts/whitelist.py",
    output:
        csv=f"{OUT}/whitelist.csv",
        protected=f"{OUT}/whitelist_protected.csv",
        stats=f"{REP}/whitelist_stats.tsv",
    log:
        f"{OUT}/logs/whitelist.log",
    shell:
        r"""
        python {input.script} -i {input.csv} -o {output.csv} --protected {output.protected} \
            --stats {output.stats} --whitelist-dir {F[whitelist_dir]} --hotspots {input.hotspots} > {log} 2>&1
        """


rule merge_report:
    input:
        merged=rules.merge.output.csv,
        removed=rules.merge.output.removed,
        rmd="scripts/merge_report.Rmd",
    output:
        html=f"{REP}/merge_report.html",
        m2=f"{OUT}/merged_mutect2_only.csv",
        vd=f"{OUT}/merged_vardict_only.csv",
    log:
        f"{OUT}/logs/merge_report.log",
    shell:
        r"""
        Rscript -e 'rmarkdown::render("{input.rmd}", output_file = "'$PWD'/{output.html}",
            intermediates_dir = "'$PWD'/{REP}", knit_root_dir = "'$PWD'",
            params = list(merged = "'$PWD'/{input.merged}", removed = "'$PWD'/{input.removed}",
                          out_dir = "'$PWD'/{OUT}"))' > {log} 2>&1
        """


rule sequence_report:
    input:
        csv=expand(f"{OUT}/{{caller}}/sequence.csv", caller=F["callers"]),
        removed=expand(f"{OUT}/{{caller}}/sequence_removed.csv", caller=F["callers"]),
        rmd="scripts/sequence_report.Rmd",
    output:
        html=f"{REP}/sequence_report.html",
    log:
        f"{OUT}/logs/sequence_report.log",
    params:
        callers=",".join(F["callers"]),
    shell:
        r"""
        Rscript -e 'rmarkdown::render("{input.rmd}", output_file = "'$PWD'/{output.html}",
            intermediates_dir = "'$PWD'/{REP}", knit_root_dir = "'$PWD'",
            params = list(dir = "'$PWD'/{OUT}", callers = "{params.callers}", min_ad_snv = {F[min_ad_snv]},
                          min_ad_indel = {F[min_ad_indel]}, strand_min = {F[strand_min]}))' > {log} 2>&1
        """
