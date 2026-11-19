# Gray Pre-Cancer somatic pipeline

Workflow for calling somatic variants on the Gray Pre-Cancer cohort (119 samples, IDT xGen Target HYB PreCancer panel) with two tumor-only callers, Mutect2 and VarDict, and curating them into a final callset.

The final calls are split into three sets:

1. **CH**: calls that match the CHIP whitelist.
2. **Somatic main**: pathogenic / likely pathogenic by AutoGVP or ClinVar, or oncogenic by OncoKB.
3. **Somatic secondary**: everything else that passes filtering.

| Stage | Rule file / script | Status |
|---|---|---|
| 1. Alignment + QC | `rules/alignment.smk` | written |
| 2a. Calling: Mutect2 | `rules/mutect2.smk` | written |
| 2b. Calling: VarDict | `rules/vardict.smk` | written |
| 3. Annotation | `rules/annotate.smk` | written (run once per caller) |
| 4. Filtering: sequence → merge → whitelist | `rules/filter.smk` | written |
| 5. Recurrence | `scripts/recurrence.R` | exploratory (RStudio); cutoff not final |
| 6. Germline, functional + final sets | `scripts/germline.R` | exploratory (RStudio) |
| 7. Downstream analysis | `scripts/clinical.R`, `metrics.R`, `top_genes.R` | exploratory (RStudio) |

### Things to clarify
- Recurrence cutoff (fraction of samples) for stage 5.
- Germline VAF cutoff (currently `VAF_MAX = 0.35`) and whether high-VAF CH clones should be exempt beyond the whitelist.
- Panel of normals: the current run is no-PON (`pon.use_pon: false`, label `no_pon`). Decide whether to build or obtain one.
- Stages 5–6 should become Snakemake rules once the cutoffs are settled. `scripts/recurrence.py` is from the older `somatic_filter` pipeline and still imports `chip_whitelist`, which no longer exists.
- While LPC is down, `scripts/local.sh` records what was run locally in place of the Snakemake rules.

## Layout

```
config/      config.yaml (all stages), fastq.yml (alignment inputs)
rules/       Snakemake rule files, one per stage
scripts/     helper scripts (Python for the pipeline, R for exploration + reports)
input/       sample sheets, panel gene list, protected variants
whitelist/   CHIP gene list, variant whitelists, rules tsv, hotspot rules (the newest file by date is used)
profiles/    LSF profile
```

## Running

Every rule file reads `config/config.yaml` and is run on its own, in stage order:

```bash
snakemake -s rules/alignment.smk --configfile config/config.yaml --profile profiles/lsf
snakemake -s rules/mutect2.smk   --configfile config/config.yaml --profile profiles/lsf
snakemake -s rules/vardict.smk   --configfile config/config.yaml --profile profiles/lsf
snakemake -s rules/annotate.smk  --configfile config/config.yaml --config caller=mutect2
snakemake -s rules/annotate.smk  --configfile config/config.yaml --config caller=vardict
snakemake -s rules/filter.smk    --configfile config/config.yaml -j 2
```

Dry-run first with `-n` after any config change.

---

## 1. Alignment + QC (`alignment.smk`)

Inputs come from `config/fastq.yml`, which lists every run / lane / index per sample.

- `fastqc` on every FASTQ.
- `bwa-mem2 mem` per lane, then add read groups (ID = run.lane, PU = run.lane.index, SM = sample).
- Merge lanes, `MarkDuplicates`.
- BQSR: `BaseRecalibrator` on the targets (+50 bp) using known indels, dbSNP and Mills / 1000G; then `ApplyBQSR`.
- `ValidateSamFile`, then the final BAM per sample. `input/bam.table` lists them all.

QC on the final BAM, collected into one `multiqc` report:

| Check | Tool |
|---|---|
| Alignment summary | `CollectAlignmentSummaryMetrics` |
| Insert size | `CollectInsertSizeMetrics` |
| Hybrid-capture metrics (on-target rate, coverage) | `CollectHsMetrics` with the bait + target interval lists |
| Coverage per target, at ≥ 50 / 100 / 250 / 500 / 1000x | `mosdepth`, after deduplication and on the raw reads |

## 2a. Mutect2 (`mutect2.smk`)

This stage follows the GATK somatic short-variant best practices in tumor-only mode, unscattered, on the targets (+50 bp).

**Calling.** Mutect2 with gnomAD as the germline resource, plus a panel of normals if `pon.use_pon` is true. It writes F1R2 counts for the orientation model.

**Panel of normals (optional).** When `use_pon` is true and no `existing_pon_vcf` is given, each normal is called (`--max-mnp-distance 0`) and filtered. Then GenomicsDBImport and CreateSomaticPanelOfNormals build the panel.

**Orientation artifact model.** LearnReadOrientationModel learns how often each base change is supported by only one read orientation, the signature of damage artifacts.

**Contamination.** GetPileupSummaries at common germline SNPs, then CalculateContamination for a contamination estimate plus tumor segmentation.

**Filtering.** FilterMutectCalls with the contamination table, segmentation and orientation priors. Output: `mutect2/calls/{label}/{sample}/{sample}.{label}.filtered.vcf.gz`, where `label` comes from `analysis.label` so runs with and without a PON don't overwrite each other.

## 2b. VarDict (`vardict.smk`)

VarDictJava per sample on the targets (+50 bp), down to `vardict_af_thr` (0.001), then `teststrandbias.R` and `var2vcf_valid.pl`. Output: one VCF per sample.

## 3. Annotation (`annotate.smk`)

Run once per caller (`--config caller=mutect2` or `caller=vardict`).

- Drop VarDict's symbolic structural-variant calls (`<DEL>`, `<DUP>`, `<INV>`), then `bcftools norm` (split multiallelics, left-align), set ID = `CHROM_POS_REF_ALT` and sort.
  - The SV calls go straight into one table for manual review instead.
- VEP 116 (`--everything --canonical --mane`) with REVEL, SpliceAI, AlphaMissense, MaveDB, LOFTEE and ClinVar + AutoGVP.
- `vep_vcf_parser.py` (cohort mode) flattens each VCF to CSV, one row per sample and variant. It adds per-strand alt read counts (`Sample.AltFwd` / `Sample.AltRev`, from Mutect2 `SB` or VarDict `ALD`) and the raw FORMAT string.
- Samples are concatenated into `{caller}/all_samples.{caller}[.{label}].vep.csv`, then annotated with OncoKB (`oncokb.tumor_type: BRCA`).

## 4. Filtering (`filter.smk`)

### 4a. Sequence filter (`sequence.py`, per caller)

Removed rows go to `sequence_removed.csv` with `removed_at` set to the step that dropped them.

| Step | Keep if |
|---|---|
| s01 | gene is on the panel gene list, `Bio.type == protein_coding`, AltFrac / Depth / AltDepth present, and on target |
| s02 | `Sample.Depth ≥ min_dp` (100) |
| s03 | REF and ALT ≤ `max_len` bp (100) |
| s04 | Mutect2 only: alt reads on both orientations (F1R2 ≥ 1 and F2R1 ≥ 1) |

Flags added (not removed here):

| Flag | Meaning |
|---|---|
| `pass_minad` | alt depth ≥ `min_ad_snv` (5) for SNVs, ≥ `min_ad_indel` (10) for indels |
| `pass_strand` | alt forward fraction within [`strand_min`, 1 − `strand_min`] (0.1–0.9); NA without strand counts |

Step counts go to `reports/filter/sequence_stats.{caller}.tsv`.

### 4b. Merge (`merge.py`)

One row per call (sample + variant + gene) across both callers.

- Kept when called by both callers (`both_callers: true`) and at least one caller passes each of:
  - `pass_filter`: Mutect2 has no FILTER tag apart from PASS and `mutect2_ignore_tags` (germline, clustered_events, multiallelic); VarDict FILTER is PASS.
  - `pass_minad` and `pass_strand` from 4a.
- Depth, alt depth and VAF are averaged over the callers. Per-caller values are kept as `Mutect2.*` / `VarDict.*`, and annotations come from Mutect2 when it has the call.
- Everything else goes to `merged_removed.csv`.

### 4c. Whitelist (`whitelist.py`)

Matches each call against the files in `whitelist/` and adds flags. Nothing is removed.

| Flag | Meaning |
|---|---|
| `wl_mis` / `wl_lof` / `wl_splice` | matches the missense, nonsense/frameshift or splice whitelist |
| `wl_exception` | matches a gene-specific rule (exon or protein range) |
| `wl_rule` | which rule(s) from `*_chip_whitelist_rules.tsv` matched |
| `wl_transcript` | whether the MANE transcript matches the gene list's accession (`.` when there's none) |
| `wl_whitelist` | CHIP candidate |
| `wl_manualreview` | possible candidate that needs manual review |
| `wl_hotspot` / `wl_hotspot_rule` | matches a tumor hotspot rule (`*_hotspot_rules.tsv`: gene + variant class + protein range) |
| `protected` | `wl_whitelist`, `wl_manualreview` or `wl_hotspot`; exempt from the recurrence and germline filters |

Protected rows are also written to `whitelist_protected.csv`.

### Reports

Not part of the default target. They need pandoc ≥ 2.8:

```bash
snakemake -s rules/filter.smk --configfile config/config.yaml -j 1 sequence_report merge_report
```

- `sequence_report.html`: what each sequence step removes, per caller.
- `merge_report.html`: caller overlap. Also writes `merged_mutect2_only.csv` and `merged_vardict_only.csv`.

## 5. Recurrence (`recurrence.R`, exploratory)

Picks a cutoff for removing variants shared by many samples, which are more likely artifacts than real clones.

- sample frequency = samples carrying the exact variant / all 119 samples.
- Variants found in every sample are removed, even when protected.
- Compares candidate cutoffs (rows and variants removed, protected rows, median VAF), and how recurrence changes when only calls above an alt-depth cutoff count as carriers.
- Output: `data/recurrence_calls.csv` for stage 6.

## 6. Germline, functional + final sets (`germline.R`, exploratory)

| Step | Removed if | Protected exempt? |
|---|---|---|
| Germline | gnomAD max AF > 0.005 (exome + genome; overall and AFR, AMR, EAS, MID, NFE, SAS, excluding the bottlenecked / mixed groups), or VAF > 0.35 | yes |
| Functional | no MANE Select transcript, or LOFTEE LoF level 3 / 4 | no |

The remaining calls are split into the final sets:

- `ch`: CHIP whitelist calls.
- `somatic_main`: AutoGVP or ClinVar pathogenic / likely pathogenic, or OncoKB oncogenic / likely oncogenic.
- `somatic_secondary`: all others.

## 7. Downstream analysis (exploratory)

- `clinical.R`: somatic main set against clinical status (BRCA carrier status, cancer-associated collection). Denominators are all 119 sequenced samples.
- `metrics.R`: mean target coverage per sample, with outliers flagged, from the two sequencing batches' HsMetrics summaries.
- `top_genes.R`: samples with a call, per gene, in the somatic main set.

---
