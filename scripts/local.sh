#!/bin/bash
# Log of what was run locally while LPC is down (stand-in for the snakemake rules). Newest at the bottom.
#   bash scripts/local.sh        (from the repo root; reruns everything)
set -euo pipefail

# ---- 2026-11-16  rules/filter.smk: merge (--both) + whitelist
# 1. keep calls made by both callers (data/merged.csv already passes filter / minad / strand)
python3 - <<'PY'
import csv
with open('data/merged.csv',newline='') as fi, open('data/merged_both.csv','w',newline='') as fo:
    r=csv.DictReader(fi)
    w=csv.DictWriter(fo,fieldnames=r.fieldnames)
    w.writeheader()
    w.writerows(x for x in r if x['callers']=='Mutect2;VarDict')
PY

# 2. CHIP whitelist + hotspot flags; protected rows also on their own
python3 scripts/whitelist.py \
    -i data/merged_both.csv \
    -o data/whitelist.csv \
    --protected data/whitelist_protected.csv \
    --stats data/whitelist_stats.tsv \
    --whitelist-dir whitelist \
    --hotspots whitelist/20261116_hotspot_rules.tsv

# ---- 2026-11-16  recurrence exploration (interactive in RStudio; writes data/recurrence_protected_all_samples_{calls,vars}.csv, data/recurrence{,_removed}_{calls,vars}.csv, figures/recurrence/)
# Rscript scripts/recurrence.R

# ---- 2026-11-16  germline + functional filters (interactive in RStudio; reads data/recurrence_calls.csv, writes data/{germline,functional}{,_removed}_{calls,vars}.csv + final sets data/{somatic_main,somatic_secondary,ch}_{calls,vars}.csv, figures/germline/)
# Rscript scripts/germline.R

# ---- 2026-11-16  somatic main set vs clinical status (interactive in RStudio; reads data/20261116_somatic_main_{calls_review,vars}.xlsx + ss/20261014_grey_master.xlsx, writes data/somatic_main_clinical_{calls,samples,genes,vars}.csv, figures/clinical/)
# Rscript scripts/clinical.R

# ---- 2026-11-16  coverage overview (reads the two data/029-Gray_PreCancer*metrics_summary files, writes data/metrics_coverage.csv, figures/metrics/)
# Rscript scripts/metrics.R
