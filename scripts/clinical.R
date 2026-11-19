# Somatic main set vs clinical status: burden, genes and VAF by BRCA carrier status (and cancer-associated collection).
# Reads the Excel copies of the final set + the master sheet; denominators = all 119 sequenced samples, not only the 45 with a call.
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(data.table)
library(here)

TOP_GENES <- 12         # genes shown in the per-gene plots
LOW_VAF   <- 0.01
VAF_TIER  <- c(0, 0.01, 0.05, 0.10, 0.35, Inf)
theme_set(theme_bw(base_size = 12))

FIG_DIR <- here("figures/clinical")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(name, width = 8, height = 4.5) ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)

# ---- read
calls   <- read_excel(here("data/20261116_somatic_main_calls_review.xlsx"), guess_max = 5000)
vars    <- read_excel(here("data/20261116_somatic_main_vars.xlsx"), guess_max = 5000)
ss      <- read_excel(here("ss/20261014_grey_master.xlsx"), guess_max = 5000)
samples <- fread(here("input/gray_sample_ids_bams.tsv"), header = FALSE)$V1      # the 119 sequenced samples
nrow(calls); nrow(vars); length(samples)
# 73 54 119

# columns used from the master sheet; stops here if one is named differently
SS_COLS <- c("Sample_ID", "Genetic_Mutation", "Cancer_Associated_Collection", "Age", "Mean_Target_Coverage", "Case_Anaysis", "Germline_Status")
stopifnot(SS_COLS %in% names(ss))

# BRCA status from Genetic_Mutation (searched anywhere in the string, so a BRCA VUS would also land here -- check the table below)
ss <- ss %>%
  select(all_of(SS_COLS)) %>%
  mutate(Age = suppressWarnings(as.numeric(Age)), Mean_Target_Coverage = suppressWarnings(as.numeric(Mean_Target_Coverage)),
         BRCA_Status = case_when(is.na(Genetic_Mutation) ~ "Unknown",
                                 grepl("BRCA1", Genetic_Mutation) & grepl("BRCA2", Genetic_Mutation) ~ "BRCA1+BRCA2",
                                 grepl("BRCA1", Genetic_Mutation) ~ "BRCA1",
                                 grepl("BRCA2", Genetic_Mutation) ~ "BRCA2",
                                 grepl("^Normal", Genetic_Mutation) ~ "Non-carrier",
                                 grepl("^Not tested", Genetic_Mutation) ~ "Not tested",
                                 TRUE ~ "Other"),
         BRCA_Carrier = case_when(startsWith(BRCA_Status, "BRCA") ~ "Carrier",
                                  BRCA_Status == "Non-carrier" ~ "Non-carrier"))      # NA = other gene / not tested / unknown
ss %>% count(BRCA_Status, Genetic_Mutation) %>% print(n = Inf)

ss$Mean_Target_Coverage

table(ss$Case_Anaysis)
# Case Control 
# 63      56 

x <- ss |> filter(Germline_Status == "Pathogenic_Carrier")
x$Genetic_Mutation

# ---- join
sum(duplicated(ss$Sample_ID))                 # > 0 would duplicate calls in the join
setdiff(samples, ss$Sample_ID)                # sequenced, not in the master sheet
setdiff(calls$Sample.ID, samples)

calls <- calls %>%
  left_join(ss, by = c("Sample.ID" = "Sample_ID")) %>%
  mutate(vaf_tier = cut(Sample.AltFrac, VAF_TIER, right = FALSE, labels = c("< 1%", "1-5%", "5-10%", "10-35%", ">= 35%")))
nrow(calls)                                   # still 73
sum(is.na(calls$BRCA_Status))                 # calls with no master sheet row

# one row per sequenced sample, 0 calls if none
per_sample <- tibble(Sample_ID = samples) %>%
  left_join(ss, by = "Sample_ID") %>%
  left_join(calls %>% group_by(Sample_ID = Sample.ID) %>%
              summarise(n_calls = n(), n_genes = n_distinct(Gene), max_vaf = max(Sample.AltFrac), median_vaf = median(Sample.AltFrac)),
            by = "Sample_ID") %>%
  mutate(n_calls = coalesce(n_calls, 0L), n_genes = coalesce(n_genes, 0L), has_call = n_calls > 0)

fwrite(calls %>% arrange(Sample.ID, ID), here("data/somatic_main_clinical_calls.csv"))
fwrite(per_sample %>% arrange(Sample_ID), here("data/somatic_main_clinical_samples.csv"))

# review sheet: one row per call with the clinical columns, + a PIK3CA / TP53 subset
review <- calls %>%
  select(Sample.ID, Genetic_Mutation, Age, Cancer_Associated_Collection, Gene, HGVSc, HGVSp, Sample.Depth, Sample.AltDepth, Sample.AltFrac, OncoKB.Oncogenic, ClinVar.SIG) %>%
  arrange(Sample.ID, Gene)
writexl::write_xlsx(list(all = review, PIK3CA_TP53 = filter(review, Gene %in% c("PIK3CA", "TP53")) %>% arrange(Gene, HGVSp, Sample.ID)),
                    here("data/20261116_somatic_main_clinical_review.xlsx"))
nrow(review); count(filter(review, Gene %in% c("PIK3CA", "TP53")), Gene)
# 73; PIK3CA 16, TP53 1

# ---- 0. mean target coverage; outliers = outside 1.5 x IQR of the cohort (the boxplot rule)
q <- quantile(per_sample$Mean_Target_Coverage, c(0.25, 0.5, 0.75), na.rm = TRUE)
COV_LOW  <- q[[1]] - 1.5 * (q[[3]] - q[[1]])
COV_HIGH <- q[[3]] + 1.5 * (q[[3]] - q[[1]])
per_sample <- per_sample %>%
  mutate(cov_outlier = case_when(is.na(Mean_Target_Coverage) ~ NA_character_,
                                 Mean_Target_Coverage < COV_LOW ~ "low", Mean_Target_Coverage > COV_HIGH ~ "high", TRUE ~ "no"))
summary(per_sample$Mean_Target_Coverage)
c(LOW = COV_LOW, HIGH = COV_HIGH)
per_sample %>% filter(cov_outlier %in% c("low", "high")) %>%
  select(Sample_ID, cov_outlier, Mean_Target_Coverage, BRCA_Status, n_calls) %>% arrange(Mean_Target_Coverage)
cor.test(per_sample$Mean_Target_Coverage, per_sample$n_calls, method = "spearman", exact = FALSE)      # deeper = more calls?

fig <- per_sample %>% filter(!is.na(Mean_Target_Coverage)) %>%
  ggplot(aes(reorder(Sample_ID, Mean_Target_Coverage), Mean_Target_Coverage)) +
  geom_hline(yintercept = q[[2]], linetype = "dashed", color = "grey50") +
  geom_point(aes(color = cov_outlier != "no")) +
  geom_text_repel(aes(label = if_else(cov_outlier == "no", "", Sample_ID)), size = 3, min.segment.length = 0) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "red"), guide = "none") +
  labs(title = "Mean target coverage per sample", subtitle = paste0("Median: ", round(q[[2]]), "x; Outliers (1.5*IQR)"),
       x = paste0("Samples"), y = "Mean target coverage") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
save_fig("00_mean_target_coverage", width = 6, height = 4.5)

# ---- 1. burden by status
# samples, samples with >= 1 call, calls per sample, VAF of the calls
burden <- function(status) per_sample %>%
  group_by(status = .data[[status]]) %>%
  summarise(n_samples = n(), n_with_call = sum(has_call), pct_with_call = round(100 * mean(has_call), 1),
            n_calls = sum(n_calls), mean_calls = round(mean(n_calls), 2), max_calls = max(n_calls),
            median_max_vaf = round(median(max_vaf, na.rm = TRUE), 4))

burden("BRCA_Status")
burden("BRCA_Carrier")
burden("Cancer_Associated_Collection")
per_sample %>% count(BRCA_Status, Cancer_Associated_Collection) %>%
  pivot_wider(names_from = Cancer_Associated_Collection, values_from = n, values_fill = 0)       # are the two confounded?

# carrier vs non-carrier (NA = other / not tested, dropped by both tests)
fisher.test(table(per_sample$BRCA_Carrier, per_sample$has_call))
wilcox.test(n_calls ~ BRCA_Carrier, data = per_sample)
fisher.test(table(per_sample$Cancer_Associated_Collection, per_sample$has_call))
kruskal.test(n_calls ~ Cancer_Associated_Collection, data = per_sample)

# age: older = more calls?
cor.test(per_sample$Age, per_sample$n_calls, method = "spearman", exact = FALSE)
per_sample %>% group_by(has_call) %>% summarise(n = n(), n_age = sum(!is.na(Age)), median_age = median(Age, na.rm = TRUE))

fig <- per_sample %>% count(BRCA_Status, has_call) %>%
  ggplot(aes(BRCA_Status, n, fill = has_call)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "steelblue"), labels = c("No call", ">= 1 call")) +
  labs(title = "Samples with a somatic main call by BRCA status", x = NULL, y = "Samples", fill = NULL) +
  theme(legend.position = "top")
save_fig("01_samples_with_call_by_brca", width = 6, height = 4.5)

fig <- ggplot(per_sample, aes(BRCA_Status, n_calls)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, height = 0.05, alpha = 0.5) +
  labs(title = "Somatic main calls per sample by BRCA status", x = NULL, y = "Calls per sample")
save_fig("02_calls_per_sample_by_brca", width = 6, height = 4.5)

fig <- ggplot(filter(per_sample, !is.na(Age)), aes(Age, n_calls, color = BRCA_Status)) +
  geom_jitter(width = 0, height = 0.08, alpha = 0.7) +
  labs(title = "Somatic main calls per sample by age", x = "Age", y = "Calls per sample", color = NULL)
save_fig("03_calls_per_sample_by_age", width = 7, height = 4.5)

# ---- 2. genes
genes <- calls %>%
  group_by(Gene) %>%
  summarise(n_calls = n(), n_variants = n_distinct(ID), n_samples = n_distinct(Sample.ID),
            median_vaf = round(median(Sample.AltFrac), 4), min_vaf = round(min(Sample.AltFrac), 4), max_vaf = round(max(Sample.AltFrac), 4)) %>%
  arrange(desc(n_samples), desc(n_calls))
print(genes, n = Inf)
fwrite(genes, here("data/somatic_main_clinical_genes.csv"))
top_genes <- head(genes$Gene, TOP_GENES)

# samples with a call in the gene, per status, as a count and as % of that status's samples
gene_by <- function(status) {
  denom <- per_sample %>% count(status = .data[[status]], name = "n_status")
  calls %>%
    distinct(Gene, Sample.ID, status = .data[[status]]) %>%
    count(Gene, status, name = "n_samples") %>%
    complete(Gene, status = denom$status, fill = list(n_samples = 0)) %>%
    left_join(denom, by = "status") %>%
    mutate(pct = round(100 * n_samples / n_status, 1), Gene = factor(Gene, levels = genes$Gene)) %>%
    arrange(Gene, status)
}
gene_brca <- gene_by("BRCA_Status")
gene_brca %>% mutate(cell = paste0(n_samples, " (", pct, "%)"), status = paste0(status, " n=", n_status)) %>%
  select(Gene, status, cell) %>% pivot_wider(names_from = status, values_from = cell) %>% print(n = Inf)
gene_by("Cancer_Associated_Collection") %>% mutate(cell = paste0(n_samples, " (", pct, "%)"), status = paste0(status, " n=", n_status)) %>%
  select(Gene, status, cell) %>% pivot_wider(names_from = status, values_from = cell) %>% print(n = Inf)

# per gene, carrier vs non-carrier: Fisher on samples with / without a call in the gene (small counts, not corrected)
n_carrier <- table(per_sample$BRCA_Carrier)
gene_fisher <- calls %>%
  filter(!is.na(BRCA_Carrier)) %>%
  distinct(Gene, Sample.ID, BRCA_Carrier) %>%
  count(Gene, BRCA_Carrier) %>%
  pivot_wider(names_from = BRCA_Carrier, values_from = n, values_fill = 0) %>%
  bind_rows(tibble(Carrier = integer(), `Non-carrier` = integer())) %>%      # both columns exist even if one group has no calls
  mutate(across(c(Carrier, `Non-carrier`), ~ coalesce(.x, 0L))) %>%
  rowwise() %>%
  mutate(p = fisher.test(matrix(c(Carrier, n_carrier[["Carrier"]] - Carrier,
                                  `Non-carrier`, n_carrier[["Non-carrier"]] - `Non-carrier`), 2))$p.value) %>%
  ungroup() %>% arrange(p)
n_carrier
print(gene_fisher, n = Inf)

fig <- gene_brca %>% filter(Gene %in% top_genes) %>%
  ggplot(aes(Gene, pct, fill = status)) +
  geom_col(position = position_dodge(preserve = "single")) +
  labs(title = "Top genes by BRCA status", subtitle = "% of each status's samples with a call in the gene", x = NULL, y = "% of samples", fill = NULL) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
save_fig("04_top_genes_by_brca", width = 9, height = 5)

fig <- calls %>% filter(Gene %in% top_genes) %>%
  distinct(Gene, Sample.ID, BRCA_Status) %>%
  mutate(Gene = factor(Gene, levels = rev(top_genes))) %>%
  ggplot(aes(Gene, fill = BRCA_Status)) +
  geom_bar() + coord_flip() +
  labs(title = "Top genes", x = NULL, y = "Samples with a call", fill = NULL) +
  theme(legend.position = "top")
save_fig("05_top_genes_n_samples", width = 7, height = 5)

# ---- 3. VAF
summary(calls$Sample.AltFrac)
calls %>% count(vaf_tier)
calls %>% count(BRCA_Status, vaf_tier) %>% pivot_wider(names_from = vaf_tier, values_from = n, values_fill = 0)
calls %>% group_by(BRCA_Status) %>%
  summarise(n_calls = n(), median_vaf = median(Sample.AltFrac), min_vaf = min(Sample.AltFrac), max_vaf = max(Sample.AltFrac),
            n_low = sum(Sample.AltFrac < LOW_VAF))

# calls >= 35% VAF are only here because they are protected; in a carrier's own gene they may be the germline variant
calls %>% filter(Sample.AltFrac >= 0.35) %>% count(Gene, BRCA_Status, germline)

# per sample (max VAF), since calls from one sample are not independent
wilcox.test(max_vaf ~ BRCA_Carrier, data = per_sample)

fig <- ggplot(calls, aes(Sample.AltFrac, fill = BRCA_Status)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = LOW_VAF, linetype = "dashed") +
  scale_x_log10() +
  labs(title = "VAF of the somatic main calls", subtitle = paste0("Line = VAF ", LOW_VAF), x = "VAF (log)", y = "Calls", fill = NULL) +
  theme(legend.position = "top")
save_fig("06_vaf_by_brca_hist", width = 7, height = 4.5)

fig <- ggplot(calls, aes(BRCA_Status, Sample.AltFrac)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, alpha = 0.5) +
  geom_hline(yintercept = LOW_VAF, linetype = "dashed") +
  scale_y_log10() +
  labs(title = "VAF by BRCA status", x = NULL, y = "VAF (log)")
save_fig("07_vaf_by_brca_box", width = 6, height = 4.5)

fig <- calls %>% filter(Gene %in% top_genes) %>%
  mutate(Gene = factor(Gene, levels = top_genes)) %>%
  ggplot(aes(Gene, Sample.AltFrac, color = BRCA_Status)) +
  geom_jitter(width = 0.2, alpha = 0.8) +
  geom_hline(yintercept = LOW_VAF, linetype = "dashed") +
  scale_y_log10() +
  labs(title = "VAF by gene (top genes)", x = NULL, y = "VAF (log)", color = NULL) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
save_fig("08_vaf_by_gene", width = 9, height = 5)

# ---- 4. variant level: who carries each variant
var_status <- calls %>%
  group_by(ID) %>%
  summarise(n_carrier = sum(BRCA_Carrier == "Carrier", na.rm = TRUE), n_noncarrier = sum(BRCA_Carrier == "Non-carrier", na.rm = TRUE),
            n_other = sum(is.na(BRCA_Carrier)), statuses = paste(sort(unique(BRCA_Status)), collapse = ";"))
vars <- vars %>% left_join(var_status, by = "ID") %>% relocate(ID, Gene, HGVSp, n_calls, median_vaf, n_carrier, n_noncarrier, n_other, statuses)
fwrite(vars %>% arrange(ID), here("data/somatic_main_clinical_vars.csv"))
vars %>% filter(n_calls > 1) %>% select(Gene, HGVSp, n_calls, median_vaf, n_carrier, n_noncarrier, n_other) %>% arrange(desc(n_calls))
