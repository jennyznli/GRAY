# Sequencing metrics overview: mean target coverage per sample, outliers flagged.
# Reads the two HsMetrics summaries (one per sequencing batch), not the master sheet.
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(data.table)
library(here)

theme_set(theme_bw(base_size = 12))
FIG_DIR <- here("figures/metrics")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(name, width = 8, height = 4.5) ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)

# ---- read
samples <- fread(here("input/gray_sample_ids_bams.tsv"), header = FALSE)$V1      # the 119 sequenced samples
# the xlsx stores its numbers as text, so read everything as text and convert
metrics <- bind_rows(
  `20240129` = fread(here("data/029-Gray_PreCancer.Target-HYB-PreCancer.20240129.metrics_summary.csv"), colClasses = "character") %>% as_tibble(),
  `20240923` = read_excel(here("data/029-Gray_PreCancer.Target-HYB-PreCancer.20240923.metrics_summary.xlsx"), col_types = "text"),
  .id = "batch") %>%
  filter(SAMPLE %in% samples) %>%
  mutate(across(-c(batch, SAMPLE, BAIT_SET), as.numeric))
nrow(metrics); sum(duplicated(metrics$SAMPLE)); setdiff(samples, metrics$SAMPLE)
# 119 0 character(0)

# ---- outliers: outside 1.5 x IQR of the cohort (the boxplot rule)
q <- quantile(metrics$MEAN_TARGET_COVERAGE, c(0.25, 0.5, 0.75))
LOW  <- q[[1]] - 1.5 * (q[[3]] - q[[1]])
HIGH <- q[[3]] + 1.5 * (q[[3]] - q[[1]])
metrics <- metrics %>%
  mutate(outlier = case_when(MEAN_TARGET_COVERAGE < LOW ~ "low", MEAN_TARGET_COVERAGE > HIGH ~ "high", TRUE ~ "no"),
         label = if_else(outlier == "no", "", SAMPLE))

summary(metrics$MEAN_TARGET_COVERAGE)
metrics %>% group_by(batch) %>%
  summarise(n = n(), median = median(MEAN_TARGET_COVERAGE), min = min(MEAN_TARGET_COVERAGE), max = max(MEAN_TARGET_COVERAGE))
c(LOW = LOW, HIGH = HIGH)
metrics %>% filter(outlier != "no") %>%
  select(SAMPLE, batch, outlier, MEAN_TARGET_COVERAGE, TOTAL_READS, PCT_DUP, PCT_SELECTED_BASES, PCT_TARGET_BASES_100X) %>%
  arrange(MEAN_TARGET_COVERAGE)

fwrite(metrics %>% select(-label) %>% arrange(MEAN_TARGET_COVERAGE), here("data/metrics_coverage.csv"))

# ---- 1. per sample, ranked
fig <- ggplot(metrics, aes(reorder(SAMPLE, MEAN_TARGET_COVERAGE), MEAN_TARGET_COVERAGE, color = batch)) +
  geom_hline(yintercept = q[[2]], linetype = "dashed", color = "grey50") +
  geom_point() +
  geom_text_repel(aes(label = label), size = 3, color = "black", min.segment.length = 0, show.legend = FALSE) +
  labs(title = "Mean target coverage per sample", subtitle = paste0("Line = cohort median ", round(q[[2]]), "x; labelled = outliers (1.5 x IQR)"),
       x = paste0("Samples (n = ", nrow(metrics), "), ranked"), y = "Mean target coverage", color = "Batch") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "top")
save_fig("01_mean_target_coverage_ranked", width = 8, height = 4.5)

# ---- 2. distribution
fig <- ggplot(metrics, aes(MEAN_TARGET_COVERAGE, fill = batch)) +
  geom_histogram(binwidth = 100, boundary = 0) +
  geom_vline(xintercept = q[[2]], linetype = "dashed") +
  labs(title = "Mean target coverage", subtitle = paste0("Line = cohort median ", round(q[[2]]), "x"), x = "Mean target coverage", y = "Samples", fill = "Batch") +
  theme(legend.position = "top")
save_fig("02_mean_target_coverage_hist", width = 6, height = 4.5)

# ---- 3. by batch
fig <- ggplot(metrics, aes(batch, MEAN_TARGET_COVERAGE)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = outlier != "no"), position = position_jitter(width = 0.15, seed = 1), alpha = 0.7) +
  geom_text_repel(aes(label = label), position = position_jitter(width = 0.15, seed = 1), size = 3, min.segment.length = 0) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "red"), guide = "none") +
  labs(title = "Mean target coverage by batch", x = "Batch", y = "Mean target coverage")
save_fig("03_mean_target_coverage_by_batch", width = 5, height = 4.5)
