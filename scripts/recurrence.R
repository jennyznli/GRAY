# Recurrence exploration on the both-caller calls (whitelist.py output) -> pick cutoffs for recurrence.py.
# Run section by section in RStudio. sample_freq = samples carrying the exact variant / all samples.
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(here)

# candidate sample frequency cutoffs
CUTS     <- c(0.03, 0.06, 0.10, 0.25, 0.50)      
VAF_TIER <- c(0, 0.01, 0.05, 0.10, 0.35, Inf)     
MIN_ADS  <- c(5, 7, 10, 20)
theme_set(theme_bw(base_size = 12))

FIG_DIR <- here("figures/recurrence")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(name, width = 8, height = 4.5) ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)

# ---- read
calls <- fread(here("data/whitelist.csv"), colClasses = "character", na.strings = NULL) %>%
  as_tibble() %>%
  mutate(across(c(Sample.Depth, Sample.AltDepth, Sample.AltFrac), ~ suppressWarnings(as.numeric(.x))),
         type = if_else(Variant.Class %in% c("insertion", "deletion", "indel"), "indel", "SNV / MNV"),
         vaf_tier = cut(Sample.AltFrac, VAF_TIER, right = FALSE,
                        labels = c("< 1%", "1-5%", "5-10%", "10-35%", ">= 35%")))
dim(calls)
N <- n_distinct(calls$Sample.ID)

table(calls$vaf_tier)
# < 1%   1-5%  5-10% 10-35% >= 35% 
# 9093  25246   5960   5476  63249 
length(unique(calls$ID))
# 5835 unique variants

# every merged call (either caller), for the per-caller carrier counts
all_calls <- fread(here("data/merged.csv"), select = c("Sample.ID", "ID", "callers"), colClasses = "character") %>%
  as_tibble()
dim(all_calls)
# 919323

# ---- carriers per variant: both callers, and each caller on its own
bins <- c("1 sample", paste0("<= ", 100 * CUTS[1], "%"),
          paste0(100 * head(CUTS, -1), "-", 100 * CUTS[-1], "%"), paste0("> ", 100 * tail(CUTS, 1), "%"), "100%")

# count how many were in the original merged, not both 
per_caller <- function(name) {
  d <- filter(all_calls, grepl(name, callers))
  d %>% distinct(ID, Sample.ID) %>% count(ID, name = "n") %>% mutate(freq = n / n_distinct(d$Sample.ID))
}

# annotate original recurrence
# sample_freq is defined by n_carriers / N in the INTERSECTION cohort 
variants <- calls %>%
  group_by(ID, type) %>%
  summarise(n_carriers = n_distinct(Sample.ID), median_vaf = median(Sample.AltFrac, na.rm = TRUE),
            protected = any(protected == "TRUE"), .groups = "drop") %>%
  left_join(per_caller("Mutect2") %>% rename(n_carriers_mutect2 = n, freq_mutect2 = freq), by = "ID") %>%
  left_join(per_caller("VarDict") %>% rename(n_carriers_vardict = n, freq_vardict = freq), by = "ID") %>%
  mutate(sample_freq = n_carriers / N,
         # highest sample frequency in either caller's own calls
         max_freq = pmax(freq_mutect2, freq_vardict),     
         # + 1e-9: 0.06 * N must land in "<= 6%", not the next bin (float error)
         rec_bin = factor(case_when(n_carriers == 1 ~ bins[1],
                                    n_carriers == N ~ "100%",
                                    TRUE ~ bins[-1][findInterval(sample_freq - 1e-9, CUTS) + 1]), levels = bins))
calls <- calls %>% left_join(select(variants, ID, n_carriers, sample_freq, max_freq, rec_bin,
                                    n_carriers_mutect2, n_carriers_vardict), by = "ID")

dim(calls)
# 109024

# ---- variants in every sample: removed even when protected (protected ones written out to check)
variants %>% filter(n_carriers == N) %>% count(protected)
# protected     n
# <lgl>     <int>
# 1 FALSE       125
# 2 TRUE          3

calls %>% filter(n_carriers == N) %>%
  fwrite(here("data/recurrence_remove_Nsample_calls.csv"))
variants %>% filter(n_carriers == N) %>%
  # annotation = the variant's first call row, minus the per-sample / per-caller columns
  left_join(calls %>% distinct(ID, .keep_all = TRUE) %>%
              select(-any_of(setdiff(names(variants), "ID")),
                     -starts_with(c("Sample.", "Mutect2.", "VarDict.")),
                     -callers, -starts_with("pass_"), -vaf_tier), by = "ID") %>%
  fwrite(here("data/recurrence_remove_Nsamples_vars.csv"))

calls %>% filter(n_carriers == N, protected == "TRUE") %>%
  fwrite(here("data/recurrence_remove_Nsample_protected_calls.csv"))
variants %>% filter(n_carriers == N, protected == "TRUE") %>%
  # annotation = the variant's first call row, minus the per-sample / per-caller columns
  left_join(calls %>% distinct(ID, .keep_all = TRUE) %>%
              select(-any_of(setdiff(names(variants), "ID")), -starts_with(c("Sample.", "Mutect2.", "VarDict.")),
                     -callers, -starts_with("pass_"), -vaf_tier), by = "ID") %>%
  fwrite(here("data/recurrence_remove_Nsample_protected_vars.csv"))

calls_x <- filter(calls, n_carriers < N)             # everything below excludes the 100% variants
variants_x <- filter(variants, n_carriers < N)
dim(calls_x)
# 93792
dim(variants_x)
# 5707 

# ---- 1. sample frequency
fig <- ggplot(variants_x, aes(sample_freq)) +
  geom_histogram(binwidth = 0.01, boundary = 0, fill = "steelblue") +
  geom_vline(xintercept = CUTS, linetype = "dashed") +
  # scale_y_log10() +
  labs(title = "Unique variants by fraction of samples carrying them",
       x = "Sample frequency", y = "Unique variants")
save_fig("01_sample_freq", width = 7, height = 4.5)

# ---- 2. what each cutoff removes (protected rows are exempt, shown for scale)
cut_table <- function(d) bind_rows(lapply(CUTS, function(p) {
  hit <- d$sample_freq > p + 1e-9
  tibble(cutoff = p, min_carriers = floor(p * N + 1e-9) + 1, variants = n_distinct(d$ID[hit]),
         rows = sum(hit), pct_rows = round(100 * mean(hit), 1), protected_rows = sum(hit & d$protected == "TRUE"),
         median_vaf = median(d$Sample.AltFrac[hit], na.rm = TRUE))
}))
cut_table(calls_x)
# cutoff min_carriers variants  rows pct_rows protected_rows median_vaf
# <dbl>        <dbl>    <int> <int>    <dbl>          <int>      <dbl>
# 1   0.03            4     2131 89002     94.9            904      0.476
# 2   0.06            8     1690 86699     92.4            870      0.476
# 3   0.1            12     1448 84445     90              870      0.476
# 4   0.25           30     1003 75874     80.9            796      0.476
# 5   0.5            60      658 61086     65.1            587      0.477

calls_x %>% group_by(vaf_tier) %>% group_modify(~ cut_table(.x)) %>% print(n = Inf)
# vaf_tier cutoff min_carriers variants  rows pct_rows protected_rows median_vaf
# <fct>     <dbl>        <dbl>    <int> <int>    <dbl>          <int>      <dbl>
#   1 < 1%       0.03            4      481  7995     90.3            360     0.0073
# 2 < 1%       0.06            8      395  7635     86.3            346     0.0074
# 3 < 1%       0.1            12      340  7231     81.7            346     0.0074
# 4 < 1%       0.25           30      217  5566     62.9            294     0.0077
# 5 < 1%       0.5            60      135  3537     40              145     0.0079
# 6 1-5%       0.03            4      528 19439     98.6            140     0.0184
# 7 1-5%       0.06            8      470 19265     97.8            132     0.0184
# 8 1-5%       0.1            12      427 19078     96.8            132     0.0185
# 9 1-5%       0.25           30      330 18224     92.5            127     0.0187
# 10 1-5%       0.5            60      240 16036     81.4            102     0.0192
# 11 5-10%      0.03            4      190  4547     98.5              1     0.0702
# 12 5-10%      0.06            8      179  4515     97.8              0     0.0703
# 13 5-10%      0.1            12      169  4485     97.2              0     0.0703
# 14 5-10%      0.25           30      146  4312     93.4              0     0.0701
# 15 5-10%      0.5            60      117  3872     83.9              0     0.0701
# 16 10-35%     0.03            4      264  3896     96.9              0     0.159 
# 17 10-35%     0.06            8      244  3813     94.9              0     0.157 
# 18 10-35%     0.1            12      218  3707     92.2              0     0.155 
# 19 10-35%     0.25           30      167  3208     79.8              0     0.148 
# 20 10-35%     0.5            60      106  2454     61                0     0.136 
# 21 >= 35%     0.03            4     1456 53125     93.9            403     0.499 
# 22 >= 35%     0.06            8     1135 51471     90.9            392     0.500 
# 23 >= 35%     0.1            12      963 49944     88.2            392     0.500 
# 24 >= 35%     0.25           30      664 44564     78.7            375     0.503 
# 25 >= 35%     0.5            60      423 35187     62.2            340     0.510 

# ---- 3. VAF by sample frequency bin
# all calls, so the variants in every sample (100%) get their own curve
fig <- ggplot(calls, aes(Sample.AltFrac, color = rec_bin)) +
  geom_density(adjust = 1.5, na.rm = TRUE) +
  geom_vline(xintercept = c(0.5, 1), linetype = "dotted", color = "grey40") +
  # scale_x_log10() + 
  scale_color_manual(values = setNames(c(scales::viridis_pal(end = 0.9)(length(bins) - 1), "firebrick"), bins)) +
  labs(title = "VAF by sample frequency bin", x = "VAF", y = "Density", color = "sample_freq") +
  theme(legend.position = "top")
save_fig("03_vaf_by_freq_bin", width = 5, height = 4.5)

calls_x %>% count(rec_bin, vaf_tier) %>% pivot_wider(names_from = vaf_tier, values_from = n, values_fill = 0)
# rec_bin  `< 1%` `1-5%` `5-10%` `10-35%` `>= 35%`
# <fct>     <int>  <int>   <int>    <int>    <int>
# 1 1 sample    520    164      36       58     1858
# 2 <= 3%       335    105      33       66     1615
# 3 3-6%        360    174      32       83     1654
# 4 6-10%       404    187      30      106     1527
# 5 10-25%     1665    854     173      499     5380
# 6 25-50%     2029   2188     440      754     9377
# 7 > 50%      3537  16036    3872     2454    35187

# ---- 4. sample frequency vs median VAF, one point per variant
fig <- ggplot(variants_x, aes(sample_freq, median_vaf, color = type)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_vline(xintercept = CUTS, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = head(VAF_TIER[-1], -1), linetype = "dotted", color = "grey40") +
  # scale_y_log10() +
  labs(title = "Sample frequency vs median VAF per variant",
       x = "Sample frequency", y = "Median VAF", color = NULL) +
  theme(legend.position = "top")
save_fig("04_freq_vs_median_vaf", width = 6, height = 5)

# ---- 5. recurrence within VAF tiers (the top tier is mostly germline)
fig <- ggplot(calls_x, aes(sample_freq)) +
  geom_histogram(binwidth = 0.02, boundary = 0, fill = "steelblue") +
  geom_vline(xintercept = CUTS, linetype = "dashed") +
  facet_wrap(~ vaf_tier, nrow = 1, scales = "free_y") +
  labs(title = "Calls by sample frequency of their variant, per VAF tier", x = "Sample frequency", y = "Calls")
save_fig("05_freq_by_vaf_tier", width = 15, height = 4)

# ---- 6. recurrence when only calls with alt depth >= cutoff count as carriers
by_ad <- bind_rows(lapply(MIN_ADS, function(t) {
  calls_x %>% filter(Sample.AltDepth >= t) %>%
    group_by(ID) %>% summarise(sample_freq = n_distinct(Sample.ID) / N, .groups = "drop") %>%
    mutate(min_ad = paste("AltDepth >=", t))
})) %>% mutate(min_ad = factor(min_ad, levels = paste("AltDepth >=", MIN_ADS)))

fig <- ggplot(by_ad, aes(sample_freq)) +
  geom_histogram(binwidth = 0.02, boundary = 0, fill = "steelblue") +
  geom_vline(xintercept = CUTS, linetype = "dashed") +
  scale_y_log10() + facet_wrap(~ min_ad, nrow = 1) +
  labs(title = "Unique variants by sample frequency, by alt depth cutoff", x = "Sample frequency", y = "Unique variants (log)")

save_fig("06_freq_by_min_ad", width = 12, height = 4)

by_ad %>% group_by(min_ad) %>%
  summarise(variants = n(), !!!setNames(lapply(CUTS, function(p) expr(sum(sample_freq > !!p + 1e-9))), paste0("over_", 100 * CUTS, "pct")))
# min_ad         variants over_3pct over_6pct over_10pct over_25pct over_50pct
# <fct>             <int>     <int>     <int>      <int>      <int>      <int>
# 1 AltDepth >= 5      5696      2130      1690       1447       1003        658
# 2 AltDepth >= 7      5618      2121      1684       1441       1002        658
# 3 AltDepth >= 10     5429      2085      1659       1414        988        650
# 4 AltDepth >= 20     4890      1882      1494       1279        890        575


# ---- 7. carriers per caller (each caller's own calls) vs carriers called by both
fig <- ggplot(variants_x, aes(freq_mutect2, freq_vardict, color = median_vaf)) +
  geom_abline(color = "grey60") + geom_point(alpha = 0.6, size = 0.8) +
  # geom_vline(xintercept = CUTS, linetype = "dashed", color = "grey70") +
  # geom_hline(yintercept = CUTS, linetype = "dashed", color = "grey70") +
  scale_color_viridis_c(trans = "log10", breaks = c(0.001, 0.01, 0.1, 0.5), direction = -1, end = 0.95) +
  labs(title = "Sample frequency of each variant: Mutect2 vs VarDict", x = "Mutect2 sample frequency", y = "VarDict sample frequency", color = "median VAF") +
  theme(legend.position = "top", legend.key.width = unit(1, "cm"))
save_fig("07_carriers_mutect2_vs_vardict", width = 5, height = 5)
# maybe i should filter out the recurrence by the max in one caller?? 

# variants that look rare in the intersection but are recurrent in one caller
sus <- variants_x %>%
  mutate(extra_vardict = n_carriers_vardict - n_carriers, extra_mutect2 = n_carriers_mutect2 - n_carriers) %>%
  filter(sample_freq <= CUTS[2], pmax(freq_mutect2, freq_vardict) > CUTS[2]) %>%
  arrange(desc(pmax(extra_vardict, extra_mutect2)))

# ---- 8. VAF density by alt depth and sample frequency thresholds
# per panel: three dashed references (all calls, AltDepth < / >= split), then rare (<= cutoff) vs recurrent (> cutoff)
# calls, each split by alt depth. Each curve is its own density (area 1), so heights don't show how many calls it has.
GROUPS <- c("all calls", "AltDepth < split", "AltDepth >= split", "rare, AltDepth < split", "rare, AltDepth >= split",
            "recurrent, AltDepth < split", "recurrent, AltDepth >= split")
vaf_groups <- bind_rows(lapply(MIN_ADS, function(t) bind_rows(lapply(CUTS, function(p) {
  g <- calls_x %>% filter(!is.na(Sample.AltFrac), !is.na(Sample.AltDepth)) %>%
    transmute(Sample.AltFrac, low = if_else(Sample.AltDepth < t, "<", ">="),
              rare = if_else(max_freq <= p + 1e-9, "rare", "recurrent"))
  bind_rows(mutate(g, group = "all calls"),
            mutate(g, group = paste("AltDepth", low, "split")),
            mutate(g, group = paste0(rare, ", AltDepth ", low, " split"))) %>%
    mutate(ad_split = paste("AltDepth split", t), freq_cutoff = paste0("Sample frequency cutoff ", 100 * p, "%"))
})))) %>%
  mutate(group = factor(group, levels = GROUPS),
         ad_split = factor(ad_split, levels = paste("AltDepth split", MIN_ADS)),
         freq_cutoff = factor(freq_cutoff, levels = paste0("Sample frequency cutoff ", 100 * CUTS, "%")))
vaf_groups <- vaf_groups %>% add_count(ad_split, freq_cutoff, group) %>% filter(n >= 20)      # too few calls for a density
plot_vaf_groups <- function(d, title) {
  ggplot(d, aes(Sample.AltFrac, color = group, linetype = group)) +
    geom_density(adjust = 2) +
    geom_vline(xintercept = c(0.5, 1), linetype = "dotted", color = "grey40") +
    scale_x_log10() +
    scale_linetype_manual(values = setNames(c("dashed", "dashed", "dashed", rep("solid", 4)), GROUPS)) +
    scale_color_manual(values = setNames(c("grey30", "grey60", "#8c6bb1", "#1b9e77", "#e79e29", "#2585a8", "#e7298a"), GROUPS)) +
    labs(title = title, subtitle = "rare = max per-caller sample frequency <= cutoff, recurrent = above; dotted = VAF 0.5 / 1.0",
         x = "VAF (log)", y = "Density", color = NULL, linetype = NULL) +
    theme(legend.position = "top") + guides(color = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3))
}

# all cutoffs in one grid
# plot_vaf_groups(vaf_groups, "VAF density by alt depth and sample frequency thresholds") +
#   facet_grid(freq_cutoff ~ ad_split, scales = "free_y")
# save_fig("08_vaf_by_ad_and_freq", width = 13, height = 12)
# 
# one plot + PDF per alt depth split x sample frequency cutoff -> figures/recurrence/08_vaf_by_ad_and_freq/
dir.create(file.path(FIG_DIR, "08_vaf_by_ad_and_freq"), showWarnings = FALSE)
for (t in MIN_ADS) for (p in CUTS) {
  d <- vaf_groups %>% filter(ad_split == paste("AltDepth split", t), freq_cutoff == paste0("Sample frequency cutoff ", 100 * p, "%"))
  fig <- plot_vaf_groups(d, paste0("VAF density: AltDepth split ", t, ", sample frequency cutoff ", 100 * p, "%"))
  save_fig(sprintf("08_vaf_by_ad_and_freq/ad%02d_freq%02d", t, round(100 * p)), width = 6, height = 5)
}

# calls per curve
vaf_groups %>% count(ad_split, freq_cutoff, group) %>% pivot_wider(names_from = group, values_from = n, values_fill = 0) %>% print(n = Inf)

# ---- 9. recurrence filter (local stand-in for recurrence.py)
# on max_freq (>= each cutoff):
#   REMOVE_ALL  removed even when protected
#   REMOVE      removed, unless protected -> kept, labelled kept_protected
#   FLAG        kept, labelled kept_flagged
REMOVE_ALL <- 0.40
REMOVE     <- 0.10
FLAG       <- 0.06
calls <- calls %>%
  mutate(recurrence = case_when(max_freq >= REMOVE_ALL - 1e-9 ~ "removed_high",
                                max_freq >= REMOVE - 1e-9 & protected != "TRUE" ~ "removed",
                                max_freq >= REMOVE - 1e-9 ~ "kept_protected",
                                max_freq >= FLAG - 1e-9 ~ "kept_flagged",
                                TRUE ~ "kept"))

table(calls$recurrence)
# kept        kept_flagged      kept_protected             removed removed_all_samples 
# 3538                3817                 905               80793               19971 

calls %>% count(recurrence, vaf_tier) %>% pivot_wider(names_from = vaf_tier, values_from = n, values_fill = 0)
# recurrence          `< 1%` `1-5%` `5-10%` `10-35%` `>= 35%`
# <chr>                <int>  <int>   <int>    <int>    <int>
# 1 kept                   190    191      58       86     3013
# 2 kept_flagged           287    173      35      127     3195
# 3 kept_protected         379    134       0        0      392
# 4 removed               7531  17115    3652     3314    49181
# 5 removed_all_samples    706   7633    2215     1949     7468

out <- calls %>% select(-type, -vaf_tier, -rec_bin)
# variant level: one row per variant (its first call's annotation) + summaries over its calls
to_vars <- function(d) d %>%
  group_by(ID) %>%
  mutate(n_calls = n(), median_vaf = median(Sample.AltFrac, na.rm = TRUE), median_alt_depth = median(Sample.AltDepth, na.rm = TRUE)) %>%
  slice(1) %>% ungroup() %>%
  select(-starts_with(c("Sample.", "Mutect2.", "VarDict.")), -callers, -starts_with("pass_")) %>%
  relocate(ID, recurrence, n_calls, n_carriers, n_carriers_mutect2, n_carriers_vardict, sample_freq, max_freq, median_vaf, median_alt_depth)
kept <- filter(out, startsWith(recurrence, "kept"))
dim(kept)
# 8260
removed <- filter(out, startsWith(recurrence, "removed"))

fwrite(kept, here("data/recurrence_calls.csv"))
fwrite(to_vars(kept), here("data/recurrence_vars.csv"))
fwrite(removed, here("data/recurrence_removed_calls.csv"))
fwrite(to_vars(removed), here("data/recurrence_removed_vars.csv"))

dim(kept)
# kept
dim(removed)
# 100764 

dim(to_vars(kept))
# 3622
dim(to_vars(removed))
# 2213



