# Germline + functional filters on the calls kept by recurrence.R (local stand-in for germline.py / functional.py).
# Germline   = gnomAD max AF > GNOMAD_AF, or VAF > VAF_MAX; protected calls are exempt.
# Functional = on what is left: LoF level 3 / 4 removed, MANE Select transcripts only.
# Final sets = CH (CHIP whitelist), then somatic main / secondary for the rest.
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(here)

GNOMAD_AF <- 0.005      # gnomAD exome / genome max AF above this -> germline
VAF_MAX   <- 0.35       # VAF above this -> germline
VAF_TIER  <- c(0, 0.01, 0.05, 0.10, 0.35, Inf)
theme_set(theme_bw(base_size = 12))

# saves the plot just made to figures/germline/ (width / height in inches); plots are assigned to fig so RStudio does not draw them
FIG_DIR <- here("figures/germline")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
save_fig <- function(name, width = 8, height = 4.5) ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), width = width, height = height)

# ---- read
# gnomAD max = exome + genome, overall AF and populations except the bottlenecked / mixed ones
gnomad_cols <- c(paste0("gnomADe.", c("AF", "AFR", "AMR", "EAS", "MID", "NFE", "SAS")),
                 paste0("gnomADg.", c("AF", "AFR", "AMR", "EAS", "MID", "NFE", "SAS")))
calls <- fread(here("data/recurrence_calls.csv"), colClasses = "character", na.strings = NULL) %>%
  as_tibble() %>%
  mutate(across(c(Sample.Depth, Sample.AltDepth, Sample.AltFrac, Mutect2.AltFrac, VarDict.AltFrac, all_of(gnomad_cols)),
                ~ suppressWarnings(as.numeric(.x))),
         gnomad_max = coalesce(do.call(pmax, c(across(all_of(gnomad_cols)), na.rm = TRUE)), 0),    # not in gnomAD -> 0
         vaf_tier = cut(Sample.AltFrac, VAF_TIER, right = FALSE, labels = c("< 1%", "1-5%", "5-10%", "10-35%", ">= 35%")))
nrow(calls)
# 8260

# do the two VAFs even agree? (Sample.AltFrac = their mean); calls where the callers land on different sides of VAF_MAX
cor(calls$Mutect2.AltFrac, calls$VarDict.AltFrac, use = "complete.obs")
calls %>% count(mutect2_over = Mutect2.AltFrac > VAF_MAX, vardict_over = VarDict.AltFrac > VAF_MAX)
fig <- ggplot(calls, aes(Mutect2.AltFrac, VarDict.AltFrac)) +
  geom_abline(color = "grey60") + geom_point(alpha = 0.4, size = 0.8) +
  geom_vline(xintercept = VAF_MAX, linetype = "dashed") + geom_hline(yintercept = VAF_MAX, linetype = "dashed") +
  labs(title = "VAF per caller", x = "Mutect2 VAF", y = "VarDict VAF")
save_fig("00_vaf_mutect2_vs_vardict", width = 5, height = 5)

# ---- writers: call level + variant level (one row per variant, its first call's annotation + summaries over its calls)
to_vars <- function(d) d %>%
  group_by(ID) %>%
  mutate(n_calls = n(), median_vaf = median(Sample.AltFrac, na.rm = TRUE), median_alt_depth = median(Sample.AltDepth, na.rm = TRUE)) %>%
  slice(1) %>% ungroup() %>%
  select(-starts_with(c("Sample.", "Mutect2.", "VarDict.")), -callers, -starts_with("pass_")) %>%
  relocate(ID, any_of(c("germline", "functional", "set", "OncoKB.Oncogenic")), n_calls, gnomad_max, median_vaf, median_alt_depth)
write_pair <- function(d, name) {
  d <- d %>% select(-vaf_tier) %>% arrange(Sample.ID, ID)
  fwrite(d, here(paste0("data/", name, "_calls.csv")))
  fwrite(to_vars(d) %>% arrange(ID), here(paste0("data/", name, "_vars.csv")))
}

# ---- germline
calls <- calls %>%
  mutate(germline = case_when((gnomad_max > GNOMAD_AF | Sample.AltFrac > VAF_MAX) & protected == "TRUE" ~ "kept_protected",
                              gnomad_max > GNOMAD_AF ~ "removed_gnomad",
                              Sample.AltFrac > VAF_MAX ~ "removed_vaf",
                              TRUE ~ "kept"))
table(calls$germline)
# kept kept_protected removed_gnomad    removed_vaf 
# 1419            457           5093           1291 
# so 1419 + 457 = 1876 protected

calls %>% count(germline, vaf_tier) %>% pivot_wider(names_from = vaf_tier, values_from = n, values_fill = 0)
# germline       `< 1%` `1-5%` `5-10%` `10-35%` `>= 35%`
# <chr>           <int>  <int>   <int>    <int>    <int>
# 1 kept              847    467      62       43        0
# 2 kept_protected      0      0       0        0      457
# 3 removed_gnomad      9     31      31      170     4852
# 4 removed_vaf         0      0       0        0     1291

# protected calls are exempt (kept_protected = would have been removed)
calls %>% count(germline, protected) %>% pivot_wider(names_from = protected, values_from = n, values_fill = 0)

# ---- 1. VAF by outcome
fig <- ggplot(calls, aes(Sample.AltFrac, fill = germline)) +
  geom_histogram(binwidth = 0.02, boundary = 0) +
  geom_vline(xintercept = VAF_MAX, linetype = "dashed") +
  labs(title = "VAF by germline outcome", subtitle = paste0("Line = VAF ", VAF_MAX), x = "VAF", y = "Calls", fill = NULL) +
  theme(legend.position = "top")
save_fig("01_vaf_by_outcome", width = 6, height = 4.5)

# ---- 2. gnomAD max AF vs VAF
fig <- calls %>% filter(gnomad_max > 0) %>%
  ggplot(aes(gnomad_max, Sample.AltFrac, color = germline)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = GNOMAD_AF, linetype = "dashed") + geom_hline(yintercept = VAF_MAX, linetype = "dashed") +
  scale_x_log10() +
  labs(title = "gnomAD max AF vs VAF (calls found in gnomAD)", x = "gnomAD max AF (log)", y = "VAF", color = NULL) +
  theme(legend.position = "top")
save_fig("02_gnomad_vs_vaf", width = 6, height = 5)

# what are the rare gnomad but low VAF? 

# ---- write the germline outcome, then carry on with the kept calls only
write_pair(filter(calls, startsWith(germline, "removed")), "germline_removed")
calls <- filter(calls, startsWith(germline, "kept"))          # kept + kept_protected
write_pair(calls, "germline")
nrow(calls)
# 1876

# ---- functional: LoF level 3 / 4 removed, MANE Select transcripts only
calls <- calls %>%
  mutate(functional = case_when(MANE.Select %in% c("", ".") ~ "removed_no_mane",
                                Variant.LoF_level %in% c(3, 4) ~ "removed_lof",
                                TRUE ~ "kept"))
table(calls$functional)
# kept removed_lof 
# 508         814

calls %>% count(functional, Variant.LoF_level) %>%
  pivot_wider(names_from = Variant.LoF_level, values_from = n, values_fill = 0, names_prefix = "LoF_")
table(calls$Variant.LoF_level)
# 1   2   3   4 
# 444  64 814 554 

# ---- write the functional outcome, then keep the final calls
write_pair(filter(calls, startsWith(functional, "removed")), "functional_removed")
calls <- filter(calls, functional == "kept")
write_pair(calls, "functional")
nrow(calls)
calls %>% count(vaf_tier, protected) %>% pivot_wider(names_from = protected, values_from = n, values_fill = 0, names_prefix = "protected_")

# ---- 3. VAF of the final calls by LoF level
fig <- ggplot(calls, aes(Sample.AltFrac, fill = Variant.LoF_level)) +
  geom_histogram(bins = 40) + 
  # scale_x_log10() +
  labs(title = "VAF of the final calls", x = "VAF (log)", y = "Calls", fill = "LoF level") +
  theme(legend.position = "top")
save_fig("03_final_vaf_by_lof", width = 7, height = 4.5)


# ---- final sets
# ch                = on the CHIP whitelist (wl_whitelist; manual review calls do not count) -> not in the somatic sets
# somatic main      = AutoGVP or ClinVar pathogenic / likely pathogenic, or OncoKB in ONCOKB_MAIN
# somatic secondary = all else
ONCOKB_MAIN <- c("Oncogenic", "Likely Oncogenic")
PLP <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")

# OncoKB per variant + gene, from the per-caller oncokb csvs (first annotated value)
# bc we didn't merge right 
oncokb <- bind_rows(lapply(c("data/all_samples.mutect2.no_pon.vep.oncokb.csv", "data/all_samples.vardict.vep.oncokb.csv"), function(f)
  fread(here(f), select = c("ID", "Gene", "OncoKB.Oncogenic"), colClasses = "character"))) %>%
  as_tibble() %>%
  arrange(OncoKB.Oncogenic == ".") %>%
  distinct(ID, Gene, .keep_all = TRUE)

calls <- calls %>%
  left_join(oncokb, by = c("ID", "Gene")) %>%
  mutate(OncoKB.Oncogenic = coalesce(OncoKB.Oncogenic, "not_annotated"),
         set = case_when(wl_whitelist == "TRUE" ~ "ch",
                         AutoGVP %in% PLP | ClinVar.SIG %in% PLP | OncoKB.Oncogenic %in% ONCOKB_MAIN ~ "somatic_main",
                         TRUE ~ "somatic_secondary"))
table(calls$set)
# ch      somatic_main somatic_secondary 
# 144               293                71 
calls %>% count(set, OncoKB.Oncogenic) %>% pivot_wider(names_from = set, values_from = n, values_fill = 0)
  # OncoKB.Oncogenic    ch somatic_main somatic_secondary
  # <chr>            <int>        <int>             <int>
  # 1 Likely Oncogenic   144          277                 0
  # 2 Oncogenic            0           15                 0
  # 3 Unknown              0            1                71
  
calls %>% count(set, wl_manualreview)
# set               wl_manualreview     n
# <chr>             <chr>           <int>
# 1 ch                FALSE             144
# 2 somatic_main      FALSE             292
# 3 somatic_main      TRUE                1
# 4 somatic_secondary FALSE              71

write_pair(filter(calls, set == "somatic_main"), "somatic_main")
write_pair(filter(calls, set == "somatic_secondary"), "somatic_secondary")
write_pair(filter(calls, set == "ch"), "ch")

x <- unique(calls |> filter(set == "somatic_main"))
length(unique(x$Sample.ID))
# 47 ppl left 

