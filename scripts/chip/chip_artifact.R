# ========================
# Artifact identification for CHIP variants (from 03-ch artifact.R + minad_age_association.R)
#   1. recurrent_flag   n QC-pass carriers >= min_carriers (chip_recurrence.py)
#   2. batch_flag       recurrent & BH-adjusted Fisher (carrier x batch) p < p_batch
#   3. age_flag         recurrent & not exempt & NOT (OR > 1 & p < p_age); needs --covariates
# review = recurrent_flag | wl_manualreview
# --batch-test / --age-test FALSE turn 2 / 3 off (flags -> FALSE)
# ========================
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

# ========================
# SET GLOBALS
# ========================
args <- commandArgs(trailingOnly = TRUE)
opt  <- setNames(as.list(args[c(FALSE, TRUE)]), sub("^--", "", args[c(TRUE, FALSE)]))
arg  <- function(k, default = NULL) if (!is.null(opt[[k]]) && opt[[k]] != "") opt[[k]] else default

BATCH_SIG <- as.numeric(arg("p-batch", "0.05"))
AGE_SIG   <- as.numeric(arg("p-age", "0.10"))
EXEMPT    <- strsplit(arg("exempt", ""), ",")[[1]]
BATCH_TEST <- as.logical(arg("batch-test", "TRUE"))
AGE_TEST   <- as.logical(arg("age-test", "TRUE"))

# ========================
# READ IN DATA
# ========================
rec   <- read_csv(arg("recurrence"), col_types = cols(variant_id = col_character(), Gene = col_character(),
                                                     recurrent_flag = col_logical(), wl_manualreview = col_logical()))
calls <- read_csv(arg("calls"), col_types = cols(.default = col_character()))
qc    <- read_csv(arg("sample-qc"), col_types = cols(.default = col_character())) %>%
    filter(sample_flags == "PASS")

carriers <- calls %>%
    filter(wl_whitelist == "TRUE" | wl_manualreview == "TRUE", Sample.ID %in% qc$sample_id) %>%
    distinct(variant_id = ID, Sample.ID)

# ========================
# 1. RECURRENT VARIANTS
# ========================
recurrent_var <- rec %>% filter(recurrent_flag) %>% pull(variant_id)
length(recurrent_var)

# ========================
# 2. BATCH ENRICHMENT
# ========================
batch_res <- tibble(variant_id = character(), p_batch = numeric())
if (!BATCH_TEST) {
    cat("batch test off\n")
} else if ("batch" %in% names(qc) && n_distinct(qc$batch) > 1) {
    batch_res <- tibble(variant_id = recurrent_var) %>%
        mutate(p_batch = sapply(variant_id, function(v) {
            tab <- table(qc$sample_id %in% carriers$Sample.ID[carriers$variant_id == v], qc$batch)
            tryCatch(fisher.test(tab, workspace = 2e8)$p.value,
                     error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value)
        }))
} else {
    cat("no batch column in sample_qc -- batch test skipped\n")
}
batch_res <- batch_res %>%
    mutate(p_batch_adj = p.adjust(p_batch, method = "BH"),
           batch_flag  = p_batch_adj < BATCH_SIG)

# ========================
# 3. AGE ASSOCIATION
# ========================
# carrier ~ every covariate column (sample_id dropped; batch as factor)
test_age_association <- function(v, cov) {
    df <- cov %>% mutate(carrier = as.integer(sample_id %in% carriers$Sample.ID[carriers$variant_id == v])) %>%
        select(-sample_id)
    tryCatch({
        m <- withCallingHandlers(
            glm(carrier ~ ., data = df, family = binomial()),
            warning = function(w) {
                if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
                    invokeRestart("muffleWarning")
            }
        )
        s  <- coef(summary(m))
        ci <- confint.default(m)
        tibble(variant_id = v, OR = exp(coef(m)[["age"]]), CI_lo = exp(ci["age", 1]), CI_hi = exp(ci["age", 2]),
               p_age = s["age", "Pr(>|z|)"], separation = any(fitted(m) %in% c(0, 1)) | !m$converged)
    }, error = function(e) {
        warning(sprintf("Model failed for %s: %s", v, conditionMessage(e)))
        tibble(variant_id = v, OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_, p_age = NA_real_, separation = NA)
    })
}

age_res <- tibble(variant_id = character(), OR = numeric(), CI_lo = numeric(), CI_hi = numeric(),
                  p_age = numeric(), separation = logical())
if (!AGE_TEST) {
    cat("age test off\n")
} else if (!is.null(arg("covariates"))) {
    cov <- read_tsv(arg("covariates"), col_types = cols(sample_id = col_character())) %>%
        filter(sample_id %in% qc$sample_id)
    stopifnot("age" %in% names(cov))
    if ("batch" %in% names(cov)) cov$batch <- as.factor(cov$batch)
    cat(sprintf("age test: %d QC-pass samples; covariates: %s\n", nrow(cov),
                paste(setdiff(names(cov), "sample_id"), collapse = ", ")))

    exempt_var <- rec %>% filter(sapply(strsplit(Gene, "&"), function(g) any(g %in% EXEMPT))) %>% pull(variant_id)
    age_res <- bind_rows(age_res, lapply(setdiff(recurrent_var, exempt_var), test_age_association, cov = cov))
}
age_res <- age_res %>% mutate(age_flag = !(!is.na(OR) & OR > 1 & p_age < AGE_SIG))

# ========================
# JOIN BACK - ALL FLAGS
# ========================
out <- rec %>%
    left_join(batch_res, by = "variant_id") %>%
    left_join(age_res, by = "variant_id") %>%
    mutate(across(c(batch_flag, age_flag), ~ ifelse(is.na(.x), FALSE, .x)),
           review = recurrent_flag | wl_manualreview) %>%
    arrange(desc(review), desc(n_carriers))
stopifnot(nrow(out) == nrow(rec), !any(grepl("\\.(x|y)$", names(out))))

calls_out <- calls %>%
    left_join(out %>% select(variant_id, n_carriers, recurrent_flag, p_batch_adj, batch_flag, OR, p_age, age_flag, review),
              by = c("ID" = "variant_id"))
stopifnot(nrow(calls_out) == nrow(calls))

flags <- c("recurrent_flag", "batch_flag", "age_flag", "review")
print(tibble(flag = flags,
             n_variants_flagged = sapply(flags, function(f) sum(out[[f]])),
             n_rows_flagged     = sapply(flags, function(f) sum(calls_out[[f]], na.rm = TRUE))))

write_csv(out, arg("out-variants"), na = ".")
write_csv(calls_out, arg("out-calls"), na = ".")
