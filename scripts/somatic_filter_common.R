# ========================
# Shared loader + plot helpers for the 3 somatic filter reports (sourced by each Rmd, next to this file):
#   basic_filter_report.Rmd       stage 1 (01-07)          needs PREFIX.basic_filtered.*
#   recurrent_filter_report.Rmd   stage 2 (08)             needs stages 1-2
#   annotation_filter_report.Rmd  stage 3 (09a-c) + caller overlap / intersection   needs stages 1-3
#   caller_regression_report.Rmd  per-sample concordance   needs stage 1
# The reports only VIEW the python outputs -- no threshold is re-applied here.
#
# Each stage writes PREFIX.<stage>.qc.tsv.gz keyed by row_id (input row number from stage 1). Stage 1's qc has
# every raw call; stages 2 / 3 only the rows that reached them. load_filter_data() stitches them into one
# row-per-raw-call table w/ last_step = index into the step list of the last step the row survived, so
# "rows that survived step k" is just last_step >= k (at_label()). Needs dplyr / tidyr / ggplot2 / data.table / knitr loaded.
# ========================

STAGES       <- c("basic", "recurrent", "annotation")
STAGE_SUFFIX <- c(basic = "basic_filtered", recurrent = "recurrent_filtered", annotation = "annotation_filtered")
# run_id this stage records for the stage before it (stale-chain check)
UPSTREAM_ID  <- c(recurrent = "basic_run_id", annotation = "recurrent_run_id")

stage_file <- function(prefix, stage, what = c("qc", "report")) {
    what <- match.arg(what)
    if (what == "qc") paste0(prefix, ".", STAGE_SUFFIX[[stage]], ".qc.tsv.gz")
    else paste0(prefix, ".", STAGE_SUFFIX[[stage]], "_report")
}

# "." placeholder and -1 sentinel -> NA
as_clean_num <- function(x) {
    x <- suppressWarnings(as.numeric(na_if(as.character(x), ".")))
    na_if(x, -1)
}

# all cols as character so per-file type guessing can't break bind_rows; gzip -dc (fread("x.gz") needs R.utils)
read_gz <- function(path) fread(cmd = paste("gzip -dc", shQuote(path)), colClasses = "character")

read_stage_tsv <- function(prefix, stage, file) {
    read.delim(file.path(stage_file(prefix, stage, "report"), file), colClasses = "character")
}

# VEP Variant.Class isn't in the qc files -- rebuilt from REF/ALT (VCF-style anchored alleles).
# Matched Variant.Class on every row of the merged mutect2 no_pon + vardict vep.csv (2026-11-11).
# "indel" = unanchored / complex (VarDict only so far)
var_class_of <- function(ref, alt) {
    case_when(nchar(ref) == 1 & nchar(alt) == 1          ~ "SNV",
              nchar(ref) == nchar(alt)                   ~ "substitution",
              nchar(ref) == 1 & startsWith(alt, ref)     ~ "insertion",
              nchar(alt) == 1 & startsWith(ref, alt)     ~ "deletion",
              TRUE                                       ~ "indel")
}
VAR_CLASSES <- c("SNV", "insertion", "deletion", "indel", "substitution")

# manifest (label<TAB>prefix, PREFIX = vep csv minus .csv, absolute) -> everything the reports use, for
# stages 1..upto. Callers missing any of those stages' outputs are skipped w/ a warning.
load_filter_data <- function(manifest_path, upto = 3) {
    stages   <- STAGES[seq_len(upto)]
    manifest <- read.delim(manifest_path, colClasses = "character")
    have     <- sapply(manifest$prefix, function(p) all(file.exists(sapply(stages, function(s) stage_file(p, s, "qc")))))
    if (any(!have)) warning("skipping (missing ", paste(stages, collapse = "/"), " outputs): ",
                            paste(manifest$label[!have], collapse = ", "))
    manifest <- manifest[have, , drop = FALSE]
    if (!nrow(manifest)) stop("no caller has all of the ", paste(stages, collapse = " / "), " stage outputs")
    callers  <- manifest$label

    # params per stage; a rerun stage w/o rerunning the ones after it -> stop (the stitched table would mix runs)
    run_params <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) bind_rows(lapply(stages, function(s) {
        read_stage_tsv(manifest$prefix[i], s, "params.tsv") %>% mutate(stage = s, method = manifest$label[i])
    }))))
    for (s in intersect(names(UPSTREAM_ID), stages)) {
        prev <- STAGES[match(s, STAGES) - 1]
        chk  <- run_params %>%
            filter((stage == s & param == UPSTREAM_ID[[s]]) | (stage == prev & param == "run_id")) %>%
            select(method, stage, value) %>%
            pivot_wider(names_from = stage, values_from = value)
        bad <- chk$method[chk[[s]] != chk[[prev]]]
        if (length(bad)) stop("the ", prev, " stage was rerun after the ", s, " stage for: ", paste(bad, collapse = ", "),
                              " -- rerun ", s, " (and every stage after it) so the reports don't mix runs")
    }
    run_params <- run_params %>% mutate(method = factor(method, levels = callers))

    funnel <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) bind_rows(lapply(stages, function(s) {
        read_stage_tsv(manifest$prefix[i], s, "funnel.tsv") %>% mutate(method = manifest$label[i])
    })))) %>%
        mutate(across(c(n_rows, pct_of_raw, n_unique_variants, n_PASS_rows, n_samples), as_clean_num))
    step_labels <- unique(funnel$step)
    stopifnot(all(table(funnel$step) == length(callers)))
    funnel <- funnel %>% mutate(step = factor(step, levels = step_labels), method = factor(method, levels = callers))

    qc <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
        p <- manifest$prefix[i]
        d <- read_gz(stage_file(p, "basic", "qc"))
        if (upto >= 2) {
            r <- read_gz(stage_file(p, "recurrent", "qc"))
            d <- r[d, on = "row_id"][, last_step := fifelse(is.na(last_step), i.last_step, last_step)][, i.last_step := NULL]
        } else {
            d[, `:=`(protected = NA_character_, n_carriers = NA_character_, flag_recurrent = NA_character_)]
        }
        if (upto >= 3) {
            a <- read_gz(stage_file(p, "annotation", "qc"))[, .(row_id, last_step)]
            d <- a[d, on = "row_id"][, last_step := fifelse(is.na(last_step), i.last_step, last_step)][, i.last_step := NULL]
        }
        as_tibble(d) %>% mutate(method = manifest$label[i])
    })) %>%
        mutate(method          = factor(method, levels = callers),
               last_step       = as.integer(last_step),
               n_carriers      = as_clean_num(n_carriers),
               Sample.AltFrac  = as_clean_num(Sample.AltFrac),
               Sample.Depth    = as_clean_num(Sample.Depth),
               Sample.AltDepth = as_clean_num(Sample.AltDepth),
               PASS            = FILTER == "PASS",
               # VEP's Variant.Class (parser column, in the qc file); rebuilt from REF/ALT only for older qc files
               var_class       = if ("Variant.Class" %in% names(.)) Variant.Class else var_class_of(REF, ALT),
               locus           = paste(Chr, Start, REF, ALT, sep = "_"),
               key             = paste(Sample.ID, locus, sep = "_"))

    list(manifest = manifest, CALLERS = callers, run_params = run_params, funnel = funnel,
         STEP_LABELS = step_labels, qc = qc,
         # caller overlap is only fair on samples every caller was run on
         shared_samples = Reduce(intersect, lapply(split(qc$Sample.ID, qc$method), unique)),
         caller_samples = qc %>% distinct(method, Sample.ID))
}

# ---- plot settings: `report:` block of config/somatic_filter.yaml (the python stages ignore it).
# Anything left out falls back to these defaults.
REPORT_DEFAULTS <- list(
    ad_cutoffs      = c(3, 5, 10, 30, 100),             # alt depth: density curves, sample-frequency splits, rows-left table
    vaf_cutoffs     = c(0.01, 0.03, 0.05, 0.07, 0.10),  # VAF: alt depth density by VAF cutoff
    recurrence_pcts = c(0.05, 0.10, 0.25, 0.50)         # fraction of samples carrying a variant: candidate remove_pct
                                                        # lines / table + rare-vs-recurrent splits
)
report_settings <- function(config_path) {
    rs <- REPORT_DEFAULTS
    if (nzchar(config_path) && file.exists(config_path)) {
        if (!requireNamespace("yaml", quietly = TRUE)) stop("the reports read config/somatic_filter.yaml -- install.packages('yaml')")
        cfg <- yaml::read_yaml(config_path)$report
        unknown <- setdiff(names(cfg), names(rs))
        if (length(unknown)) stop("unknown report: key(s) in ", config_path, ": ", paste(unknown, collapse = ", "),
                                  " -- allowed: ", paste(names(rs), collapse = ", "))
        for (k in names(cfg)) rs[[k]] <- unlist(cfg[[k]])
    }
    rs
}

# report tsv from one stage, every caller
read_report <- function(stage, file) {
    bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
        read_stage_tsv(manifest$prefix[i], stage, file) %>% mutate(method = manifest$label[i])
    })) %>%
        mutate(method = factor(method, levels = CALLERS))
}

# rows that survived the step w/ this exact label (last_step = 0-based index of the last step survived).
# By name, not number, so renumbering the steps can't silently shift which rows a report looks at.
at_label <- function(label) {
    i <- match(label, STEP_LABELS)
    if (is.na(i)) stop("step '", label, "' not in this run's funnel: ", paste(STEP_LABELS, collapse = ", "))
    qc %>% filter(last_step >= i - 1)
}
# the step sets the reports use
STEP_PRE   <- "03 MANE Select"         # preprocessed + on target + MANE = input to the basic thresholds
STEP_BASIC <- "07 VAF"                 # passed basic filters = input to recurrent_filter.py
STEP_REC08 <- "08 recurrent removed"   # passed 08 (before the 08b indel alt depth cutoff)
STEP_REC   <- "08b indel alt depth"    # passed recurrent stage (08 + 08b)
STEP_KEPT  <- "09c germline removed"   # kept

param_values <- function(p, stage = NULL) {
    keep <- run_params$param == p
    if (!is.null(stage)) keep <- keep & run_params$stage %in% stage
    run_params$value[keep]
}
param_value <- function(p, stage = NULL) as.numeric(param_values(p, stage)[1])
# for the Overview step lists: one value, or "a / b" if callers were run with different values
param_text  <- function(p, stage = NULL) paste(unique(param_values(p, stage)), collapse = " / ")
# values w/ commas get spaces so they wrap; fixed column widths via the .params-table CSS in the Rmd
params_table <- function(stage) {
    run_params %>%
        filter(stage == !!stage, !param %in% c("stage")) %>%
        mutate(value = gsub(",\\s*", ", ", value)) %>%
        select(method, param, value) %>%
        pivot_wider(names_from = method, values_from = value) %>%
        kable(format = "html", table.attr = 'class="table params-table"')
}

# one DT table per caller, each w/ an h4 title (no shared table to scroll through). Last value of a chunk.
per_caller_tables <- function(df, title, page_length = 10) {
    htmltools::tagList(lapply(CALLERS, function(m) {
        d <- df %>% filter(method == m) %>% select(-method)
        htmltools::tagList(htmltools::h4(paste0(title, " -- ", m)),
                           datatable(d, filter = "top", rownames = FALSE,
                                     options = list(pageLength = page_length, scrollX = TRUE)))
    }))
}

# ---- funnel
funnel_table <- function() {
    funnel %>%
        select(method, step, n_rows) %>%
        pivot_wider(names_from = method, values_from = n_rows) %>%
        kable(format.args = list(big.mark = ","))
}
funnel_plot <- function() {
    ggplot(funnel, aes(x = step, y = n_rows, group = method, color = method)) +
        geom_line() +
        geom_point() +
        scale_y_log10(labels = scales::comma) +
        labs(title = "Rows remaining at each step", x = NULL, y = "Rows (log scale)", color = NULL) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- core set, shown at every step; one facet column per caller
step_counts <- function(d) {
    d %>%
        group_by(method, .drop = FALSE) %>%
        summarise(rows             = n(),
                  unique_variants  = n_distinct(locus),
                  samples          = n_distinct(Sample.ID),
                  genes            = n_distinct(Gene),
                  median_VAF       = median(Sample.AltFrac, na.rm = TRUE),
                  median_depth     = median(Sample.Depth, na.rm = TRUE),
                  median_alt_depth = median(Sample.AltDepth, na.rm = TRUE),
                  .groups = "drop") %>%
        kable(format.args = list(big.mark = ","))
}

# funnel + step counts 
funnel_step_tables <- function(labels = STEP_LABELS, level = "###") {
    tab <- bind_rows(lapply(labels, function(lbl) {
        at_label(lbl) %>%
            group_by(method, .drop = FALSE) %>%
            summarise(rows             = n(),
                      n_unique_variants  = n_distinct(locus),
                      n_samples          = n_distinct(Sample.ID),
                      n_genes            = n_distinct(Gene),
                      median_VAF       = median(Sample.AltFrac, na.rm = TRUE),
                      median_depth     = median(Sample.Depth, na.rm = TRUE),
                      median_alt_depth = median(Sample.AltDepth, na.rm = TRUE),
                      .groups = "drop") %>%
            mutate(step = lbl)
    }))
    for (m in CALLERS) {
        cat("\n\n", level, " ", m, "\n\n", sep = "")
        print(tab %>%
                  filter(method == m) %>%
                  select(step, rows, n_unique_variants, n_samples, n_genes, median_VAF, median_depth, median_alt_depth) %>%
                  kable(format.args = list(big.mark = ",")))
    }
}

# one count table per step label (results='asis' chunk)
step_count_sections <- function(labels, level = "###") {
    for (lbl in labels) {
        cat("\n\n", level, " ", lbl, "\n\n", sep = "")
        print(step_counts(at_label(lbl)))
    }
}

plot_distributions <- function(d, label) {
    print(ggplot(d, aes(x = Sample.AltFrac)) +
        geom_histogram(bins = 50, fill = "darkorange", na.rm = TRUE) +
        facet_wrap(~ method, nrow = 1, scales = "free_y", drop = FALSE) +
        labs(title = paste("VAF -", label), x = "Sample.AltFrac (VAF)", y = "Count") +
        theme_bw())
    print(ggplot(d, aes(x = Sample.AltDepth)) +
        geom_histogram(bins = 50, fill = "darkgreen", na.rm = TRUE) +
        scale_x_log10() +
        facet_wrap(~ method, nrow = 1, scales = "free_y", drop = FALSE) +
        labs(title = paste("Alt allele depth -", label), x = "Sample.AltDepth (log scale)", y = "Count") +
        theme_bw())
    print(ggplot(d, aes(x = Sample.Depth)) +
        geom_histogram(bins = 50, fill = "purple", na.rm = TRUE) +
        scale_x_log10() +
        facet_wrap(~ method, nrow = 1, scales = "free_y", drop = FALSE) +
        labs(title = paste("Total depth -", label), x = "Sample.Depth (log scale)", y = "Count") +
        theme_bw())
    print(ggplot(d, aes(x = Variant.LoF_level)) +
        geom_bar(fill = "steelblue") +
        facet_wrap(~ method, nrow = 1, scales = "free_y", drop = FALSE) +
        labs(title = paste("LoF tier -", label), x = "Variant.LoF_level", y = "Count") +
        theme_bw())
}

# each facet has its own ranking -> key the axis by method+gene, lock order w/ factor(levels=)
plot_top_genes <- function(d, label, n_top = 25) {
    top <- d %>%
        count(method, Gene, name = "n") %>%
        group_by(method) %>%
        slice_max(n, n = n_top, with_ties = FALSE) %>%
        arrange(method, n) %>%
        ungroup() %>%
        mutate(gene_key = factor(paste(method, Gene, sep = "___"), levels = paste(method, Gene, sep = "___")))
    ggplot(top, aes(x = n, y = gene_key)) +
        geom_col(fill = "steelblue") +
        facet_wrap(~ method, nrow = 1, scales = "free") +
        scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
        labs(title = paste("Top", n_top, "genes by rows -", label), x = "Rows", y = NULL) +
        theme_bw()
}

# ---- per-sample counts by variant class (adapted from katya 01.2, but per caller, no caller averaging).
# sample universe per caller = every sample in that caller's raw input, so a sample w/ 0 calls at this
# step counts as 0 instead of dropping out of the denominator
class_per_sample <- function(d) {
    d %>%
        count(method, Sample.ID, var_class, name = "n") %>%
        right_join(caller_samples %>% crossing(var_class = VAR_CLASSES),
                   by = c("method", "Sample.ID", "var_class")) %>%
        mutate(n = coalesce(n, 0L),
               var_class = factor(var_class, levels = VAR_CLASSES))
}

class_counts_table <- function(d) {
    class_per_sample(d) %>%
        group_by(method, var_class) %>%
        summarise(median = median(n), .groups = "drop") %>%
        pivot_wider(names_from = var_class, values_from = median) %>%
        kable(caption = "Median calls per sample by variant class",
              format.args = list(big.mark = ","))
}

# one point per sample; free y per caller since VarDict is orders of magnitude above Mutect2
plot_class_per_sample <- function(d, label) {
    ggplot(class_per_sample(d), aes(x = var_class, y = n)) +
        geom_boxplot(outlier.shape = NA, fill = "grey90") +
        geom_jitter(width = 0.2, height = 0, size = 0.7, alpha = 0.5, colour = "steelblue") +
        facet_wrap(~ method, nrow = 1, scales = "free_y", drop = FALSE) +
        labs(title = paste("Calls per sample by variant class -", label),
             subtitle = "one point per sample; free y per caller",
             x = NULL, y = "Calls per sample") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# samples w/ the most calls per caller -- outliers to check (cf. MGHO-B1B2 in the MNV pilot)
top_samples_table <- function(d, n_top = 5) {
    class_per_sample(d) %>%
        group_by(method, Sample.ID) %>%
        summarise(total = sum(n), .groups = "drop_last") %>%
        mutate(x_median = round(total / median(total), 1)) %>%
        slice_max(total, n = n_top, with_ties = FALSE) %>%
        ungroup() %>%
        kable(caption = paste("Top", n_top, "samples by total calls per caller (x_median = fold over that caller's median sample)"),
              format.args = list(big.mark = ","))
}

# ---- caller overlap: same sample + locus + allele, shared samples only
plot_overlap <- function(d, label) {
    if (length(CALLERS) < 2) return(invisible(NULL))
    keys <- lapply(split(d$key[d$Sample.ID %in% shared_samples], d$method[d$Sample.ID %in% shared_samples]), unique)
    ggVennDiagram(keys) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        scale_x_continuous(expand = expansion(mult = 0.2)) +
        labs(title = paste("Caller overlap -", label),
             subtitle = paste0("Same sample and SNP (chr_start_REF_ALT); ", length(shared_samples), " shared samples"))
}

# ---- sample frequency of each variant among a set of calls (per caller), for the explorations
variant_freq_of <- function(d) {
    totals <- d %>% group_by(method) %>% summarise(total_samples = n_distinct(Sample.ID), .groups = "drop")
    d %>%
        distinct(method, Sample.ID, locus) %>%
        count(method, locus, name = "sample_count") %>%
        left_join(totals, by = "method") %>%
        mutate(sample_fraction = sample_count / total_samples)
}

# ---- FILTER tags (non-PASS rows): frequency per caller + one UpSet per caller (tags differ between callers)
plot_filter_freq <- function(d, label) {
    ff <- d %>%
        filter(!PASS) %>%
        select(method, FILTER) %>%
        separate_rows(FILTER, sep = ";") %>%
        count(method, FILTER, name = "n") %>%
        arrange(method, n) %>%
        mutate(tag_key = factor(paste(method, FILTER, sep = "___"), levels = paste(method, FILTER, sep = "___")))
    if (!nrow(ff)) return(invisible(NULL))
    ggplot(ff, aes(x = n, y = tag_key)) +
        geom_col(fill = "steelblue") +
        facet_wrap(~ method, nrow = 1, scales = "free") +
        scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
        labs(title = paste("FILTER reason frequency -", label, ""), x = "Rows flagged", y = NULL) +
        theme_bw()
}

# built from combination counts (FILTER strings are already tag;tag;...) instead of a row x tag matrix
plot_filter_upset <- function(d, label) {
    for (m in levels(d$method)) {
        combos <- d %>% filter(method == m, !PASS) %>% count(FILTER)
        if (nrow(combos) == 0) next
        tags <- unique(unlist(strsplit(combos$FILTER, ";")))
        if (length(tags) < 2) next
        print(upset(fromExpression(setNames(combos$n, gsub(";", "&", combos$FILTER, fixed = TRUE))),
                    nsets = length(tags), nintersects = 25, order.by = "freq",
                    mainbar.y.label = paste("Rows -", m, "-", label), sets.x.label = "Total flagged"))
        grid::grid.text(paste0("FILTER tag combinations - ", m, " (", label, ", non-PASS rows)"),
                        x = 0.5, y = 0.98, just = "top", gp = grid::gpar(fontsize = 13, fontface = "bold"))
    }
}
