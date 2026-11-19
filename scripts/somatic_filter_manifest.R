# config/somatic_filter.yaml callers: -> manifest for the filter reports (label<TAB>prefix, absolute prefixes),
# keeping only callers whose basic-stage qc file exists. Called by somatic_filter_report.bsub.
#   Rscript scripts/somatic_filter_manifest.R <config.yaml> <working dir> <out manifest.tsv>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: Rscript somatic_filter_manifest.R <config.yaml> <working dir> <out.tsv>")
cfg_path <- args[1]; wd <- args[2]; out <- args[3]
if (!requireNamespace("yaml", quietly = TRUE)) stop("needs the R yaml package -- install.packages('yaml')")

callers <- yaml::read_yaml(cfg_path)$callers
if (!length(callers)) stop(cfg_path, " has no callers: block")
abs_path <- function(p) if (startsWith(p, "/")) p else file.path(wd, p)

rows <- lapply(names(callers), function(key) {
    c <- callers[[key]]
    prefix <- if (!is.null(c$output_prefix)) c$output_prefix else sub("\\.csv$", "", c$input)
    prefix <- abs_path(prefix)
    qc <- paste0(prefix, ".basic_filtered.qc.tsv.gz")
    if (!file.exists(qc)) {
        message("WARNING: skipping ", key, " -- ", qc, " missing (run somatic_filter.bsub CALLER=", key, ")")
        return(NULL)
    }
    data.frame(label = if (!is.null(c$label)) c$label else key, prefix = prefix)
})
m <- do.call(rbind, rows)
if (is.null(m) || !nrow(m)) stop("no caller in ", cfg_path, " has basic-stage outputs yet")
if (anyDuplicated(m$label)) stop("duplicate caller labels in ", cfg_path, ": ", paste(m$label[duplicated(m$label)], collapse = ", "))
write.table(m, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(m), "\n", sep = "")    # caller count, read by the bsub
