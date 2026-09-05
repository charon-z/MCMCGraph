#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args) >= 1L) args[1L] else "results/runs"
output_dir <- if (length(args) >= 2L) args[2L] else "results"

csv_files <- list.files(run_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(csv_files)) stop("No benchmark result CSV files found in ", run_dir)

rows <- lapply(csv_files, utils::read.csv, stringsAsFactors = FALSE)
raw <- do.call(rbind, rows)

peak_memory_kb <- vapply(raw$run_id, function(run_id) {
  path <- file.path(run_dir, paste0(run_id, ".time.txt"))
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  match <- grep("Maximum resident set size", lines, value = TRUE)
  if (!length(match)) return(NA_real_)
  as.numeric(trimws(sub("^[^:]+:", "", match[1L])))
}, numeric(1))
raw$peak_ram_gb <- peak_memory_kb / 1024^2

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  raw, file.path(output_dir, "computational_performance_raw.csv"),
  row.names = FALSE
)

# `interaction()` drops rows containing NA. Comparator methods do not have an
# MCMC iteration count, so represent that value explicitly in the grouping key
# instead of silently omitting their results from the summary table.
niter_key <- ifelse(is.na(raw$niter), "not_applicable", as.character(raw$niter))
group_key <- interaction(
  raw$method, raw$dataset, raw$n, raw$dimensions, raw$J, niter_key,
  drop = TRUE, lex.order = TRUE
)
groups <- split(raw, group_key)
summary_rows <- lapply(groups, function(group) {
  data.frame(
    method = group$method[1L],
    dataset = group$dataset[1L],
    n = group$n[1L],
    dimensions = group$dimensions[1L],
    J = group$J[1L],
    niter = group$niter[1L],
    repeats = nrow(group),
    median_elapsed_seconds = stats::median(group$elapsed_seconds),
    min_elapsed_seconds = min(group$elapsed_seconds),
    max_elapsed_seconds = max(group$elapsed_seconds),
    median_peak_ram_gb = stats::median(group$peak_ram_gb, na.rm = TRUE),
    max_peak_ram_gb = max(group$peak_ram_gb, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[order(
  summary_table$dataset, summary_table$n, summary_table$J,
  summary_table$method
), , drop = FALSE]
row.names(summary_table) <- NULL

utils::write.csv(
  summary_table,
  file.path(output_dir, "computational_performance_summary.csv"),
  row.names = FALSE
)
print(summary_table)
