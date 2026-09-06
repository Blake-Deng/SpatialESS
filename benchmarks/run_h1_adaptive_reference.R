#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(SpatialESS))

source(file.path(project_dir, "tests", "testthat", "helper-stratified-fixture.R"))
n_replicates <- as.integer(Sys.getenv("SPATIALESS_ADAPTIVE_REPLICATES", "10"))
nperm <- 999L
out_dir <- Sys.getenv(
  "SPATIALESS_ADAPTIVE_OUT",
  file.path(project_dir, "benchmarks", "results", "h1_adaptive_reference")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_one <- function(index) {
  fixture <- make_stratified_analytic_fixture(seed = 1000L + index)
  expression <- fixture$prepared
  fixed_start <- proc.time()[["elapsed"]]
  fixed <- permutation_cellchat_group_support(
    expression, fixture$group, fixture$support, fixture$components,
    fixture$blocks, nperm = nperm, seed = 2000L + index,
    finite_correction = TRUE, tail = "greater_equal", verbose = FALSE
  )
  fixed_seconds <- proc.time()[["elapsed"]] - fixed_start
  adaptive_start <- proc.time()[["elapsed"]]
  adaptive <- experimental_adaptive_permutation_cellchat_group_support(
    expression, fixture$group, fixture$support, fixture$components,
    fixture$blocks, min_permutations = 100L, batch_size = 100L,
    max_permutations = nperm, seed = 2000L,
    confidence_level = 0.95, tail = "greater_equal"
  )
  adaptive_seconds <- proc.time()[["elapsed"]] - adaptive_start
  fixed_pairs <- fixed$group_pairs
  adaptive_pairs <- adaptive$group_pairs
  fixed_key <- paste(fixed_pairs$lr_index, fixed_pairs$support_pair_index)
  adaptive_key <- paste(adaptive_pairs$lr_index, adaptive_pairs$support_pair_index)
  match_index <- match(adaptive_key, fixed_key)
  data.frame(
    replicate = index,
    fixed_seconds = fixed_seconds,
    adaptive_seconds = adaptive_seconds,
    adaptive_median_permutations = stats::median(
      adaptive_pairs$permutations_used
    ),
    adaptive_mean_permutations = mean(adaptive_pairs$permutations_used),
    adaptive_max_permutations = max(adaptive_pairs$permutations_used),
    adaptive_record_permutation_evaluations =
      adaptive$diagnostics$record_permutation_evaluations,
    adaptive_full_matrix_evaluations =
      adaptive$diagnostics$full_matrix_record_permutation_evaluations,
    adaptive_work_fraction =
      adaptive$diagnostics$record_permutation_evaluations /
      adaptive$diagnostics$full_matrix_record_permutation_evaluations,
    threshold_agreement = mean(
      (adaptive_pairs$pvalue <= 0.05) == (fixed_pairs$pvalue[match_index] <= 0.05)
    ),
    pvalue_mean_abs_difference = mean(abs(
      adaptive_pairs$pvalue - fixed_pairs$pvalue[match_index]
    )),
    stringsAsFactors = FALSE
  )
}

message(sprintf("Running %d fixed/adaptive comparisons...", n_replicates))
start <- proc.time()[["elapsed"]]
records <- do.call(rbind, lapply(seq_len(n_replicates), run_one))
wall_seconds <- proc.time()[["elapsed"]] - start
summary <- data.frame(
  replicates = n_replicates, fixed_nperm = nperm,
  adaptive_min = 100L, adaptive_batch = 100L, adaptive_max = nperm,
  median_fixed_seconds = stats::median(records$fixed_seconds),
  median_adaptive_seconds = stats::median(records$adaptive_seconds),
  median_adaptive_permutations = stats::median(
    records$adaptive_median_permutations
  ),
  mean_adaptive_permutations = mean(records$adaptive_mean_permutations),
  max_adaptive_permutations = max(records$adaptive_max_permutations),
  mean_adaptive_work_fraction = mean(records$adaptive_work_fraction),
  mean_threshold_agreement = mean(records$threshold_agreement),
  mean_abs_pvalue_difference = mean(records$pvalue_mean_abs_difference),
  wall_seconds = wall_seconds
)
write.table(records, file.path(out_dir, "records.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_dir, "summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
