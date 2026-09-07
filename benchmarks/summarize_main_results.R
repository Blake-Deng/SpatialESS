#!/usr/bin/env Rscript

# Copy small benchmark summaries into the repository. Large expression objects
# and full RDS outputs stay outside Git and are referenced by provenance paths.

script_file <- tryCatch(sys.frames()[[1L]]$ofile, error = function(e) NULL)
repo_root <- if (!is.null(script_file) && nzchar(script_file)) {
  normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = FALSE)
}
out_dir <- file.path(repo_root, "results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_first <- function(paths) {
  for (path in paths) {
    if (file.exists(path)) return(read.delim(path, check.names = FALSE))
  }
  NULL
}

copy_table <- function(name, paths) {
  value <- read_first(paths)
  if (!is.null(value)) {
    write.table(value, file.path(out_dir, name), sep = "\t",
                quote = FALSE, row.names = FALSE)
    message("wrote ", name)
  } else {
    message("source not found; retained packaged table: ", name)
  }
}

gse250346 <- Sys.getenv(
  "SPATIALESS_GSE250346_ROOT",
  "/data/dzf/GSE250346/benchmarks/multisample_v3_exact"
)
gse250346_lmm <- Sys.getenv(
  "SPATIALESS_GSE250346_LMM_ROOT",
  "/data/dzf/GSE250346/benchmarks/multisample/slice_glm_vs_patient_lmm"
)
gse313006 <- Sys.getenv(
  "SPATIALESS_GSE313006_ROOT",
  "/data/dzf/GSE313006/benchmarks/spatial_v3_headtohead"
)
cosmx <- Sys.getenv(
  "SPATIALESS_COSMX_ROOT",
  "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead"
)

copy_table("gse250346_exact_comparison.tsv", c(
  file.path(gse250346, "exact_comparison_summary.tsv")
))
copy_table("gse250346_lmm_model_summary.tsv", c(
  file.path(gse250346_lmm, "model_summary.tsv")
))
copy_table("gse250346_lmm_fit_status.tsv", c(
  file.path(gse250346_lmm, "lmm_status_summary.tsv")
))
copy_table("gse250346_lmm_failure_summary.tsv", c(
  file.path(gse250346_lmm, "lmm_failure_summary.tsv")
))
copy_table("gse250346_patient_robustness.tsv", c(
  file.path(gse250346, "patient_biology_robustness_20260904",
            "leave_one_patient_out_summary.tsv")
))
copy_table("gse250346_aggregation_sensitivity.tsv", c(
  file.path(gse250346, "patient_biology_robustness_20260904",
            "aggregation_sensitivity_summary.tsv")
))
copy_table("gse250346_min_cells_sensitivity.tsv", c(
  file.path(gse250346, "min_cells_sr_sensitivity_20260904",
            "threshold_patient_glm_comparison.tsv")
))
copy_table("cosmx_ladder_method_results.tsv", c(
  file.path(cosmx, "ladder_method_results.tsv")
))
copy_table("cosmx_ladder_accuracy_results.tsv", c(
  file.path(cosmx, "ladder_accuracy_results.tsv")
))
copy_table("gse313006_fidelity_20perm.tsv", c(
  file.path(gse313006, "cells001000_exact_perm00020", "comparison.tsv")
))
copy_table("gse313006_fidelity_100perm.tsv", c(
  file.path(gse313006, "cells001000_exact_perm00100", "comparison.tsv")
))

message("Summary refresh complete.")

