#!/usr/bin/env Rscript

root <- "/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902"
out_dir <- file.path(root, "results", "official_rctd_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

record_key <- function(x) paste(
  x$sender_group_name, x$receiver_group_name, x$interaction_name, sep = "\r"
)

compare_records <- function(reference, query, method) {
  reference_key <- record_key(reference)
  query_key <- record_key(query)
  all_keys <- union(reference_key, query_key)
  ref_index <- match(all_keys, reference_key)
  query_index <- match(all_keys, query_key)
  ref_probability <- reference$probability[ref_index]
  query_probability <- query$probability[query_index]
  ref_pvalue <- reference$pvalue[ref_index]
  query_pvalue <- query$pvalue[query_index]
  ref_probability[is.na(ref_probability)] <- 0
  query_probability[is.na(query_probability)] <- 0
  ref_pvalue[is.na(ref_pvalue)] <- 1
  query_pvalue[is.na(query_pvalue)] <- 1
  probability_diff <- abs(ref_probability - query_probability)
  pvalue_diff <- abs(ref_pvalue - query_pvalue)
  shared <- intersect(reference_key, query_key)
  data.frame(
    method = method,
    official_active_records = nrow(reference),
    spatialess_active_records = nrow(query),
    shared_active_records = length(shared),
    active_jaccard = length(shared) / length(all_keys),
    max_abs_probability_difference = max(probability_diff),
    probability_exact_fraction = mean(probability_diff == 0),
    max_abs_pvalue_difference = max(pvalue_diff),
    pvalue_exact_fraction = mean(pvalue_diff == 0),
    significance_disagreements = sum(
      (ref_pvalue < 0.05) != (query_pvalue < 0.05)
    ),
    stringsAsFactors = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
method <- if (length(args) >= 1L) args[[1L]] else "sparkle"
run_tag <- if (length(args) >= 2L) args[[2L]] else ""
if (!method %in% c("raw", "sparkle")) {
  stop("Usage: compare_official_rctd.R raw|sparkle [run_tag]", call. = FALSE)
}
suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""

comparisons <- do.call(rbind, lapply(method, function(method) {
  official <- readRDS(file.path(
    root, "results", paste0("official_rctd_", method, suffix), "records.rds"
  ))
  spatialess <- readRDS(file.path(
    root, "results", paste0("spatialess_rctd_", method, suffix), "records.rds"
  ))
  compare_records(official, spatialess, method)
}))
write.table(
  comparisons, file.path(
    out_dir,
    paste0("official_vs_spatialess_", method, suffix, ".tsv")
  ),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(comparisons)
