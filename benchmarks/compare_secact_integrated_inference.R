#!/usr/bin/env Rscript

root <- Sys.getenv(
  "SECACT_INTEGRATED_OUT",
  "/data/dzf/SecAct_2026/results/spatialess_integrated_stage2"
)
cells <- as.integer(Sys.getenv("SECACT_INTEGRATED_CELLS", "1000"))
nperm <- as.integer(Sys.getenv("SECACT_INTEGRATED_NPERM", "20"))
alpha <- as.numeric(Sys.getenv("SECACT_INTEGRATED_ALPHA", "0.05"))
if (is.na(cells) || cells < 20L || is.na(nperm) || nperm < 1L ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
  stop("Invalid comparison configuration.", call. = FALSE)
}
run_dir <- function(mode) file.path(root, sprintf(
  "cells%06d_full_%s_perm%05d", cells, mode, nperm
))
permutation <- readRDS(file.path(run_dir("permutation"), "records.rds"))
hybrid <- readRDS(file.path(run_dir("guarded_hybrid"), "records.rds"))
key <- function(x) paste(
  x$mechanism, x$global_lr_index, x$sender_group, x$receiver_group,
  sep = ""
)
permutation_key <- key(permutation)
hybrid_key <- key(hybrid)
if (anyDuplicated(permutation_key) || anyDuplicated(hybrid_key)) {
  stop("Comparison keys are not unique.", call. = FALSE)
}
all_key <- union(permutation_key, hybrid_key)
permutation_index <- match(all_key, permutation_key)
hybrid_index <- match(all_key, hybrid_key)
permutation_p <- permutation$pvalue[permutation_index]
hybrid_p <- hybrid$pvalue[hybrid_index]
permutation_p[is.na(permutation_p)] <- 1
hybrid_p[is.na(hybrid_p)] <- 1
permutation_significant <- permutation_p <= alpha
hybrid_significant <- hybrid_p <= alpha
comparison <- data.frame(
  cells = cells,
  nperm = nperm,
  alpha = alpha,
  permutation_active = nrow(permutation),
  hybrid_active = nrow(hybrid),
  common_active = sum(!is.na(permutation_index) & !is.na(hybrid_index)),
  active_jaccard = sum(!is.na(permutation_index) & !is.na(hybrid_index)) /
    length(all_key),
  permutation_significant = sum(permutation_significant),
  hybrid_significant = sum(hybrid_significant),
  significance_jaccard = if (any(permutation_significant |
                                 hybrid_significant)) {
    sum(permutation_significant & hybrid_significant) /
      sum(permutation_significant | hybrid_significant)
  } else 1,
  significance_threshold_disagreements = sum(
    permutation_significant != hybrid_significant
  ),
  hybrid_unconfirmed_significant = sum(
    hybrid$pvalue <= alpha & !hybrid$permutation_confirmed
  ),
  screened_nonsignificant = sum(
    hybrid$inference_mode ==
      "experimental_hybrid_screened_nonsignificant"
  ),
  permutation_evaluated = sum(hybrid$permutation_confirmed),
  max_abs_pvalue_difference = max(abs(permutation_p - hybrid_p)),
  stringsAsFactors = FALSE
)
out_dir <- file.path(root, sprintf(
  "cells%06d_full_inference_comparison_perm%05d", cells, nperm
))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  comparison, file.path(out_dir, "summary.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE
)
disagreement <- data.frame(
  key = all_key,
  permutation_pvalue = permutation_p,
  hybrid_pvalue = hybrid_p,
  permutation_significant = permutation_significant,
  hybrid_significant = hybrid_significant,
  stringsAsFactors = FALSE
)
disagreement <- disagreement[
  disagreement$permutation_significant != disagreement$hybrid_significant,
  , drop = FALSE
]
write.table(
  disagreement, file.path(out_dir, "threshold_disagreements.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE
)
if (comparison$hybrid_unconfirmed_significant != 0L) {
  stop("Guarded hybrid contains an unconfirmed discovery.", call. = FALSE)
}
print(comparison)
