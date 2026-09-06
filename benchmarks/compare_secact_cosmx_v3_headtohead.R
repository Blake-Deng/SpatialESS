#!/usr/bin/env Rscript

root <- Sys.getenv(
  "SECACT_V3_OUT",
  "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead"
)
cells <- as.integer(Sys.getenv("SECACT_V3_CELLS", "250"))
panel <- Sys.getenv("SECACT_V3_PANEL", "small")
nperm <- as.integer(Sys.getenv("SECACT_V3_NPERM", "20"))
dir_ess <- file.path(root, sprintf(
  "cells%06d_%s_spatialess_perm%05d", cells, panel, nperm
))
dir_v3 <- file.path(root, sprintf(
  "cells%06d_%s_spatialcellchat_v3_perm%05d", cells, panel, nperm
))
ess <- readRDS(file.path(dir_ess, "records.rds"))
v3 <- readRDS(file.path(dir_v3, "records.rds"))

key_columns <- c("interaction_name", "sender_group_name", "receiver_group_name")
key <- function(x) do.call(paste, c(x[key_columns], sep = "\r"))
jaccard <- function(a, b) {
  total <- union(a, b)
  if (!length(total)) return(1)
  length(intersect(a, b)) / length(total)
}
safe_cor <- function(x, y, method) {
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y, method = method)
}

ess$key <- key(ess)
v3$key <- key(v3)
joined <- merge(
  ess[, c(key_columns, "key", "probability", "pvalue", "pathway_name", "annotation")],
  v3[, c("key", "probability", "pvalue")],
  by = "key", suffixes = c("_ess", "_v3")
)
probability_difference <- abs(joined$probability_ess - joined$probability_v3)
pvalue_difference <- abs(joined$pvalue_ess - joined$pvalue_v3)
sig_ess <- ess$key[ess$pvalue <= 0.05]
sig_v3 <- v3$key[v3$pvalue <= 0.05]

summary <- data.frame(
  cells = cells,
  panel = panel,
  nperm = nperm,
  ess_active = nrow(ess),
  v3_active = nrow(v3),
  common_active = nrow(joined),
  active_jaccard = jaccard(ess$key, v3$key),
  max_abs_probability_difference = if (nrow(joined)) max(probability_difference) else NA_real_,
  mean_abs_probability_difference = if (nrow(joined)) mean(probability_difference) else NA_real_,
  probability_pearson = safe_cor(joined$probability_ess, joined$probability_v3, "pearson"),
  probability_spearman = safe_cor(joined$probability_ess, joined$probability_v3, "spearman"),
  max_abs_pvalue_difference = if (nrow(joined)) max(pvalue_difference) else NA_real_,
  exact_pvalue_fraction = if (nrow(joined)) mean(pvalue_difference == 0) else NA_real_,
  pvalue_spearman = safe_cor(joined$pvalue_ess, joined$pvalue_v3, "spearman"),
  ess_significant = length(sig_ess),
  v3_significant = length(sig_v3),
  significance_jaccard = jaccard(sig_ess, sig_v3),
  significance_threshold_disagreements = length(setdiff(union(sig_ess, sig_v3), intersect(sig_ess, sig_v3)))
)

aggregate_network <- function(x, suffix) {
  value <- aggregate(probability ~ sender_group_name + receiver_group_name, x, sum)
  names(value)[names(value) == "probability"] <- paste0("weight_", suffix)
  value
}
network <- merge(
  aggregate_network(ess, "ess"),
  aggregate_network(v3, "v3"),
  by = c("sender_group_name", "receiver_group_name"),
  all = TRUE
)

aggregate_pathway <- function(x, suffix) {
  value <- aggregate(probability ~ pathway_name, x, sum)
  value$rank <- rank(-value$probability, ties.method = "min")
  names(value)[names(value) == "probability"] <- paste0("weight_", suffix)
  names(value)[names(value) == "rank"] <- paste0("rank_", suffix)
  value
}
pathway <- merge(
  aggregate_pathway(ess, "ess"),
  aggregate_pathway(v3, "v3"),
  by = "pathway_name", all = TRUE
)

out_dir <- file.path(root, sprintf(
  "cells%06d_%s_comparison_perm%05d", cells, panel, nperm
))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(summary, file.path(out_dir, "comparison_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(joined, gzfile(file.path(out_dir, "shared_records.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(network, file.path(out_dir, "aggregate_network_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pathway, file.path(out_dir, "pathway_rank_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summary)

