#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))))
root <- Sys.getenv(
  "SPATIALESS_PSORIASIS_OUT",
  file.path(project_dir, "benchmarks", "results", "spatialcellchat_v3_psoriasis")
)
spots <- as.integer(Sys.getenv("SPATIAL_BENCH_SPOTS", "100"))
panel <- Sys.getenv("SPATIAL_BENCH_PANEL", "small")
nperm <- as.integer(Sys.getenv("SPATIAL_BENCH_NPERM", "20"))
dir_ess <- file.path(root, sprintf("spots%04d_%s_spatialess_perm%05d", spots, panel, nperm))
dir_v3 <- file.path(root, sprintf("spots%04d_%s_spatialcellchat_v3_perm%05d", spots, panel, nperm))
ess <- readRDS(file.path(dir_ess, "records.rds"))
v3 <- readRDS(file.path(dir_v3, "records.rds"))

key_columns <- c("interaction_name", "sender_group_name", "receiver_group_name")
key <- function(x) do.call(paste, c(x[key_columns], sep = "|"))
ess$key <- key(ess)
v3$key <- key(v3)
common <- intersect(ess$key, v3$key)
joined <- merge(
  ess[, c(key_columns, "key", "probability", "pvalue", "pathway_name", "annotation")],
  v3[, c("key", "probability", "pvalue")], by = "key", suffixes = c("_ess", "_v3")
)

safe_cor <- function(x, y, method) {
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y, method = method)
}
jaccard <- function(a, b) {
  u <- union(a, b)
  if (!length(u)) return(1)
  length(intersect(a, b)) / length(u)
}
summary <- data.frame(
  spots = spots, panel = panel, nperm = nperm,
  ess_active = nrow(ess), v3_active = nrow(v3), common_active = length(common),
  active_jaccard = jaccard(ess$key, v3$key),
  shared_probability_spearman = safe_cor(joined$probability_ess, joined$probability_v3, "spearman"),
  shared_pvalue_spearman = safe_cor(joined$pvalue_ess, joined$pvalue_v3, "spearman"),
  ess_significant = sum(ess$pvalue <= 0.05),
  v3_significant = sum(v3$pvalue <= 0.05),
  significant_jaccard = jaccard(ess$key[ess$pvalue <= 0.05], v3$key[v3$pvalue <= 0.05])
)

aggregate_pathway <- function(x, method) {
  z <- aggregate(probability ~ pathway_name, x, sum)
  z$rank <- rank(-z$probability, ties.method = "min")
  names(z)[names(z) == "probability"] <- paste0("weight_", method)
  names(z)[names(z) == "rank"] <- paste0("rank_", method)
  z
}
pathway <- merge(aggregate_pathway(ess, "ess"), aggregate_pathway(v3, "v3"),
                 by = "pathway_name", all = TRUE)
pathway$known_psoriasis_axis <- pathway$pathway_name %in% c("TNF", "CXCL", "IL17")
pathway <- pathway[order(pathway$rank_ess, pathway$rank_v3, na.last = TRUE), ]

out_dir <- file.path(root, sprintf("spots%04d_%s_comparison_perm%05d", spots, panel, nperm))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(summary, file.path(out_dir, "comparison_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(joined, file.path(out_dir, "shared_records.tsv.gz"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(pathway, file.path(out_dir, "pathway_rank_comparison.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(summary)
print(pathway[pathway$known_psoriasis_axis, , drop = FALSE])

