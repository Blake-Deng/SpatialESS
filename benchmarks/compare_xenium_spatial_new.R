#!/usr/bin/env Rscript

root <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/xenium_spatial_new_20260731"
n_cells <- Sys.getenv("XENIUM_SPATIAL_COMPARE_LABEL", "n005000")
ess <- readRDS(file.path(root, paste0(n_cells, "_spatialess"), "records.rds"))
v3 <- readRDS(file.path(root, paste0(n_cells, "_spatialcellchat_v3"), "records.rds"))

key <- function(x) paste(x$sender_group_name, x$receiver_group_name,
                         x$interaction_name, sep = "\r")
ke <- key(ess); kv <- key(v3)
common <- intersect(ke, kv)
ie <- match(common, ke); iv <- match(common, kv)
prob_spearman <- if (length(common) > 1L) {
  cor(ess$probability[ie], v3$probability[iv], method = "spearman")
} else NA_real_
p_spearman <- if (length(common) > 1L) {
  cor(ess$pvalue[ie], v3$pvalue[iv], method = "spearman")
} else NA_real_
sig_e <- ke[ess$pvalue <= 0.05]
sig_v <- kv[v3$pvalue <= 0.05]
sig_union <- union(sig_e, sig_v)
pathways <- c("CXCL", "TNF", "IL17")
pathway_summary <- do.call(rbind, lapply(pathways, function(pathway) {
  data.frame(
    pathway = pathway,
    ess_records = sum(ess$pathway_name == pathway, na.rm = TRUE),
    v3_records = sum(v3$pathway_name == pathway, na.rm = TRUE),
    ess_significant = sum(ess$pathway_name == pathway & ess$pvalue <= 0.05, na.rm = TRUE),
    v3_significant = sum(v3$pathway_name == pathway & v3$pvalue <= 0.05, na.rm = TRUE)
  )
}))
summary <- data.frame(
  subset = n_cells,
  ess_active = length(ke), v3_active = length(kv),
  common_active = length(common),
  active_jaccard = length(common) / length(union(ke, kv)),
  probability_spearman = prob_spearman,
  pvalue_spearman = p_spearman,
  significant_jaccard = length(intersect(sig_e, sig_v)) / length(sig_union)
)
out_dir <- file.path(root, paste0(n_cells, "_comparison"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(summary, file.path(out_dir, "comparison_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pathway_summary, file.path(out_dir, "pathway_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summary)
print(pathway_summary)
