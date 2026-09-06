#!/usr/bin/env Rscript

root <- "/home/dzf/cellchat_acceleration/SpatialESS"
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
coords_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/Xenium_Prime_Human_Ovary_FF_cells.csv.gz"
n_cells <- as.integer(Sys.getenv("XENIUM_V3_COMPAT_NCELLS", "250"))
panel <- Sys.getenv("XENIUM_V3_COMPAT_PANEL", "full")
if (!panel %in% c("small", "full")) stop("Invalid panel")
base_label <- paste0(sprintf("n%06d", n_cells), "_", panel)
reference_path <- file.path(
  root, "benchmarks", "results", "xenium_spatial_new_20260731",
  paste0(base_label, "_spatialcellchat_v3"), "records.rds"
)
out_dir <- file.path(
  root, "benchmarks", "results", "xenium_v3_compat_20260801", base_label
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(prepared_path)
expression_all <- object@data.signaling
cell_ids <- colnames(expression_all)
coords_table <- read.csv(
  gzfile(coords_path),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
matched <- match(cell_ids, coords_table$cell_id)
stopifnot(!anyNA(matched))
coords_all <- as.matrix(coords_table[matched, c("x_centroid", "y_centroid")])
rownames(coords_all) <- cell_ids
center <- colMeans(coords_all)
selection <- order((coords_all[, 1] - center[1])^2 +
                   (coords_all[, 2] - center[2])^2)[seq_len(n_cells)]
expression <- expression_all[, selection, drop = FALSE]
coordinates <- coords_all[selection, , drop = FALSE]
group <- droplevels(object@idents[selection])
lr <- object@LR$LRsig
if (panel == "small") {
  pathway <- as.character(lr$pathway_name)
  focus <- which(pathway %in% c("CXCL", "TNF", "IL17", "TGFb", "MIF"))
  keep <- c(focus, setdiff(seq_len(nrow(lr)), focus))
  lr <- lr[keep[seq_len(min(50L, length(keep)))], , drop = FALSE]
}

start <- proc.time()[["elapsed"]]
result <- spatial_v3_compatible(
  expression = expression,
  coordinates = coordinates,
  group = group,
  lr = lr,
  complex = object@DB$complex,
  cofactor = object@DB$cofactor,
  ratio = 1, tol = 5,
  interaction.range = 30, scale.distance = 20,
  contact.range = 10,
  Kh = 0.5, n = 1, use.AGAN = TRUE,
  contact.dependent = TRUE,
  min.percent = 0, min.cells.sr = 1,
  nboot = 20, seed.use = 20260731
)
elapsed <- proc.time()[["elapsed"]] - start
candidate <- result$records
reference <- readRDS(reference_path)

key <- function(x) paste(
  x$sender_group_name, x$receiver_group_name, x$interaction_name, sep = "\r"
)
candidate_key <- key(candidate)
reference_key <- key(reference)
common <- intersect(candidate_key, reference_key)
candidate_index <- match(common, candidate_key)
reference_index <- match(common, reference_key)
probability_difference <- abs(
  candidate$probability[candidate_index] - reference$probability[reference_index]
)
pvalue_difference <- abs(
  candidate$pvalue[candidate_index] - reference$pvalue[reference_index]
)
summary <- data.frame(
  candidate_active = length(candidate_key),
  reference_active = length(reference_key),
  common_active = length(common),
  active_jaccard = length(common) / length(union(candidate_key, reference_key)),
  max_probability_difference = if (length(common)) max(probability_difference) else NA_real_,
  probability_spearman = if (length(common) > 1L) cor(
    candidate$probability[candidate_index], reference$probability[reference_index],
    method = "spearman"
  ) else NA_real_,
  max_pvalue_difference = if (length(common)) max(pvalue_difference) else NA_real_,
  pvalue_exact_fraction = if (length(common)) mean(pvalue_difference == 0) else NA_real_,
  elapsed_seconds = elapsed
)
saveRDS(result, file.path(out_dir, "spatialess_v3_compatible.rds"), compress = FALSE)
write.table(summary, file.path(out_dir, "comparison_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
