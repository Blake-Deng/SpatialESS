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
out_root <- file.path(root, "benchmarks", "results", "xenium_spatial_new_20260731")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

n_cells <- as.integer(Sys.getenv("XENIUM_SPATIAL_NCELLS", "5000"))
nperm <- as.integer(Sys.getenv("XENIUM_SPATIAL_NPERM", "20"))
panel <- Sys.getenv("XENIUM_SPATIAL_PANEL", "full")
seed <- as.integer(Sys.getenv("XENIUM_SPATIAL_SEED", "20260731"))
radius <- as.numeric(Sys.getenv("XENIUM_SPATIAL_RADIUS_UM", "30"))
block_size <- as.numeric(Sys.getenv("XENIUM_SPATIAL_BLOCK_UM", "100"))
if (!panel %in% c("small", "full")) stop("Invalid XENIUM_SPATIAL_PANEL")
if (n_cells < 20L || nperm < 1L || radius <= 0 || block_size <= 0) {
  stop("Invalid XENIUM_SPATIAL_NCELLS, NPERM, RADIUS_UM or BLOCK_UM")
}

object <- readRDS(prepared_path)
expression_all <- object@data.signaling
cell_id_all <- colnames(expression_all)
group_all <- object@idents
names(group_all) <- cell_id_all

coords_all <- read.csv(
  gzfile(coords_path),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
matched <- match(cell_id_all, coords_all$cell_id)
stopifnot(!anyNA(matched))
coords_all <- as.matrix(coords_all[matched, c("x_centroid", "y_centroid")])
rownames(coords_all) <- cell_id_all

center <- colMeans(coords_all)
selection_order <- order((coords_all[, 1] - center[1])^2 +
                         (coords_all[, 2] - center[2])^2)
if (n_cells >= ncol(expression_all)) {
  idx <- seq_len(ncol(expression_all))
  run_label <- "full"
} else {
  idx <- selection_order[seq_len(n_cells)]
  run_label <- sprintf("n%06d", n_cells)
}
run_label <- paste0(run_label, "_", panel)

expression <- expression_all[, idx, drop = FALSE]
coords <- coords_all[idx, , drop = FALSE]
group <- droplevels(group_all[idx])
names(group) <- colnames(expression)
lr <- object@LR$LRsig
if (panel == "small") {
  pathway <- if ("pathway_name" %in% colnames(lr)) as.character(lr$pathway_name) else rep("", nrow(lr))
  focus <- which(pathway %in% c("CXCL", "TNF", "IL17", "TGFb", "MIF"))
  keep <- c(focus, setdiff(seq_len(nrow(lr)), focus))
  lr <- lr[keep[seq_len(min(50L, length(keep)))], , drop = FALSE]
}
db <- object@DB

out_dir <- file.path(out_root, paste0(run_label, "_spatialess"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workflow_start <- proc.time()[["elapsed"]]

graph_start <- proc.time()[["elapsed"]]
graph <- build_radius_graph_csr(
  coords, radius = radius, sample_id = rep("xenium_ovary", ncol(expression)),
  weight = "gaussian", scale = radius / 2, store_distance = FALSE,
  max_edges = 5e8
)
graph_seconds <- proc.time()[["elapsed"]] - graph_start
stopifnot(isTRUE(attr(validate_spatial_csr_graph(graph), "valid")))

components <- prepare_cellchat_lr_components(
  lr, rownames(expression), db$complex, db$cofactor
)
prepared <- prepare_sparse_trimean(
  expression, genes = components$genes, normalize = TRUE
)
support_start <- proc.time()[["elapsed"]]
support <- build_group_support_csr(graph, group, max_group_pairs = 1e7)
support_seconds <- proc.time()[["elapsed"]] - support_start
blocks <- build_spatial_blocks(
  coords, block_size = block_size,
  sample_id = rep("xenium_ovary", ncol(expression)),
  compartment = rep("tissue", ncol(expression))
)

score_start <- proc.time()[["elapsed"]]
result <- permutation_cellchat_group_support(
  prepared, group, support, components, blocks,
  nperm = nperm, seed = seed,
  finite_correction = FALSE, tail = "strict_greater",
  max_score_entries = 1e8, verbose = FALSE
)
score_seconds <- proc.time()[["elapsed"]] - score_start
records <- result$group_pairs
if (nrow(records)) {
  index <- match(records$interaction_name, lr$interaction_name)
  records$pathway_name <- as.character(lr$pathway_name[index])
  records$annotation <- as.character(lr$annotation[index])
}

saveRDS(records, file.path(out_dir, "records.rds"), compress = FALSE)
write.table(records, gzfile(file.path(out_dir, "records.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(
  method = "SpatialESS",
  dataset = "Xenium Prime fresh-frozen human ovary",
  cells = ncol(expression), genes = nrow(expression),
  groups = nlevels(group), lr = nrow(lr), nperm = nperm,
  radius_um = radius, block_size_um = block_size,
  graph_edges = length(graph$neighbors),
  supported_group_pairs = nrow(support$group_pairs),
  active_group_lr_records = nrow(records),
  significant_p_0_05 = sum(records$pvalue <= 0.05),
  graph_seconds = graph_seconds,
  support_seconds = support_seconds,
  score_seconds = score_seconds,
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start
)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
