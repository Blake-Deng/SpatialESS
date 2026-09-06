#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_n5000_grid4.rds"
coords_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/Xenium_Prime_Human_Ovary_FF_cells.csv.gz"
out_dir <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/g1_xenium4995"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(prepared_path)
expression <- object@data.signaling
cell_id <- colnames(expression)

coord_all <- read.csv(
  gzfile(coords_path),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
matched <- match(cell_id, coord_all$cell_id)
stopifnot(!anyNA(matched))
coords <- as.matrix(coord_all[matched, c("x_centroid", "y_centroid")])
rownames(coords) <- cell_id
rm(coord_all)
gc()

graph_start <- proc.time()[["elapsed"]]
graph <- build_radius_graph(coords, radius = 30, weight = "gaussian", scale = 15)
graph_seconds <- proc.time()[["elapsed"]] - graph_start
stopifnot(isTRUE(attr(validate_spatial_graph(graph), "valid")))

lr <- object@LR$LRsig
simple <- lr$ligand %in% rownames(expression) & lr$receptor %in% rownames(expression)
lr_simple <- lr[simple, , drop = FALSE]
support <- Matrix::rowSums(expression > 0)
priority <- support[lr_simple$ligand] * support[lr_simple$receptor]
selected <- order(priority, decreasing = TRUE)[seq_len(min(5L, nrow(lr_simple)))]
lr_use <- lr_simple[selected, , drop = FALSE]

score_start <- proc.time()[["elapsed"]]
scores <- score_simple_lr_edges(expression, graph, lr_use, Kh = 0.5, n = 1,
                                normalize = TRUE, max_output_edges = 5e6)
score_seconds <- proc.time()[["elapsed"]] - score_start

summary <- data.frame(
  cells = ncol(expression), genes = nrow(expression),
  graph_edges = nrow(graph$edges), radius_um = 30,
  tested_lr = nrow(scores$lr), supported_lr_edges = nrow(scores$edges),
  graph_seconds = graph_seconds, score_seconds = score_seconds,
  graph_invariants_pass = isTRUE(attr(validate_spatial_graph(graph), "valid"))
)
write.table(summary, file.path(out_dir, "g1_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(scores$diagnostics, file.path(out_dir, "g1_lr_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(summary = summary, diagnostics = scores$diagnostics,
             edge_head = head(scores$edges, 1000), lr = scores$lr),
        file.path(out_dir, "g1_result_compact.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
print(scores$diagnostics)
