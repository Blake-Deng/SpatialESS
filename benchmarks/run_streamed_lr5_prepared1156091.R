#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
graph_path <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/csr_prepared1156091_xenium30/xenium_prepared1156091_radius30_csr.rds"
out_dir <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/streamed_lr5_prepared1156091"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading exact prepared object and CSR graph...")
object <- readRDS(prepared_path)
graph <- readRDS(graph_path)
expression <- object@data.signaling
stopifnot(identical(colnames(expression), graph$cells$cell_id))
stopifnot(length(graph$cells$group_code) == ncol(expression))

lr <- object@LR$LRsig
simple <- lr$ligand %in% rownames(expression) & lr$receptor %in% rownames(expression)
lr_simple <- lr[simple, , drop = FALSE]
support <- Matrix::rowSums(expression > 0)
priority <- as.double(support[lr_simple$ligand]) * as.double(support[lr_simple$receptor])
selected <- order(priority, decreasing = TRUE)[seq_len(min(5L, nrow(lr_simple)))]
lr_use <- lr_simple[selected, , drop = FALSE]

message("Streaming five LR pairs over the full graph...")
start <- proc.time()[["elapsed"]]
result <- aggregate_simple_lr_csr(
  expression, graph, lr_use, group = graph$cells$group_code,
  Kh = 0.5, n = 1, normalize = TRUE,
  min_expression = 0, max_group_pairs = 1e7
)
score_seconds <- proc.time()[["elapsed"]] - start

diagnostics <- result$diagnostics
diagnostics$visited_fraction <- diagnostics$visited_out_edges / diagnostics$total_graph_edges
diagnostics$supported_fraction_of_visited <-
  ifelse(diagnostics$visited_out_edges > 0,
         diagnostics$supported_edges / diagnostics$visited_out_edges, 0)
summary <- data.frame(
  cells = ncol(expression),
  genes = nrow(expression),
  groups = max(graph$cells$group_code),
  graph_edges = length(graph$neighbors),
  tested_lr = nrow(result$lr),
  visited_out_edges = sum(diagnostics$visited_out_edges),
  supported_lr_edges = sum(diagnostics$supported_edges),
  sparse_lr_group_pairs = nrow(result$group_pairs),
  score_seconds = score_seconds,
  aggregate_object_mib = as.numeric(object.size(result)) / 1024^2
)

write.table(summary, file.path(out_dir, "streamed_lr5_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diagnostics, file.path(out_dir, "streamed_lr5_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(result, file.path(out_dir, "streamed_lr5_result.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
print(diagnostics)
