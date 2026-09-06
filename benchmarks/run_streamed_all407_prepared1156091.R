#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
graph_path <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/csr_prepared1156091_xenium30/xenium_prepared1156091_radius30_csr.rds"
out_dir <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/streamed_all407_prepared1156091"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading exact prepared object and CSR graph...")
object <- readRDS(prepared_path)
graph <- readRDS(graph_path)
expression <- object@data.signaling
stopifnot(identical(colnames(expression), graph$cells$cell_id))

lr <- object@LR$LRsig
simple <- lr$ligand %in% rownames(expression) & lr$receptor %in% rownames(expression)
lr_simple <- lr[simple, , drop = FALSE]
stopifnot(nrow(lr_simple) == 407L)

chunk_size <- 25L
chunk_id <- split(seq_len(nrow(lr_simple)),
                  ceiling(seq_len(nrow(lr_simple)) / chunk_size))
manifest <- vector("list", length(chunk_id))
all_diagnostics <- vector("list", length(chunk_id))

for (chunk in seq_along(chunk_id)) {
  index <- chunk_id[[chunk]]
  message(sprintf("Chunk %02d/%02d: LR %d-%d", chunk, length(chunk_id),
                  min(index), max(index)))
  start <- proc.time()[["elapsed"]]
  result <- aggregate_simple_lr_csr(
    expression, graph, lr_simple[index, , drop = FALSE],
    group = graph$cells$group_code,
    Kh = 0.5, n = 1, normalize = TRUE,
    min_expression = 0, max_group_pairs = 1e7
  )
  seconds <- proc.time()[["elapsed"]] - start
  result$diagnostics$global_lr_index <- index[result$diagnostics$lr_index]
  chunk_path <- file.path(out_dir, sprintf("streamed_simple_lr_chunk_%02d.rds", chunk))
  saveRDS(result, chunk_path, compress = FALSE)
  manifest[[chunk]] <- data.frame(
    chunk = chunk,
    first_lr = min(index),
    last_lr = max(index),
    lr_count = length(index),
    supported_edges = sum(result$diagnostics$supported_edges),
    sparse_group_pairs = nrow(result$group_pairs),
    core_seconds = seconds,
    result_mib = as.numeric(object.size(result)) / 1024^2,
    path = chunk_path
  )
  all_diagnostics[[chunk]] <- result$diagnostics
  write.table(do.call(rbind, manifest[seq_len(chunk)]),
              file.path(out_dir, "chunk_manifest.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  rm(result)
  gc()
}

manifest <- do.call(rbind, manifest)
diagnostics <- do.call(rbind, all_diagnostics)
diagnostics$visited_fraction <- diagnostics$visited_out_edges /
  diagnostics$total_graph_edges
diagnostics$supported_fraction_of_visited <- ifelse(
  diagnostics$visited_out_edges > 0,
  diagnostics$supported_edges / diagnostics$visited_out_edges, 0
)
summary <- data.frame(
  cells = ncol(expression),
  genes = nrow(expression),
  groups = max(graph$cells$group_code),
  graph_edges = length(graph$neighbors),
  total_lr_records = nrow(lr),
  completed_simple_lr = nrow(lr_simple),
  remaining_complex_or_missing_lr = sum(!simple),
  visited_out_edges = sum(diagnostics$visited_out_edges),
  supported_lr_edges = sum(diagnostics$supported_edges),
  sparse_lr_group_pairs = sum(manifest$sparse_group_pairs),
  total_core_seconds = sum(manifest$core_seconds),
  maximum_chunk_result_mib = max(manifest$result_mib)
)
write.table(manifest, file.path(out_dir, "chunk_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diagnostics, file.path(out_dir, "all407_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_dir, "all407_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
