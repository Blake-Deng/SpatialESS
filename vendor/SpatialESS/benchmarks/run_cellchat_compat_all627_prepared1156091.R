#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
graph_path <- file.path(
  "benchmarks", "results", "csr_prepared1156091_xenium30",
  "xenium_prepared1156091_radius30_csr.rds"
)
out_dir <- file.path(
  "benchmarks", "results", "cellchat_compat_all627_prepared1156091"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading the canonical prepared object and CSR graph...")
workflow_start <- proc.time()[["elapsed"]]
object <- readRDS(prepared_path)
graph <- readRDS(graph_path)
expression <- object@data.signaling
stopifnot(identical(colnames(expression), graph$cells$cell_id))
stopifnot(length(graph$cells$group_code) == ncol(expression))
lr <- object@LR$LRsig

components <- prepare_cellchat_lr_components(
  lr, rownames(expression), object@DB$complex, object@DB$cofactor
)
stopifnot(nrow(lr) == 627L, length(components$genes) == 534L)
saveRDS(components, file.path(out_dir, "lr_components_all627.rds"),
        compress = FALSE)

group_expression_path <- file.path(out_dir, "group_trimean_534genes_2357groups.rds")
if (file.exists(group_expression_path)) {
  message("Loading checkpointed group triMeans...")
  group_expression <- readRDS(group_expression_path)
  trimean_seconds <- NA_real_
} else {
  message("Computing exact sparse type-7 triMeans for 534 referenced genes...")
  stage_start <- proc.time()[["elapsed"]]
  group_expression <- summarize_cellchat_trimean(
    expression, graph$cells$group_code, components, normalize = TRUE
  )
  trimean_seconds <- proc.time()[["elapsed"]] - stage_start
  saveRDS(group_expression, group_expression_path, compress = FALSE)
}
colnames(group_expression$average) <- graph$cells$group_levels
group_expression$group_levels <- graph$cells$group_levels
saveRDS(group_expression, group_expression_path, compress = FALSE)


group_support_path <- file.path(out_dir, "group_support_2357groups.rds")
if (file.exists(group_support_path)) {
  message("Loading checkpointed group support...")
  group_support <- readRDS(group_support_path)
  support_seconds <- NA_real_
} else {
  message("Collapsing 24.4 million cell edges to sparse group support...")
  stage_start <- proc.time()[["elapsed"]]
  group_support <- build_group_support_csr(
    graph, graph$cells$group_code, max_group_pairs = 1e7
  )
  support_seconds <- proc.time()[["elapsed"]] - stage_start
  saveRDS(group_support, group_support_path, compress = FALSE)
}

message("Building an independent CellChat-helper reference for all 627 LR records...")
reference_start <- proc.time()[["elapsed"]]
average <- group_expression$average
reference_ligand <- CellChat::computeExpr_LR(
  as.character(lr$ligand), average, object@DB$complex
)
reference_receptor <- CellChat::computeExpr_LR(
  as.character(lr$receptor), average, object@DB$complex
)
reference_co_a <- CellChat::computeExpr_coreceptor(
  object@DB$cofactor, average, lr, type = "A"
)
reference_co_i <- CellChat::computeExpr_coreceptor(
  object@DB$cofactor, average, lr, type = "I"
)
reference_receptor <- reference_receptor * reference_co_a / reference_co_i
reference_agonist <- matrix(1, nrow = nrow(lr), ncol = ncol(average))
reference_antagonist <- matrix(1, nrow = nrow(lr), ncol = ncol(average))
agonist_index <- which(!is.na(lr$agonist) & nzchar(as.character(lr$agonist)))
antagonist_index <- which(!is.na(lr$antagonist) &
                          nzchar(as.character(lr$antagonist)))
for (i in agonist_index) {
  reference_agonist[i, ] <- CellChat::computeExpr_agonist(
    average, lr, object@DB$cofactor, i, Kh = 0.5, n = 1
  )
}
for (i in antagonist_index) {
  reference_antagonist[i, ] <- CellChat::computeExpr_antagonist(
    average, lr, object@DB$cofactor, i, Kh = 0.5, n = 1
  )
}
reference_seconds <- proc.time()[["elapsed"]] - reference_start

chunk_size <- 25L
chunk_indices <- split(seq_len(nrow(lr)),
                       ceiling(seq_len(nrow(lr)) / chunk_size))
manifest <- vector("list", length(chunk_indices))
maximum_probability_difference <- 0

for (chunk in seq_along(chunk_indices)) {
  index <- chunk_indices[[chunk]]
  chunk_path <- file.path(out_dir, sprintf("compat_lr_chunk_%02d.rds", chunk))
  if (file.exists(chunk_path)) {
    message(sprintf("Chunk %02d/%02d already exists; validating checkpoint.",
                    chunk, length(chunk_indices)))
    result <- readRDS(chunk_path)
    score_seconds <- result$benchmark$core_seconds
  } else {
    message(sprintf("Scoring chunk %02d/%02d: LR %d-%d",
                    chunk, length(chunk_indices), min(index), max(index)))
    chunk_components <- prepare_cellchat_lr_components(
      lr[index, , drop = FALSE], rownames(expression),
      object@DB$complex, object@DB$cofactor
    )
    stage_start <- proc.time()[["elapsed"]]
    result <- score_cellchat_group_support(
      group_expression, group_support, chunk_components,
      Kh = 0.5, n = 1, max_output_records = 2e7
    )
    score_seconds <- proc.time()[["elapsed"]] - stage_start
    result$group_pairs$global_lr_index <- index[result$group_pairs$lr_index]
    result$diagnostics$global_lr_index <- index[result$diagnostics$lr_index]
    result$benchmark <- list(core_seconds = score_seconds,
                             first_lr = min(index), last_lr = max(index))
    saveRDS(result, chunk_path, compress = FALSE)
  }

  pairs <- result$group_pairs
  if (nrow(pairs)) {
    global_lr <- pairs$global_lr_index
    sender <- pairs$sender_group
    receiver <- pairs$receiver_group
    product <- reference_ligand[cbind(global_lr, sender)] *
      reference_receptor[cbind(global_lr, receiver)]
    expected <- product / (0.5 + product) *
      reference_agonist[cbind(global_lr, sender)] *
      reference_agonist[cbind(global_lr, receiver)] *
      reference_antagonist[cbind(global_lr, sender)] *
      reference_antagonist[cbind(global_lr, receiver)]
    difference <- max(abs(pairs$molecular_probability - expected))
    maximum_probability_difference <- max(maximum_probability_difference,
                                          difference)
  } else {
    difference <- 0
  }
  manifest[[chunk]] <- data.frame(
    chunk = chunk, first_lr = min(index), last_lr = max(index),
    lr_count = length(index), emitted_group_pairs = nrow(pairs),
    core_seconds = score_seconds,
    max_cellchat_probability_difference = difference,
    result_mib = as.numeric(object.size(result)) / 1024^2,
    path = chunk_path
  )
  write.table(
    do.call(rbind, manifest[seq_len(chunk)]),
    file.path(out_dir, "chunk_manifest.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE
  )
  rm(result, pairs)
  gc()
}

manifest <- do.call(rbind, manifest)
summary <- data.frame(
  cells = ncol(expression), signaling_genes = nrow(expression),
  referenced_genes = length(components$genes),
  groups = length(group_support$group_levels),
  cell_graph_edges = length(graph$neighbors),
  supported_group_pairs = nrow(group_support$group_pairs),
  lr_records = nrow(lr),
  complex_lr_records = sum(lengths(components$ligand) > 1L |
                             lengths(components$receptor) > 1L),
  emitted_lr_group_pairs = sum(manifest$emitted_group_pairs),
  trimean_seconds = trimean_seconds,
  group_support_seconds = support_seconds,
  reference_helper_seconds = reference_seconds,
  score_core_seconds = sum(manifest$core_seconds),
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start,
  max_cellchat_probability_difference = maximum_probability_difference,
  maximum_chunk_result_mib = max(manifest$result_mib)
)
write.table(manifest, file.path(out_dir, "chunk_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
