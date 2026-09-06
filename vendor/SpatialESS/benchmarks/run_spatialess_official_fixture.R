#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))
out_dir <- file.path(project_dir, "benchmarks", "results", "official_spatial_fixture")

suppressPackageStartupMessages({
  library(CellChat)
  library(SpatialESS)
})
fixture <- readRDS(file.path(out_dir, "fixture.rds"))
data(CellChatDB.human, package = "CellChat")

workflow_start <- proc.time()[["elapsed"]]
components <- prepare_cellchat_lr_components(
  fixture$lr, rownames(fixture$expression),
  CellChatDB.human$complex, CellChatDB.human$cofactor
)
graph <- build_radius_graph_csr(
  fixture$coordinates_um,
  radius = fixture$interaction_length_um,
  sample_id = rep("sample_1", ncol(fixture$expression)),
  weight = "binary", store_distance = FALSE
)
graph$cells$group_code <- as.integer(fixture$group)
graph$cells$group_levels <- levels(fixture$group)
support <- build_group_support_csr(graph, fixture$group)
prepared <- prepare_sparse_trimean(
  fixture$expression, genes = components$genes, normalize = TRUE
)
group_expression <- summarize_prepared_trimean(
  prepared, fixture$group, group_levels = support$group_levels
)
probability_start <- proc.time()[["elapsed"]]
observed <- score_cellchat_group_support(
  group_expression, support, components, Kh = 0.5, n = 1
)
probability_seconds <- proc.time()[["elapsed"]] - probability_start

global_blocks <- build_spatial_blocks(
  fixture$coordinates_um, block_size = 1000,
  sample_id = rep("sample_1", ncol(fixture$expression)),
  compartment = rep("tissue", ncol(fixture$expression)), origin = c(0, 0)
)
global_permutation <- permutation_cellchat_group_support(
  prepared, fixture$group, support, components, global_blocks,
  nperm = fixture$nboot, seed = fixture$seed,
  finite_correction = FALSE, tail = "strict_greater", verbose = FALSE
)
spatial_blocks <- build_spatial_blocks(
  fixture$coordinates_um, block_size = 75,
  sample_id = rep("sample_1", ncol(fixture$expression)),
  compartment = rep("tissue", ncol(fixture$expression)), origin = c(0, 0)
)
calibrated_permutation <- permutation_cellchat_group_support(
  prepared, fixture$group, support, components, spatial_blocks,
  nperm = 100, seed = fixture$seed,
  finite_correction = TRUE, tail = "greater_equal", verbose = FALSE
)
summary <- data.frame(
  method = "SpatialESS",
  cells = ncol(fixture$expression),
  groups = nlevels(fixture$group),
  lr = nrow(fixture$lr),
  graph_edges = length(graph$neighbors),
  supported_group_pairs = nrow(support$group_pairs),
  probability_seconds = probability_seconds,
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start,
  active_probabilities = nrow(observed$group_pairs),
  global_strict_significant_p_0_05 = sum(global_permutation$group_pairs$pvalue <= 0.05),
  spatial_calibrated_significant_p_0_05 = sum(calibrated_permutation$group_pairs$pvalue <= 0.05)
)
saveRDS(
  list(graph = graph, support = support, components = components,
       group_expression = group_expression, observed = observed,
       global_permutation = global_permutation,
       calibrated_permutation = calibrated_permutation),
  file.path(out_dir, "spatialess.rds"), compress = FALSE
)
write.table(summary, file.path(out_dir, "spatialess_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "spatialess_sessionInfo.txt"))
print(summary)

