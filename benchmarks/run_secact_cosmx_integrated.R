#!/usr/bin/env Rscript

project_dir <- normalizePath(file.path(dirname(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)), ".."))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(SpatialESS)
})

root <- Sys.getenv(
  "SECACT_V3_OUT",
  "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead"
)
out_root <- Sys.getenv(
  "SECACT_INTEGRATED_OUT",
  "/data/dzf/SecAct_2026/results/spatialess_integrated_stage2"
)
cells <- as.integer(Sys.getenv("SECACT_INTEGRATED_CELLS", "1000"))
nperm <- as.integer(Sys.getenv("SECACT_INTEGRATED_NPERM", "20"))
seed <- as.integer(Sys.getenv("SECACT_INTEGRATED_SEED", "20260804"))
inference <- Sys.getenv("SECACT_INTEGRATED_INFERENCE", "permutation")
contact_radius <- as.numeric(Sys.getenv("SECACT_INTEGRATED_CONTACT_UM", "15"))
diffusion_radius <- as.numeric(Sys.getenv("SECACT_INTEGRATED_DIFFUSION_UM", "35"))
block_size <- as.numeric(Sys.getenv("SECACT_INTEGRATED_BLOCK_UM", "100"))
if (is.na(cells) || cells < 20L || is.na(nperm) || nperm < 1L ||
    is.na(seed) || !inference %in% c("permutation", "guarded_hybrid") ||
    any(!is.finite(c(contact_radius, diffusion_radius, block_size))) ||
    any(c(contact_radius, diffusion_radius, block_size) <= 0)) {
  stop("Invalid integrated benchmark configuration.", call. = FALSE)
}

fixture <- readRDS(file.path(root, "secact_cosmx_v3_common_fixture.rds"))
if (cells > ncol(fixture$expression)) stop("Requested cells exceed fixture.")
idx <- fixture$selection_order[seq_len(cells)]
expression <- fixture$expression[, idx, drop = FALSE]
coordinates <- fixture$coordinates_um[idx, , drop = FALSE]
group <- droplevels(fixture$group[idx])
names(group) <- colnames(expression)
sample_id <- rep("CancerousLiver", cells)
compartment <- rep("all", cells)

run_dir <- file.path(out_root, sprintf(
  "cells%06d_full_%s_perm%05d", cells, inference, nperm
))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
start <- proc.time()[["elapsed"]]
result <- spatialess(
  expression = expression,
  coordinates = coordinates,
  group = group,
  lr = fixture$lr_full,
  complex = fixture$db$complex,
  cofactor = fixture$db$cofactor,
  sample_id = sample_id,
  compartment = compartment,
  contact_radius = contact_radius,
  diffusion_radius = diffusion_radius,
  block_size = block_size,
  inference = inference,
  nperm = nperm,
  seed = seed,
  Kh = 0.5,
  n = 1,
  alpha = 0.05,
  screening_margin = 0.01,
  finite_correction = TRUE,
  tail = "greater_equal",
  graph_weight = "binary",
  max_edges = 5e8,
  max_group_pairs = 1e7,
  max_score_entries = 1e8,
  max_records = 1e8,
  verbose = FALSE
)
elapsed <- proc.time()[["elapsed"]] - start

records <- result$records
records$cells <- cells
records$nperm <- nperm
saveRDS(records, file.path(run_dir, "records.rds"), compress = FALSE)
write.table(
  records, gzfile(file.path(run_dir, "records.tsv.gz")),
  sep = "	", quote = FALSE, row.names = FALSE
)
write.table(
  result$diagnostics, file.path(run_dir, "mechanism_diagnostics.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE
)
summary <- data.frame(
  method = "SpatialESS integrated",
  dataset = fixture$parameters$dataset,
  cells = cells,
  signaling_genes = nrow(expression),
  groups = nlevels(group),
  lr = nrow(fixture$lr_full),
  contact_lr = sum(fixture$lr_full$annotation == "Cell-Cell Contact"),
  diffusion_lr = sum(fixture$lr_full$annotation != "Cell-Cell Contact"),
  inference = inference,
  nperm = nperm,
  contact_radius_um = contact_radius,
  diffusion_radius_um = diffusion_radius,
  block_size_um = block_size,
  contact_edges = if ("contact" %in% names(result$graphs)) {
    length(result$graphs$contact$neighbors)
  } else 0,
  diffusion_edges = if ("diffusion" %in% names(result$graphs)) {
    length(result$graphs$diffusion$neighbors)
  } else 0,
  spatial_blocks = length(result$blocks$block_levels),
  active_records = nrow(records),
  significant_records_0_05 = sum(records$pvalue <= 0.05),
  unconfirmed_significant_records = sum(
    records$pvalue <= 0.05 & !records$permutation_confirmed
  ),
  method_end_to_end_seconds = elapsed,
  stringsAsFactors = FALSE
)
write.table(
  summary, file.path(run_dir, "summary.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
writeLines(c(
  sprintf("fixture=%s", normalizePath(file.path(
    root, "secact_cosmx_v3_common_fixture.rds"
  ))),
  sprintf("selection=%s", fixture$parameters$selection),
  "null=complete cell profiles permuted within 100 um spatial blocks",
  "pvalue=(greater_or_equal_null_count + 1) / (nperm + 1)",
  "all discoveries require fixed-permutation confirmation"
), file.path(run_dir, "provenance.txt"))
print(summary)
