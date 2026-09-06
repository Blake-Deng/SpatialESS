#!/usr/bin/env Rscript

project_dir <- normalizePath(file.path(dirname(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)), ".."))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(SpatialESS)
})

prepared_path <- Sys.getenv(
  "XENIUM_INTEGRATED_PREPARED",
  "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
)
coordinates_path <- Sys.getenv(
  "XENIUM_INTEGRATED_COORDINATES",
  "/data/dzf/Xenium_Prime_Human_Ovary_FF/Xenium_Prime_Human_Ovary_FF_cells.csv.gz"
)
out_root <- Sys.getenv(
  "XENIUM_INTEGRATED_OUT",
  "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/spatialess_integrated_stage2"
)
nperm <- as.integer(Sys.getenv("XENIUM_INTEGRATED_NPERM", "20"))
seed <- as.integer(Sys.getenv("XENIUM_INTEGRATED_SEED", "20260804"))
contact_radius <- as.numeric(Sys.getenv("XENIUM_INTEGRATED_CONTACT_UM", "15"))
diffusion_radius <- as.numeric(Sys.getenv("XENIUM_INTEGRATED_DIFFUSION_UM", "35"))
block_size <- as.numeric(Sys.getenv("XENIUM_INTEGRATED_BLOCK_UM", "100"))
if (is.na(nperm) || nperm < 1L || is.na(seed) ||
    any(!is.finite(c(contact_radius, diffusion_radius, block_size))) ||
    any(c(contact_radius, diffusion_radius, block_size) <= 0)) {
  stop("Invalid Xenium integrated benchmark configuration.", call. = FALSE)
}

run_dir <- file.path(out_root, sprintf(
  "cells1156091_grid50_permutation_perm%05d", nperm
))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
workflow_start <- proc.time()[["elapsed"]]

object <- readRDS(prepared_path)
expression <- object@data.signaling
cell_ids <- colnames(expression)
coordinates_table <- read.csv(
  gzfile(coordinates_path),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
matched <- match(cell_ids, coordinates_table$cell_id)
if (anyNA(matched)) stop("Prepared cells are missing from coordinates.")
coordinates <- as.matrix(
  coordinates_table[matched, c("x_centroid", "y_centroid")]
)
rownames(coordinates) <- cell_ids
group <- droplevels(object@idents)
rm(coordinates_table, matched)
gc()

result <- spatialess(
  expression = expression,
  coordinates = coordinates,
  group = group,
  lr = object@LR$LRsig,
  complex = object@DB$complex,
  cofactor = object@DB$cofactor,
  sample_id = rep("Xenium_Ovary_FF", ncol(expression)),
  compartment = rep("all", ncol(expression)),
  contact_radius = contact_radius,
  diffusion_radius = diffusion_radius,
  block_size = block_size,
  inference = "permutation",
  nperm = nperm,
  seed = seed,
  Kh = 0.5,
  n = 1,
  alpha = 0.05,
  finite_correction = TRUE,
  tail = "greater_equal",
  graph_weight = "binary",
  max_edges = 2e8,
  max_group_pairs = 1e7,
  max_score_entries = 1e8,
  max_records = 2e7,
  verbose = FALSE
)

records <- result$records
saveRDS(records, file.path(run_dir, "records.rds"), compress = FALSE)
write.table(
  records, gzfile(file.path(run_dir, "records.tsv.gz")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  result$diagnostics, file.path(run_dir, "mechanism_diagnostics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
summary <- data.frame(
  method = "SpatialESS integrated",
  dataset = "10x Xenium Prime fresh-frozen human ovary",
  cells = ncol(expression),
  signaling_genes = nrow(expression),
  groups = nlevels(group),
  lr = nrow(object@LR$LRsig),
  contact_lr = sum(object@LR$LRsig$annotation == "Cell-Cell Contact"),
  diffusion_lr = sum(object@LR$LRsig$annotation != "Cell-Cell Contact"),
  inference = "permutation", nperm = nperm,
  contact_radius_um = contact_radius,
  diffusion_radius_um = diffusion_radius,
  block_size_um = block_size,
  contact_edges = length(result$graphs$contact$neighbors),
  diffusion_edges = length(result$graphs$diffusion$neighbors),
  spatial_blocks = length(result$blocks$block_levels),
  active_records = nrow(records),
  significant_records_0_05 = sum(records$pvalue <= 0.05),
  unconfirmed_significant_records = sum(
    records$pvalue <= 0.05 & !records$permutation_confirmed
  ),
  method_end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start,
  stringsAsFactors = FALSE
)
write.table(summary, file.path(run_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
writeLines(c(
  sprintf("prepared=%s", normalizePath(prepared_path)),
  sprintf("coordinates=%s", normalizePath(coordinates_path)),
  "null=complete cell profiles permuted within 100 um spatial blocks",
  "pvalue=(greater_or_equal_null_count + 1) / (nperm + 1)",
  "all discoveries require fixed-permutation confirmation"
), file.path(run_dir, "provenance.txt"))
print(summary)
