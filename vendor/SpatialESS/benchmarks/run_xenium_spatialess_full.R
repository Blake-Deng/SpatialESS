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
nboot <- as.integer(Sys.getenv("XENIUM_V3_COMPAT_NBOOT", "0"))
out_dir <- file.path(
  root, "benchmarks", "results", "xenium_spatialess_20260801",
  sprintf("full_full_perm%05d", nboot)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

workflow_start <- proc.time()[["elapsed"]]
object <- readRDS(prepared_path)
expression <- object@data.signaling
cell_ids <- colnames(expression)
coords_table <- read.csv(
  gzfile(coords_path),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
matched <- match(cell_ids, coords_table$cell_id)
stopifnot(!anyNA(matched))
coordinates <- as.matrix(coords_table[matched, c("x_centroid", "y_centroid")])
rownames(coordinates) <- cell_ids
group <- droplevels(object@idents)

core_start <- proc.time()[["elapsed"]]
result <- spatialess(
  expression = expression,
  coordinates = coordinates,
  group = group,
  lr = object@LR$LRsig,
  complex = object@DB$complex,
  cofactor = object@DB$cofactor,
  ratio = 1, tol = 5,
  interaction.range = 30, scale.distance = 20,
  contact.range = 10,
  Kh = 0.5, n = 1, use.AGAN = TRUE,
  contact.dependent = TRUE,
  min.percent = 0, min.cells.sr = 1,
  nboot = nboot, seed.use = 20260731,
  max_edges = 1e8, max_output_records = 2e7,
  max_active_edges_per_lr = 1e8
)
core_seconds <- proc.time()[["elapsed"]] - core_start
records <- result$records
saveRDS(records, file.path(out_dir, "records.rds"), compress = FALSE)
write.table(records, gzfile(file.path(out_dir, "records.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(
  method = "SpatialESS",
  cells = ncol(expression), genes = nrow(expression),
  groups = nlevels(group), lr = nrow(result$lr), nboot = nboot,
  graph_edges = length(result$graph$neighbors),
  active_group_lr_records = nrow(records),
  significant_p_0_05 = if (nboot > 0L) sum(records$pvalue <= 0.05) else NA_integer_,
  core_seconds = core_seconds,
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start
)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
