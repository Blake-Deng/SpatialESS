#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(SpatialESSV3Exact)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1L]] %in% c("raw", "sparkle")) {
  stop("Usage: run_spatialess_rctd.R raw|sparkle", call. = FALSE)
}
method <- args[[1L]]
root <- "/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902"
input_dir <- file.path(root, "inputs", "spatialess_rctd")
nboot <- as.integer(Sys.getenv("RNA_LEAKAGE_NBOOT", "100"))
seed <- as.integer(Sys.getenv("RNA_LEAKAGE_SEED", "20260903"))
cell_limit <- as.integer(Sys.getenv("RNA_LEAKAGE_CELL_LIMIT", "0"))
lr_limit <- as.integer(Sys.getenv("RNA_LEAKAGE_LR_LIMIT", "0"))
run_tag <- Sys.getenv("RNA_LEAKAGE_RUN_TAG", "")
if (is.na(nboot) || nboot < 1L || is.na(seed) || is.na(cell_limit) ||
    cell_limit < 0L || is.na(lr_limit) || lr_limit < 0L) {
  stop("Invalid benchmark environment variables.", call. = FALSE)
}
suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""
out_dir <- file.path(root, "results", paste0("spatialess_rctd_", method, suffix))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(
  input_dir,
  if (method == "raw") "raw_log1p_cp10k.mtx.gz" else
    "sparkle_log1p_cp10k.mtx.gz"
)
genes <- readLines(file.path(input_dir, "genes.tsv"))
cells <- readLines(file.path(input_dir, "cells.tsv"))
metadata <- read.delim(
  file.path(input_dir, "cell_metadata.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
expression <- methods::as(readMM(gzfile(matrix_file)), "dgCMatrix")
if (!identical(dim(expression), c(length(genes), length(cells)))) {
  stop("Matrix dimensions do not match genes.tsv and cells.tsv.", call. = FALSE)
}
if (!identical(metadata$cell, cells)) {
  stop("Metadata and matrix cell order differ.", call. = FALSE)
}
rownames(expression) <- genes
colnames(expression) <- cells
if (cell_limit > 0L) {
  if (cell_limit > ncol(expression)) stop("Cell limit exceeds input size.")
  keep <- seq_len(cell_limit)
  expression <- expression[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  cells <- cells[keep]
}
coordinates <- as.matrix(metadata[, c("x_um", "y_um"), drop = FALSE])
storage.mode(coordinates) <- "double"
rownames(coordinates) <- cells
group <- droplevels(factor(metadata$group, levels = sort(unique(metadata$group))))
names(group) <- cells

lr <- readRDS(file.path(root, "results", "spatialess_raw", "computable_lr.rds"))
if (lr_limit > 0L) {
  if (lr_limit > nrow(lr)) stop("LR limit exceeds available interactions.")
  lr <- lr[seq_len(lr_limit), , drop = FALSE]
}
data(CellChatDB.human, package = "CellChat")
message(sprintf(
  "Running RCTD-labelled %s: %d genes, %d cells, %d groups, %d LR",
  method, nrow(expression), ncol(expression), nlevels(group), nrow(lr)
))
start <- proc.time()[["elapsed"]]
result <- spatial_v3_exact(
  expression = expression,
  coordinates = coordinates,
  group = group,
  lr = lr,
  complex = CellChatDB.human$complex,
  cofactor = CellChatDB.human$cofactor,
  ratio = 1,
  tol = 5,
  interaction.range = 30,
  contact.range = 10,
  scale.distance = 20,
  min.percent = 0,
  min.cells.sr = 1L,
  nboot = nboot,
  seed.use = seed
)
elapsed <- proc.time()[["elapsed"]] - start

saveRDS(result$records, file.path(out_dir, "records.rds"), compress = FALSE)
summary <- data.frame(
  method = method,
  label_source = "RAW_RCTD_shared",
  genes = nrow(expression),
  cells = ncol(expression),
  groups = nlevels(group),
  lr = nrow(lr),
  nboot = nboot,
  seed = seed,
  graph_edges = length(result$graph$neighbors),
  active_records = nrow(result$records),
  significant_records = sum(result$records$pvalue < 0.05),
  inference_seconds = elapsed,
  stringsAsFactors = FALSE
)
write.table(
  summary, file.path(out_dir, "summary.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
print(summary)

