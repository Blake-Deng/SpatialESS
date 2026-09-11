#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SpatialCellChat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1L]] %in% c("raw", "sparkle")) {
  stop("Usage: run_official_rctd.R raw|sparkle", call. = FALSE)
}

method <- args[[1L]]
root <- "/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902"
input_dir <- file.path(root, "inputs", "spatialess_rctd")
nboot <- as.integer(Sys.getenv("RNA_LEAKAGE_NBOOT", "100"))
seed <- as.integer(Sys.getenv("RNA_LEAKAGE_SEED", "20260903"))
cell_limit <- as.integer(Sys.getenv("RNA_LEAKAGE_CELL_LIMIT", "0"))
lr_limit <- as.integer(Sys.getenv("RNA_LEAKAGE_LR_LIMIT", "0"))
run_tag <- Sys.getenv("RNA_LEAKAGE_RUN_TAG", "")
if (is.na(nboot) || nboot < 1L || is.na(seed) ||
    is.na(cell_limit) || cell_limit < 0L || is.na(lr_limit) || lr_limit < 0L) {
  stop("Invalid benchmark environment variables.", call. = FALSE)
}

suffix <- if (nzchar(run_tag)) paste0("_", run_tag) else ""
out_dir <- file.path(root, "results", paste0("official_rctd_", method, suffix))
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
meta <- data.frame(labels = group, row.names = cells)

lr <- readRDS(file.path(root, "results", "spatialess_raw", "computable_lr.rds"))
if (lr_limit > 0L) {
  if (lr_limit > nrow(lr)) stop("LR limit exceeds available interactions.")
  lr <- lr[seq_len(lr_limit), , drop = FALSE]
}
data(CellChatDB.human, package = "CellChat")
options(future.globals.maxSize = Inf)

message(sprintf(
  paste(
    "Official SpatialCellChat V3 RCTD-labelled %s:",
    "%d genes, %d cells, %d groups, %d LR, %d permutations"
  ),
  method, nrow(expression), ncol(expression), nlevels(group), nrow(lr), nboot
))
workflow_start <- proc.time()[["elapsed"]]
chat <- createSpatialCellChat(
  object = expression,
  meta = meta,
  group.by = "labels",
  datatype = "spatial",
  coordinates = coordinates,
  spatial.factors = list(ratio = 1, tol = 5),
  do.sparse = TRUE
)
chat@DB <- CellChatDB.human
chat <- subsetData(chat)
chat@LR$LRsig <- lr

probability_start <- proc.time()[["elapsed"]]
chat <- computeCommunProb(
  chat,
  LR.use = lr,
  raw.use = TRUE,
  Kh = 0.5,
  n = 1,
  distance.use = TRUE,
  interaction.range = 30,
  scale.distance = 20,
  use.AGAN = TRUE,
  contact.dependent = TRUE,
  contact.range = 10,
  contact.dependent.forced = FALSE
)
probability_seconds <- proc.time()[["elapsed"]] - probability_start

permutation_start <- proc.time()[["elapsed"]]
chat <- computeAvgCommunProb(
  chat,
  avg.type = "avg",
  min.percent = 0,
  min.cells.sr = 1L,
  do.permutation = TRUE,
  nboot = nboot,
  seed.use = seed,
  colocalization.use = FALSE
)
permutation_seconds <- proc.time()[["elapsed"]] - permutation_start

prob <- chat@net$prob
pvalue <- chat@net$pval
active <- which(prob > 0, arr.ind = TRUE)
records <- data.frame(
  sender_group_name = dimnames(prob)[[1L]][active[, 1L]],
  receiver_group_name = dimnames(prob)[[2L]][active[, 2L]],
  interaction_name = dimnames(prob)[[3L]][active[, 3L]],
  probability = prob[active],
  pvalue = pvalue[active],
  stringsAsFactors = FALSE
)
end_to_end_seconds <- proc.time()[["elapsed"]] - workflow_start

saveRDS(records, file.path(out_dir, "records.rds"), compress = FALSE)
summary <- data.frame(
  method = method,
  engine = "official_spatialcellchat_v3",
  label_source = "RAW_RCTD_shared",
  genes = nrow(expression),
  cells = ncol(expression),
  groups = nlevels(group),
  lr = nrow(lr),
  nboot = nboot,
  seed = seed,
  active_records = nrow(records),
  significant_records = sum(records$pvalue < 0.05),
  probability_seconds = probability_seconds,
  permutation_seconds = permutation_seconds,
  end_to_end_seconds = end_to_end_seconds,
  stringsAsFactors = FALSE
)
write.table(
  summary, file.path(out_dir, "summary.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
print(summary)
