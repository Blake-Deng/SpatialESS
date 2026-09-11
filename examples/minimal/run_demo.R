#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
h5ad <- file.path(root, "examples", "minimal", "demo_sparkle.h5ad")
output_dir <- file.path(root, "examples", "minimal", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

python <- Sys.getenv("PYTHON", unset = Sys.which("python3"))
suppressPackageStartupMessages(library(SpatialESSSPARKLE))

lr <- data.frame(
  ligand = "L", receptor = "R", interaction_name = "L_R",
  pathway_name = "demo", annotation = "Secreted Signaling",
  agonist = "", antagonist = "", co_A_receptor = "",
  co_I_receptor = "", row.names = "L_R", stringsAsFactors = FALSE
)
complex <- data.frame(subunit_1 = character(), row.names = character())
cofactor <- data.frame(cofactor1 = character(), row.names = character())

result <- spatialess_from_h5ad(
  h5ad = h5ad,
  group_col = "annotation",
  coordinate_cols = c("x_um", "y_um"),
  lr = lr,
  complex = complex,
  cofactor = cofactor,
  gene_list = c("L", "R"),
  normalization = "log1p_cp10k",
  python = python,
  output_dir = file.path(output_dir, "bundle"),
  keep_bundle = TRUE,
  ratio = 1,
  tol = 0,
  interaction.range = 30,
  contact.range = 10,
  scale.distance = 1,
  min.percent = 0,
  min.cells.sr = 1L,
  nboot = 10L,
  seed.use = 20260904L
)
saveRDS(result, file.path(output_dir, "demo_result.rds"), compress = FALSE)
utils::write.table(
  result$records,
  file.path(output_dir, "demo_records.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(result$bridge)
print(result$records)
