#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    paste(
      "Usage: run_spatialess_from_h5ad.R",
      "<corrected.h5ad> <group_col> <x_col> <y_col> <output.rds>"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(CellChat)
  library(SpatialESSSPARKLE)
})

data(CellChatDB.human, package = "CellChat")
result <- spatialess_from_h5ad(
  h5ad = args[[1L]],
  group_col = args[[2L]],
  coordinate_cols = args[3:4],
  lr = CellChatDB.human$interaction,
  complex = CellChatDB.human$complex,
  cofactor = CellChatDB.human$cofactor,
  normalization = "log1p_cp10k",
  nboot = 100L,
  seed.use = 1L
)
saveRDS(result, args[[5L]], compress = FALSE)
print(result$bridge)
print(utils::head(result$records))
