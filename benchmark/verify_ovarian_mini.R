#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath(".")
f <- file.path(root,"validation/fixtures")
suppressPackageStartupMessages(library(SpatialESSSPARKLE))
db <- readRDS(file.path(f,"database.rds"))
genes <- readLines(file.path(f,"genes.tsv"))
for(condition in c("raw","sparkle")) {
  result <- spatialess_from_h5ad(
    h5ad=file.path(f,paste0(condition,"_1000.h5ad")),
    metadata_file=file.path(f,paste0(condition,"_1000_metadata.tsv")),
    group_col="group",coordinate_cols=c("x_um","y_um"),gene_list=genes,
    lr=db$interaction,complex=db$complex,cofactor=db$cofactor,
    python=Sys.getenv("PYTHON",unset=Sys.which("python3")),
    normalization="log1p_cp10k",lr_missing="filter",ratio=1,tol=5,
    interaction.range=30,contact.range=10,scale.distance=20,
    min.percent=0,min.cells.sr=1L,nboot=100L,seed.use=20260903L)
  expected <- readRDS(file.path(f,paste0(condition,"_1000_expected_records.rds")))
  stopifnot(identical(result$records,expected),result$bridge$computable_lr==1828L,
            result$bridge$excluded_lr==111L,!is.null(result$lr_filter$missing_genes))
  cat("PASS:",condition,"real ovarian 1k H5AD -> bridge -> bundled SpatialESS; records identical to old release.\n")
}
