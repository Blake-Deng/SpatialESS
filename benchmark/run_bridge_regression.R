#!/usr/bin/env Rscript
# Installed old/new libraries are selected by the invoking isolated process.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop("Usage: run_bridge_regression.R ROOT CONDITION CELLS MODE H5AD METADATA")
root <- normalizePath(args[1]); condition <- args[2]; cells <- as.integer(args[3]); mode <- args[4]
suppressPackageStartupMessages(library(SpatialESSSPARKLE))
suppressPackageStartupMessages(library(SpatialESS))
stopifnot(mode %in% c("old_bridge", "new_bridge", "direct_main"))
out <- file.path(root, "runs", paste(condition,cells,sep="_"), mode)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
db <- readRDS(file.path(root,"ovarian_historical_database.rds"))
source <- Sys.getenv("OVARIAN_SOURCE", "/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902")
genes <- readLines(file.path(source,"inputs/spatialess_rctd/genes.tsv"))
params <- list(ratio=1, tol=5, interaction.range=30, contact.range=10,
               scale.distance=20, min.percent=0, min.cells.sr=1L,
               nboot=100L, seed.use=20260903L)
gc(); start <- proc.time()[["elapsed"]]
if (mode != "direct_main") {
  a <- list(h5ad=args[5], metadata_file=args[6], group_col="group", coordinate_cols=c("x_um","y_um"),
    lr=db$interaction, complex=db$complex, cofactor=db$cofactor, gene_list=genes,
    python=Sys.getenv("PYTHON"), output_dir=file.path(out,"bundle"), keep_bundle=TRUE,
    normalization="log1p_cp10k")
  if (mode == "new_bridge") a$lr_missing <- "filter"
  result <- do.call(spatialess_from_h5ad,c(a,params))
  b <- read_spatialess_bundle(file.path(out,"bundle"))
} else {
  b <- read_spatialess_bundle(file.path(dirname(out),"new_bridge","bundle"))
  filter <- filter_cellchat_lr(db$interaction,rownames(b$expression),db$complex,db$cofactor)
  result <- do.call(spatialess,c(list(expression=b$expression,coordinates=b$coordinates,group=b$group,
    lr=filter$lr,complex=db$complex,cofactor=db$cofactor),params))
}
elapsed <- proc.time()[["elapsed"]] - start
stopifnot(ncol(b$expression)==cells)
saveRDS(result$records,file.path(out,"records.rds"),compress=FALSE)
saveRDS(.Random.seed,file.path(out,"rng.rds"))
saveRDS(list(expression=b$expression,coordinates=b$coordinates,group=b$group,parameters=params),
        file.path(out,"input.rds"),compress=FALSE)
saveRDS(list(lr=result$lr,parameters=result$parameters,diagnostics=result$diagnostics),file.path(out,"contract.rds"))
if (!is.null(result$lr_filter)) saveRDS(result$lr_filter,file.path(out,"lr_filter.rds"))
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
s <- data.frame(condition=condition,cells=cells,mode=mode,
  bridge_version=as.character(packageVersion("SpatialESSSPARKLE")),engine_version=as.character(packageVersion("SpatialESS")),
  genes=nrow(b$expression),lr=nrow(result$lr),nboot=params$nboot,seed=params$seed.use,
  pipeline_seconds=elapsed,active=nrow(result$records),significant=sum(result$records$pvalue<0.05),status="completed")
write.table(s,file.path(out,"summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
print(s)
