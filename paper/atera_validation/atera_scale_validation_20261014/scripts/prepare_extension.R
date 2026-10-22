args <- commandArgs(TRUE)
root <- normalizePath(args[1])
base <- normalizePath(args[2])
suppressPackageStartupMessages({library(Matrix); library(data.table); library(jsonlite)})
.libPaths(c(file.path(root, 'lib'), .libPaths()))
stopifnot(as.character(packageVersion('SpatialESS')) == '0.1.4')
dir.create(file.path(root, 'reports'), recursive=TRUE, showWarnings=FALSE)

for (tissue in c('breast', 'cervical')) {
  src_dir <- file.path(base, 'inputs', tissue)
  dst_dir <- file.path(root, 'inputs', tissue)
  dir.create(dst_dir, recursive=TRUE, showWarnings=FALSE)
  source_rds <- file.path(src_dir, 'central_25000.rds')
  target_rds <- file.path(dst_dir, 'central_10000.rds')
  if (!file.exists(target_rds)) {
    x <- readRDS(source_rds)
    keep <- seq_len(10000L)
    y <- x
    y$expression <- x$expression[, keep, drop=FALSE]
    y$coordinates <- x$coordinates[keep, , drop=FALSE]
    y$group <- droplevels(x$group[keep])
    names(y$group) <- colnames(y$expression)
    saveRDS(y, target_rds, compress=FALSE)
    rm(x, y)
    gc()
  }
  actual <- readRDS(target_rds)
  reference <- readRDS(source_rds)
  previous <- readRDS(file.path(src_dir, 'central_5000.rds'))
  keep <- seq_len(10000L)
  stopifnot(identical(actual$expression, reference$expression[,keep,drop=FALSE]),
            identical(actual$coordinates, reference$coordinates[keep,,drop=FALSE]),
            identical(as.character(actual$group), as.character(reference$group[keep])),
            identical(colnames(actual$expression)[seq_len(5000L)], colnames(previous$expression)),
            identical(actual$expression[,seq_len(5000L),drop=FALSE], previous$expression),
            identical(actual$lr, reference$lr), nrow(actual$expression)==18028L,
            nrow(actual$lr)==3081L, actual$nboot==100L, actual$seed.use==1L)
  constants <- setdiff(names(reference), c('expression','coordinates','group'))
  stopifnot(identical(actual[constants], reference[constants]))
  meta <- fread(file.path(src_dir, 'cells.tsv'))
  coords <- as.matrix(meta[,.(x_centroid,y_centroid)])
  center <- colMeans(coords)
  indices <- order(rowSums(sweep(coords,2,center)^2), seq_len(nrow(meta)))[keep]
  stopifnot(identical(meta$cell_id[indices], colnames(actual$expression)))
  barcodes <- copy(meta[indices]); barcodes[,subset_order:=seq_len(.N)]
  fwrite(barcodes, file.path(dst_dir, 'central_10000_barcodes.tsv'), sep='\t')
  fwrite(barcodes[,.(cells=.N),by=group], file.path(dst_dir, 'central_10000_groups.tsv'), sep='\t')
  write_json(list(status='PASS',tissue=tissue,N=10000L,genes=18028L,lr=3081L,
                  nboot=100L,seed_use=1L,central_center=center,
                  exact_prefix_of_central_25000=TRUE,exact_extension_of_central_5000=TRUE,
                  full_coordinate_center_rule_verified=TRUE,
                  sampling='Nearest squared Euclidean distance to full eligible mean; original index breaks ties; no randomness.',
                  input_md5=unname(tools::md5sum(target_rds))),
             file.path(root,'reports',paste0(tissue,'_10000_input_validation.json')),
             pretty=TRUE,auto_unbox=TRUE,digits=NA)
  cat('prepared', tissue, target_rds, file.info(target_rds)$size, '\n')
  rm(actual,reference,previous,meta,coords,barcodes)
  gc()
}
writeLines(capture.output(sessionInfo()),file.path(root,'reports/preparation_sessionInfo.txt'))
