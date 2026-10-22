args <- commandArgs(TRUE)
root <- normalizePath(args[1])
.libPaths(c(file.path(root, 'lib_candidate'), .libPaths()))
suppressPackageStartupMessages({library(Matrix); library(SpatialESS)})
base <- Sys.getenv('SPATIALESS_WORKSPACE', '/home/dzf/cellchat_acceleration')
compat <- file.path(base, 'releases/SpatialESS_compatibility_20260911')
rows <- list()
emit <- function(id, a, source, historical='', official='', repeats=1L, family=id) {
  stopifnot(ncol(a$expression)==nrow(a$coordinates), length(a$group)==ncol(a$expression))
  path <- file.path(root, 'inputs', paste0(id,'.rds'))
  saveRDS(a, paste0(path,'.tmp'), compress=FALSE)
  stopifnot(file.rename(paste0(path,'.tmp'), path))
  rows[[id]] <<- data.frame(case=id,family=family,input=path,source=source,
    historical=historical,official=official,N=ncol(a$expression),genes=nrow(a$expression),
    nboot=a$nboot,seed=a$seed.use,repeats=repeats,input_md5=unname(tools::md5sum(path)))
  write.table(do.call(rbind,rows), file.path(root,'inputs/manifest.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
  cat('PREPARED',id,'\n'); flush.console()
}
spatial_args <- function(f, n, nboot, seed, lr=f$lr) {
  ix <- f$selection_order[seq_len(n)]
  list(expression=f$expression[,ix,drop=FALSE],coordinates=f$coordinates_um[ix,,drop=FALSE],
    group=droplevels(f$group[ix]),lr=lr,complex=f$db$complex,cofactor=f$db$cofactor,
    ratio=f$spatial_factors$ratio,tol=f$spatial_factors$tol,
    interaction.range=f$parameters$interaction_range_um,scale.distance=f$parameters$scale_distance,
    contact.range=f$parameters$contact_range_um,min.percent=0,min.cells.sr=1L,
    nboot=nboot,seed.use=seed)
}
cosmx <- '/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead'
path <- file.path(cosmx,'secact_cosmx_v3_common_fixture.rds'); f <- readRDS(path)
for (n in c(250L,1000L,5000L,10000L,25000L,50000L,100000L,200000L,443515L)) {
  old <- file.path(cosmx,sprintf('cells%06d_full_spatialess_perm00020/records.rds',n))
  official <- file.path(cosmx,sprintf('cells%06d_full_spatialcellchat_v3_perm00020/records.rds',n))
  emit(paste0('cosmx_',n),spatial_args(f,n,20L,20260803L,f$lr_full),path,old,
    if(file.exists(official)) official else '',3L,'cosmx')
}
rm(f); gc()
for (case in c('center_1000','center_5000','stratified_1000','stratified_5000')) {
  d <- file.path(compat,'validation/dual_track',case); path <- file.path(d,'fixture.rds'); f <- readRDS(path)
  for (filter in c('default','relaxed')) {
    settings <- if(filter=='default') list(min.percent=.1,min.cells.sr=5L) else list(min.percent=0,min.cells.sr=1L)
    a <- c(list(expression=f$expression,coordinates=f$coordinates,group=f$group,lr=f$lr,
      complex=f$db$complex,cofactor=f$db$cofactor),f$parameters,settings)
    emit(paste('cervical',case,filter,sep='_'),a,path,file.path(d,'results',filter,'ess.rds'),
      file.path(d,'results',filter,'official_records.rds'),3L,'cervical')
  }
}
db <- readRDS(file.path(base,'releases/SpatialESS_main_SPARKLE_20260910/ovarian_historical_database.rds'))
for (condition in c('raw','sparkle')) for(n in c(1000L,5000L,16247L)) {
  d <- file.path(compat,'validation/bridge/runs',paste(condition,n,sep='_'),'direct_main')
  path <- file.path(d,'input.rds'); f <- readRDS(path)
  a <- c(f[c('expression','coordinates','group')],list(lr=readRDS(file.path(d,'contract.rds'))$lr,
    complex=db$complex,cofactor=db$cofactor),f$parameters)
  emit(paste('ovarian',condition,n,sep='_'),a,path,file.path(d,'records.rds'),'',3L,'ovarian')
}
for (sample in c('6_BCC','35_BCC','4_SCC','23_SCC')) {
  d <- file.path('/data/dzf/GSE306130/benchmarks/spatial_v3_headtohead',sample)
  path <- file.path(d,'exact_fixture.rds'); f <- readRDS(path)
  for(n in if(sample=='35_BCC') c(1000L,5000L,10000L,25000L,50000L) else 1000L) {
    emit(paste('gse306130',sample,n,sep='_'),spatial_args(f,n,20L,20260828L),path,
      file.path(d,sprintf('cells%06d_exact_perm00020/records.rds',n)),
      file.path(d,sprintf('cells%06d_official_perm00020/records.rds',n)),1L,'gse306130')
  }
}
d <- '/data/dzf/GSE313006/benchmarks/spatial_v3_headtohead'
path <- file.path(d,'gse313006_tma1_exact_fixture.rds'); f <- readRDS(path)
for(b in c(20L,100L)) emit(paste0('gse313006_1000_boot',b),spatial_args(f,1000L,b,20260820L),path,
  file.path(d,sprintf('cells001000_exact_perm%05d/records.rds',b)),
  file.path(d,sprintf('cells001000_official_perm%05d/records.rds',b)),1L,'gse313006')
rm(f); gc()

# Preserve both the paper's normalized 30-um workflow and the collaborator's
# raw-count/default 250-um workflow; their outputs are not interchangeable.
suppressPackageStartupMessages(library(SeuratObject))
path <- '/data/qqin/testess/workspace/GSE250346_lung_fibrosis/GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds'
obj <- readRDS(path); meta <- obj@meta.data
x <- methods::as(obj@assays$RNA@data,'dgCMatrix'); counts <- methods::as(obj@assays$RNA@counts,'dgCMatrix')
stopifnot(identical(colnames(x),rownames(meta)))
db_env <- new.env(); load('/data/qqin/testess/workspace/_shared/cellchatdb/CellChatDB.human.rda',envir=db_env)
db <- db_env$CellChatDB.human
old_root <- '/data/dzf/GSE250346/benchmarks/multisample'
lr <- readRDS(file.path(old_root,'computable_lr.rds'))
samples <- sort(unique(as.character(meta$sample)))
write.table(unique(meta[,c('sample','patient','disease_status','sample_affect','tma','run')]),
  file.path(root,'inputs/lung_metadata.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
for(i in seq_along(samples)) {
  sample <- samples[i]; ix <- which(as.character(meta$sample)==sample & !is.na(meta$final_CT) &
    nzchar(as.character(meta$final_CT)) & is.finite(meta$x_centroid) & is.finite(meta$y_centroid))
  coord <- as.matrix(meta[ix,c('x_centroid','y_centroid')]); storage.mode(coord) <- 'double'
  group <- droplevels(factor(as.character(meta$final_CT[ix]))); names(group) <- colnames(x)[ix]
  a <- list(expression=x[,ix,drop=FALSE],coordinates=coord,group=group,
    lr=lr,complex=db$complex,cofactor=db$cofactor,ratio=1,tol=5,interaction.range=30,
    contact.range=10,scale.distance=20,min.percent=0,min.cells.sr=1L,nboot=20L,seed.use=20260831L)
  emit(paste0('lung_paper_',sample),a,path,file.path(old_root,'full_sample_results',paste0(sample,'.rds')),
    '',1L,'lung_paper')
  a$expression <- counts[,ix,drop=FALSE]; a$lr <- db$interaction; a$lr_missing <- 'filter'
  a$interaction.range <- 250; a$scale.distance <- .01; a$min.percent <- .1; a$min.cells.sr <- 5L
  a$nboot <- 100L; a$seed.use <- 1L+i
  # Leave scale.distance at the actual package default (0.01).
  a$scale.distance <- NULL
  emit(paste0('lung_teacher_',sample),a,path,'','',if(i %in% c(6L,33L,45L)) 3L else 1L,'lung_teacher')
}
writeLines(capture.output(sessionInfo()),file.path(root,'logs/prepare_session.txt'))
writeLines('complete',file.path(root,'inputs/PREPARED'))
