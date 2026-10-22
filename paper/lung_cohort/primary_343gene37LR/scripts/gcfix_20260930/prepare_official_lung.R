args <- commandArgs(TRUE); root <- normalizePath(args[1])
suppressPackageStartupMessages({library(Matrix); library(CellChat)})
env <- new.env(); data('CellChatDB.human',package='CellChat',envir=env); db <- env$CellChatDB.human
old <- '/data/dzf/GSE250346/benchmarks/multisample_v3_exact'
m <- read.delim(file.path(root,'inputs/manifest.tsv')); m <- m[m$family=='lung_paper',]
rows <- list()
for(i in seq_len(nrow(m))) {
  a <- readRDS(m$input[i]); sample <- sub('lung_paper_','',m$case[i])
  genes <- intersect(CellChat::extractGene(db),rownames(a$expression))
  stopifnot(length(genes)==99L,identical(a$complex,db$complex),identical(a$cofactor,db$cofactor))
  a$expression <- a$expression[genes,,drop=FALSE]
  a$lr <- readRDS(file.path(old,'computable_lr.rds'))
  id <- paste0('lung_official_',sample); path <- file.path(root,'inputs',paste0(id,'.rds'))
  saveRDS(a,path,compress=FALSE)
  ref <- file.path(old,'official_independent',sample,'records.rds'); kind <- 'unmodified_official'
  if(!file.exists(ref)) {
    ref <- file.path(old,'official_patched_retry',sample,'records.rds'); kind <- 'dims_only_official'
  }
  if(!file.exists(ref)) {ref <- '';kind <- 'official_failed_no_completed_reference'}
  rows[[i]] <- data.frame(case=id,N=ncol(a$expression),genes=nrow(a$expression),input=path,
    input_md5=unname(tools::md5sum(path)),historical=file.path(old,'full_sample_results',paste0(sample,'.rds')),
    official=ref,official_kind=kind,source=m$input[i])
}
write.table(do.call(rbind,rows),file.path(root,'inputs/official_lung_manifest.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'logs/official_lung_prepare_session.txt'))
