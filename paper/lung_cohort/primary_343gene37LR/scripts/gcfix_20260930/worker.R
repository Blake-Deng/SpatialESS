args <- commandArgs(TRUE)
root <- normalizePath(args[1]); case <- args[2]; variant <- args[3]; repetition <- as.integer(args[4])
.libPaths(c(file.path(root,paste0('lib_',variant)),.libPaths()))
suppressPackageStartupMessages({library(Matrix); library(SpatialESS)})
out <- file.path(root,'results',case,paste0(variant,'_',repetition))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
a <- readRDS(file.path(root,'inputs',paste0(case,'.rds')))
cat('START',case,variant,repetition,format(Sys.time()),'\n'); flush.console()
gc()
tm <- system.time(result <- do.call(spatialess,a))
rng <- .Random.seed
g <- result$graph
stopifnot(length(g$offsets)==ncol(a$expression)+1L,tail(g$offsets,1)==length(g$neighbors),
  length(g$neighbors)==length(g$distances),all(diff(g$offsets)>=0))
saveRDS(result$records,file.path(out,'records.rds'),compress=FALSE)
saveRDS(rng,file.path(out,'rng.rds'))
saveRDS(result[c('lr','parameters','diagnostics','lr_filter')],file.path(out,'contract.rds'))
if(!is.null(result$lr_filter)) {
  write.table(result$lr_filter$audit,file.path(out,'lr_audit.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
  write.table(result$lr_filter$missing_genes,file.path(out,'missing_genes.tsv'),sep='\t',row.names=FALSE,quote=TRUE)
}
s <- data.frame(case=case,variant=variant,repetition=repetition,version=as.character(packageVersion('SpatialESS')),
  cells=ncol(a$expression),genes=nrow(a$expression),groups=nlevels(factor(a$group)),lr=nrow(result$lr),
  nboot=a$nboot,seed=a$seed.use,inference_seconds=unname(tm['elapsed']),
  graph_edges=length(g$neighbors),graph_bytes=as.numeric(object.size(g)),
  active=nrow(result$records),significant_lt005=sum(result$records$pvalue<.05),
  significant_le005=sum(result$records$pvalue<=.05),status='completed')
write.table(s,file.path(out,'summary.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
print(s)
writeLines('complete',file.path(out,'COMPLETE'))
