#!/usr/bin/env Rscript
# Serialize an already saved score matrix; no inference or model fitting.
a <- commandArgs(trailingOnly=TRUE)
if(length(a)!=2L) stop('Usage: Rscript tools/export_saved_scores.R PREPARED_RDS OUTPUT_TSV_GZ')
stopifnot(!file.exists(a[2]))
x <- readRDS(a[1]); stopifnot(all(dim(x$score)==c(45L,50040L)),!anyNA(x$score),all(x$score>=0))
con <- gzfile(a[2],'wt'); on.exit(close(con))
old <- options(digits=17); on.exit(options(old),add=TRUE)
write.table(data.frame(sample_id=rownames(x$score),x$score,check.names=FALSE),con,sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
