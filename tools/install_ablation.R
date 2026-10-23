#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop('Usage: Rscript tools/install_ablation.R REPOSITORY_ROOT')
root <- normalizePath(args[1],mustWork=TRUE)
base <- file.path(root,'paper','ablation')
for (entry in list(c('vanilla','lib_vanilla'),c('instrumented','lib'))) {
  src <- file.path(base,'source',entry[1]); lib <- file.path(base,entry[2])
  dir.create(lib,recursive=TRUE,showWarnings=FALSE)
  status <- system2(file.path(R.home('bin'),'R'),c('CMD','INSTALL',paste0('--library=',shQuote(lib)),shQuote(src)))
  if (status!=0L) stop('Installation failed: ',src)
}
cat('Installed vanilla and instrumented into the libraries expected by worker.R\n')
