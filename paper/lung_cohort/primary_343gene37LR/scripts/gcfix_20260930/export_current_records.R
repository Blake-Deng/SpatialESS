args <- commandArgs(TRUE)
root <- normalizePath(args[1]); dest <- normalizePath(args[2])
.libPaths(c(file.path(root, 'lib_candidate'), .libPaths()))
suppressPackageStartupMessages(library(Matrix))
m <- read.delim(file.path(dest, 'build_inputs.tsv'), stringsAsFactors=FALSE)
main <- file.path(dest, 'repositories/SpatialESS')
dir.create(file.path(main,'results/records'), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(main,'results/contracts'), recursive=TRUE, showWarnings=FALSE)
extract <- function(path) {
  x <- readRDS(path)
  if (is.list(x) && !is.data.frame(x)) x <- x$records
  stopifnot(is.data.frame(x)); x
}
key <- function(x) do.call(paste, c(x[c('sender_group_name','receiver_group_name','interaction_name')],sep='\r'))
parameters <- list(); comparisons <- list(); maxcols <- list()
for (i in seq_len(nrow(m))) {
  row <- m[i,]; records <- extract(file.path(row$run_dir,'records.rds'))
  contract <- readRDS(file.path(row$run_dir,'contract.rds'))
  saveRDS(records, file.path(main,'results/records',paste0(row$case,'.rds')))
  saveRDS(contract, file.path(main,'results/contracts',paste0(row$case,'.rds')))
  parameters[[i]] <- data.frame(case=row$case, parameter=names(contract$parameters),
    value=vapply(contract$parameters, function(x)paste(x,collapse=';'), character(1)))
  if (nzchar(row$official) && file.exists(row$official)) {
    ref <- extract(row$official); a <- key(ref); b <- key(records)
    stopifnot(!anyDuplicated(a), !anyDuplicated(b))
    shared <- intersect(a,b); united <- union(a,b); ai <- match(shared,a); bi <- match(shared,b)
    significant <- function(x, keys, strict) {
      ans <- setNames(rep(FALSE,length(united)),united)
      ans[keys] <- if (strict) x$pvalue < .05 else x$pvalue <= .05
      ans
    }
    jaccard <- if(length(united)) length(shared)/length(united) else 1
    pd <- if(length(shared)) max(abs(ref$probability[ai]-records$probability[bi])) else NA_real_
    pe <- if(length(shared)) mean(ref$pvalue[ai]==records$pvalue[bi]) else NA_real_
    pdelta <- if(length(shared)) max(abs(ref$pvalue[ai]-records$pvalue[bi])) else NA_real_
    lt <- sum(significant(ref,a,TRUE)!=significant(records,b,TRUE))
    le <- sum(significant(ref,a,FALSE)!=significant(records,b,FALSE))
    rho <- if(length(shared)>1L && sd(ref$probability[ai])>0 && sd(records$probability[bi])>0)
      cor(ref$probability[ai],records$probability[bi]) else NA_real_
    comparisons[[length(comparisons)+1L]] <- data.frame(case=row$case,
      comparator=row$official_kind, comparator_version='0.1.0',
      comparator_commit='be88f300464cf970b36a6165ea7d2e41a47a7511',
      ESS_version='0.1.4', reference_active=nrow(ref), ESS_active=nrow(records),
      active_jaccard=jaccard, probability_pearson=rho,
      max_probability_difference=pd, exact_pvalue_fraction=pe,
      max_pvalue_difference=pdelta, significant_disagreements_lt005=lt,
      significant_disagreements_le005=le,
      strict_fidelity_pass=jaccard==1 && !is.na(pd) && pd<=1e-15 && pe==1 && lt==0 && le==0,
      reference_file=basename(dirname(row$official)), reference_md5=unname(tools::md5sum(row$official)))
  }
}
write.table(do.call(rbind,parameters),file.path(main,'results/tables/parameters.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
write.table(do.call(rbind,comparisons),file.path(main,'results/tables/official_fidelity.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(main,'results/environment/export_session.txt'))
cat('Exported',nrow(m),'current inference configurations and',length(comparisons),'official comparisons.\n')
