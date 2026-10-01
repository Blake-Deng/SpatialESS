args <- commandArgs(TRUE)
if (length(args) != 2L) stop('Usage: Rscript run_prepared.R input.rds output_directory')
suppressPackageStartupMessages(library(SpatialESS))
stopifnot(as.character(packageVersion('SpatialESS')) == '0.1.4')
input <- normalizePath(args[1])
out <- args[2]
if (dir.exists(out) && length(list.files(out, all.files=TRUE, no..=TRUE))) {
  stop('Use an empty output directory to preserve completed results.')
}
dir.create(out, recursive=TRUE, showWarnings=FALSE)
a <- readRDS(input)
stopifnot(is.list(a), !is.null(names(a)), !is.null(a$seed.use))
elapsed <- system.time(result <- do.call(spatialess, a))
graph <- result$graph
stopifnot(length(graph$offsets) == ncol(a$expression) + 1L,
  tail(graph$offsets, 1) == length(graph$neighbors),
  length(graph$neighbors) == length(graph$distances))
saveRDS(result$records, file.path(out, 'records.rds'))
saveRDS(result[c('lr','parameters','diagnostics','lr_filter')], file.path(out, 'contract.rds'))
saveRDS(.Random.seed, file.path(out, 'rng.rds'))
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
summary <- data.frame(version=as.character(packageVersion('SpatialESS')),
  cells=ncol(a$expression), genes=nrow(a$expression), lr=nrow(result$lr),
  inference_seconds=unname(elapsed['elapsed']), graph_edges=length(graph$neighbors),
  active=nrow(result$records), significant_lt005=sum(result$records$pvalue<.05),
  input_md5=unname(tools::md5sum(input)), status='completed')
write.table(summary, file.path(out,'summary.tsv'), sep='\t', quote=TRUE, row.names=FALSE)
print(summary)
