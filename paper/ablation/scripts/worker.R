args <- commandArgs(TRUE)
if (length(args) != 6L) stop('worker.R root case variant nboot repetition phase')
root <- normalizePath(args[1]); id <- args[2]; variant <- args[3]
nboot <- as.integer(args[4]); repetition <- as.integer(args[5]); phase <- args[6]
.libPaths(c(file.path(root, if (variant == 'vanilla') 'lib_vanilla' else 'lib'), .libPaths()))
suppressPackageStartupMessages(library(SpatialESS))
suppressPackageStartupMessages(library(Matrix))
stopifnot(as.character(packageVersion('SpatialESS')) == '0.1.4')
if (requireNamespace('future', quietly = TRUE)) future::plan(future::sequential)
options(digits = 7, scipen = 0)
storage <- switch(variant, coo = 'coo', dense = 'dense', 'csr')
options(spatialess.ablation.storage = storage,
        spatialess.ablation.hash = variant != 'no_hash',
        spatialess.ablation.reuse = variant != 'no_reuse')
tag <- sprintf('%s_%s_b%d_rep%d', id, variant, nboot, repetition)
out <- file.path(root, 'results', phase, tag)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out, 'summary.tsv'))) stop('Completed output already exists: ', out)
cat('START', tag, format(Sys.time()), 'package:', find.package('SpatialESS'), '\n')
a <- readRDS(file.path(root, 'inputs', paste0(id, '.rds')))
if (phase == 'tiny') {
  count <- min(250L, ncol(a$expression))
  a$expression <- a$expression[, seq_len(count), drop = FALSE]
  a$coordinates <- a$coordinates[seq_len(count), , drop = FALSE]
  a$group <- droplevels(factor(a$group[seq_len(count)]))
  a$lr <- a$lr[seq_len(min(30L, nrow(a$lr))), , drop = FALSE]
  a$min.percent <- 0.1; a$min.cells.sr <- 5L
  a$nboot <- nboot
  set.seed(a$seed.use)
  if (requireNamespace('future.apply', quietly = TRUE)) {
    invisible(future.apply::future_lapply(1L, identity, future.seed = TRUE))
  }
  permutations <- replicate(nboot, sample.int(count, count))
  if (nboot == 1L) permutations <- matrix(permutations, nrow = count)
} else {
  a$nboot <- nboot
  permutations <- readRDS(file.path(root, 'inputs', paste0(id, '_permutations.rds')))
  permutations <- permutations[, seq_len(nboot), drop = FALSE]
}
options(spatialess.ablation.permutations = permutations)
gc()
start <- proc.time()[['elapsed']]
result <- do.call(spatialess, a)
elapsed <- proc.time()[['elapsed']] - start
g <- result$graph
if (storage == 'csr') {
  stopifnot(length(g$offsets) == ncol(a$expression) + 1L,
            tail(g$offsets, 1) == length(g$neighbors),
            length(g$neighbors) == length(g$distances), all(diff(g$offsets) >= 0))
  E <- length(g$neighbors)
  primary_bytes <- sum(vapply(g[c('offsets', 'neighbors', 'distances')],
                             function(x) length(x) * if (typeof(x) == 'integer') 4 else 8, numeric(1)))
} else if (storage == 'coo') {
  E <- length(g$neighbors)
  stopifnot(E == length(g$senders), E == length(g$distances))
  primary_bytes <- E * 16
} else {
  E <- g$directed_edges
  stopifnot(identical(dim(g$dense_distances), rep(ncol(a$expression), 2L)))
  primary_bytes <- length(g$dense_distances) * 8
}
d <- result$diagnostics
summary <- data.frame(case = id, variant = variant, phase = phase, repetition = repetition,
  cells = ncol(a$expression), genes = nrow(a$expression), lr = nrow(result$lr),
  nboot = nboot, seed = a$seed.use, directed_edges = E,
  inference_seconds = elapsed,
  graph_seconds = if (is.null(d$graph_seconds)) NA_real_ else d$graph_seconds,
  observed_seconds = if (is.null(d$observed_seconds)) NA_real_ else d$observed_seconds,
  permutation_seconds = if (is.null(d$permutation_seconds)) NA_real_ else d$permutation_seconds,
  output_seconds = if (is.null(d$output_seconds)) NA_real_ else d$output_seconds,
  primary_graph_array_bytes = primary_bytes, graph_object_bytes = as.numeric(object.size(g)),
  active_records = nrow(result$records), significant_lt005 = sum(result$records$pvalue < .05),
  significant_le005 = sum(result$records$pvalue <= .05),
  version = as.character(packageVersion('SpatialESS')), status = 'completed')
saveRDS(result$records, file.path(out, 'records.rds'), compress = FALSE)
saveRDS(result[c('parameters', 'diagnostics')], file.path(out, 'contract.rds'))
if (phase %in% c('tiny', 'preflight') && repetition == 0L && nboot <= 100L) {
  saveRDS(g, file.path(out, 'graph.rds'), compress = FALSE)
}
write.table(summary, file.path(out, 'summary.tsv'), sep = '\t', quote = TRUE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
print(summary)
cat('FINISH', tag, format(Sys.time()), '\n')
