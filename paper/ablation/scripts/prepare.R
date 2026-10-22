args <- commandArgs(TRUE)
root <- normalizePath(args[1])
release <- normalizePath(args[2])
input_root <- normalizePath(args[3])
.libPaths(c(file.path(root, 'lib'), .libPaths()))
suppressPackageStartupMessages(library(SpatialESS))
suppressPackageStartupMessages(library(Matrix))
stopifnot(as.character(packageVersion('SpatialESS')) == '0.1.4')
if (requireNamespace('future', quietly = TRUE)) future::plan(future::sequential)
table <- read.delim(file.path(release, 'results/tables/inputs.tsv'))
cases <- c('cosmx_5000', 'cosmx_10000', 'cosmx_25000',
           'cervical_center_5000_default', 'cervical_stratified_5000_default',
           'cosmx_50000')
rows <- list()
for (id in cases) {
  expected <- table[table$case == id, , drop = FALSE]
  stopifnot(nrow(expected) == 1L)
  path <- file.path(input_root, expected$input_filename)
  actual <- strsplit(system2('sha256sum', shQuote(path), stdout = TRUE), ' ')[[1]][1]
  stopifnot(actual == expected$input_sha256)
  a <- readRDS(path)
  stopifnot(ncol(a$expression) == expected$cells,
            nrow(a$expression) == expected$genes,
            a$seed.use == expected$seed, a$nboot == expected$nboot)
  target <- file.path(root, 'inputs', paste0(id, '.rds'))
  if (!file.exists(target)) stopifnot(file.copy(path, target))
  max_boot <- if (id == 'cosmx_10000') 500L else a$nboot
  set.seed(a$seed.use)
  if (requireNamespace('future.apply', quietly = TRUE)) {
    invisible(future.apply::future_lapply(1L, identity, future.seed = TRUE))
  }
  permutations <- replicate(max_boot, sample.int(ncol(a$expression), ncol(a$expression)))
  storage.mode(permutations) <- 'integer'
  saveRDS(permutations, file.path(root, 'inputs', paste0(id, '_permutations.rds')),
          compress = FALSE)
  rows[[id]] <- data.frame(case = id, input_source = path,
    input_sha256 = actual, cells = ncol(a$expression), genes = nrow(a$expression),
    groups = nlevels(droplevels(factor(a$group))), lr = nrow(a$lr),
    nboot = a$nboot, seed = a$seed.use, max_frozen_permutations = max_boot,
    effective_radius = a$interaction.range + a$tol, status = 'PASS')
  cat('PREPARED', id, ncol(a$expression), nrow(a$lr), max_boot, '\n')
}
write.table(do.call(rbind, rows), file.path(root, 'reports/input_manifest.tsv'),
            sep = '\t', quote = TRUE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, 'logs/preparation_sessionInfo.txt'))
