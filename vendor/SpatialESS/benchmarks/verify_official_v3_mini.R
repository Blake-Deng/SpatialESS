#!/usr/bin/env Rscript
# Requires a separately installed, dimension-only-patched SpatialCellChat V3.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
out <- if (length(args) >= 2L) args[2] else tempfile("spatialess-v3-mini-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
workers <- as.integer(Sys.getenv("V3_WORKERS", "1"))
stopifnot(length(workers) == 1L, !is.na(workers), workers >= 1L)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SpatialESS))
suppressPackageStartupMessages(library(SpatialCellChat))
f <- readRDS(file.path(root, "tests/testthat/fixtures/v3_percent_boundary.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(capture.output(packageDescription("SpatialCellChat")), file.path(out, "official_description.txt"))
saveRDS(f, file.path(out, "input.rds"))
chat <- createSpatialCellChat(f$expression,
  meta = data.frame(labels = f$group, row.names = colnames(f$expression)),
  group.by = "labels", datatype = "spatial", coordinates = f$coordinates,
  spatial.factors = list(ratio = 1, tol = 5), do.sparse = TRUE)
chat@DB <- f$db
chat <- subsetData(chat)
chat@LR$LRsig <- f$lr
stopifnot(identical(chat@data.signaling, f$expression))
geometry <- SpatialCellChat::computeCellDistance(f$coordinates,
  interaction.range = 30, contact.range = 10, ratio = 1, tol = 5)
if (!identical(dim(geometry$d.spatial), rep(ncol(f$expression), 2))) {
  stop("Official matrix dimension mismatch. Install the isolated dims-only reference described in docs/V3_COMPATIBILITY.md; do not change its fraction filter.")
}
future::plan(future::multisession, workers = workers)
options(future.globals.maxSize = Inf)
chat <- computeCommunProb(chat, LR.use = f$lr, tol = 5,
  interaction.range = 30, contact.range = 10, scale.distance = 1)
future::plan(future::sequential)
official <- computeAvgCommunProb(chat, avg.type = "avg", do.permutation = TRUE,
  nboot = 100L, seed.use = 1L, colocalization.use = FALSE, min.percent = .1, min.cells.sr = 5L)
ess <- do.call(spatialess, c(list(expression = f$expression, coordinates = f$coordinates,
  group = f$group, lr = f$lr, complex = f$db$complex, cofactor = f$db$cofactor,
  min.percent = .1, min.cells.sr = 5L), f$parameters))
p <- official@net[["prob", exact = TRUE]]; q <- official@net[["pval", exact = TRUE]]
stopifnot(!is.null(p), !is.null(q))
ix <- which(p > 0, arr.ind = TRUE)
a <- data.frame(sender_group_name = dimnames(p)[[1]][ix[, 1]],
  receiver_group_name = dimnames(p)[[2]][ix[, 2]], interaction_name = dimnames(p)[[3]][ix[, 3]],
  probability = p[ix], pvalue = q[ix])
keys <- c("sender_group_name", "receiver_group_name", "interaction_name")
m <- merge(a, ess$records[, c(keys, "probability", "pvalue")], by = keys,
  all = TRUE, suffixes = c("_official", "_ess"))
both <- complete.cases(m)
summary <- data.frame(official_active = nrow(a), ess_active = nrow(ess$records),
  active_jaccard = if (nrow(m)) mean(both) else 1,
  max_abs_probability_difference = if (any(both)) max(abs(m$probability_official[both] - m$probability_ess[both])) else 0,
  pvalue_exact_fraction = if (all(both)) { if (nrow(m)) mean(m$pvalue_official == m$pvalue_ess) else 1 } else NA_real_,
  significance_disagreements = sum((!is.na(m$pvalue_official) & m$pvalue_official < .05) !=
                                  (!is.na(m$pvalue_ess) & m$pvalue_ess < .05)))
summary$pass <- with(summary, active_jaccard == 1 && max_abs_probability_difference <= 1e-15 &&
  pvalue_exact_fraction == 1 && significance_disagreements == 0)
saveRDS(official, file.path(out, "official.rds"))
saveRDS(ess, file.path(out, "ess.rds"))
write.table(m, file.path(out, "records_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary)
stopifnot(summary$pass)
cat("PASS: fresh official cell stage and group permutation versus SpatialESS.\n")
