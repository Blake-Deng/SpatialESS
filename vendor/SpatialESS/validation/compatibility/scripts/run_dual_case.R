args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
case <- args[[1]]
root <- "/home/dzf/cellchat_acceleration/releases/SpatialESS_compatibility_20260911"
old <- "/home/dzf/cellchat_acceleration/validation/stratified_V3_patched_20260910"
.libPaths(c(file.path(root, "lib_new"), file.path(old, "lib_patched"), .libPaths()))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SpatialESS))
suppressPackageStartupMessages(library(SpatialCellChat))
stopifnot(as.character(packageVersion("SpatialESS")) == "0.1.3")
d <- file.path(root, "validation/dual_track", case)
out <- file.path(d, "results")
dir.create(out, showWarnings = FALSE)
options(future.globals.maxSize = Inf)
options(error = function() {
  writeLines(geterrmessage(), file.path(out, "ERROR.txt"))
  traceback(20)
  quit(status = 1L, save = "no")
})
f <- readRDS(file.path(d, "fixture.rds"))
chat <- readRDS(file.path(d, "chat_input.rds"))
p <- f$parameters
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(c("Official V3 commit be88f300464cf970b36a6165ea7d2e41a47a7511",
  "Only official source patch: dims=c(NC,NC)", paste("SpatialESS", find.package("SpatialESS")),
  "Cell stage 8 workers, group stage sequential, BLAS/OMP 1; shared-server timing",
  "Group labels are public clusters, not collaborator biological annotations"), file.path(out, "provenance.txt"))
cell_args <- p[c("tol", "interaction.range", "contact.range", "scale.distance", "Kh", "n",
                 "use.AGAN", "contact.dependent", "contact.dependent.forced")]
# Reuse a frozen cell-stage result only after checking the full input contract.
cached <- "/home/dzf/cellchat_acceleration/validation/collaborator_cervical_20260911/results/official_normalized/chat_cell_stage.rds"
if (case == "stratified_1000") {
  prior <- readRDS(cached)
  stopifnot(identical(prior@data.signaling, chat@data.signaling),
    identical(prior@images$coordinates, chat@images$coordinates),
    identical(prior@idents, chat@idents), identical(prior@DB, chat@DB),
    identical(prior@LR$LRsig, chat@LR$LRsig),
    isTRUE(all.equal(prior@options$parameter[names(cell_args)], cell_args)))
  chat <- prior
  cell_seconds <- NA_real_
  cell_source <- "frozen_20260911_checked_input_contract"
} else {
  future::plan(future::multisession, workers = 8L)
  start <- proc.time()[["elapsed"]]
  chat <- do.call(SpatialCellChat::computeCommunProb,
    c(list(object = chat, LR.use = f$lr, raw.use = TRUE, distance.use = TRUE), cell_args))
  cell_seconds <- proc.time()[["elapsed"]] - start
  cell_source <- "new_run"
}
future::plan(future::sequential)
saveRDS(chat, file.path(out, "official_cell_stage.rds"))
cell_summary <- data.frame(case = case, cells = ncol(f$expression), lr = nrow(f$lr),
  positive_cell_lr_links = sum(vapply(chat@net$tmp$prob.cell, function(x) sum(x@x > 0), numeric(1))),
  cell_seconds = cell_seconds, source = cell_source)
write.table(cell_summary, file.path(out, "cell_stage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (filter_name in c("default", "relaxed")) {
  sub <- file.path(out, filter_name)
  dir.create(sub, showWarnings = FALSE)
  settings <- if (filter_name == "default") list(min.percent = .1, min.cells.sr = 5L) else list(min.percent = 0, min.cells.sr = 1L)
  t <- proc.time()[["elapsed"]]
  official <- do.call(SpatialCellChat::computeAvgCommunProb,
    c(list(object = chat, avg.type = "avg", do.permutation = TRUE, nboot = p$nboot,
           seed.use = p$seed.use, colocalization.use = FALSE), settings))
  official_seconds <- proc.time()[["elapsed"]] - t
  saveRDS(official, file.path(sub, "official.rds"))
  t <- proc.time()[["elapsed"]]
  ess <- do.call(SpatialESS::spatialess, c(list(expression = f$expression, coordinates = f$coordinates,
    group = f$group, lr = f$lr, complex = f$db$complex, cofactor = f$db$cofactor), p, settings))
  ess_seconds <- proc.time()[["elapsed"]] - t
  saveRDS(ess, file.path(sub, "ess.rds"))
  prob <- official@net[["prob", exact = TRUE]]
  pval <- official@net[["pval", exact = TRUE]]
  ix <- which(prob > 0, arr.ind = TRUE)
  a <- data.frame(sender_group_name = dimnames(prob)[[1]][ix[, 1]],
    receiver_group_name = dimnames(prob)[[2]][ix[, 2]], interaction_name = dimnames(prob)[[3]][ix[, 3]],
    probability = prob[ix], pvalue = pval[ix])
  b <- ess$records
  keys <- c("sender_group_name", "receiver_group_name", "interaction_name")
  m <- merge(a, b[, c(keys, "probability", "pvalue")], by = keys, all = TRUE,
             suffixes = c("_official", "_ess"))
  both <- complete.cases(m)
  summary <- data.frame(case = case, filter = filter_name, cells = ncol(f$expression),
    groups = nlevels(f$group), lr = nrow(f$lr), min_percent = settings$min.percent,
    min_cells_sr = settings$min.cells.sr, official_active = nrow(a), ess_active = nrow(b),
    active_jaccard = if (nrow(m)) sum(both)/nrow(m) else 1,
    max_abs_probability_difference = if (any(both)) max(abs(m$probability_official[both] - m$probability_ess[both])) else NA_real_,
    pvalue_exact_fraction = if (all(both) && any(both)) mean(m$pvalue_official == m$pvalue_ess) else NA_real_,
    official_significant_lt_005 = sum(a$pvalue < .05), ess_significant_lt_005 = sum(b$pvalue < .05),
    significance_disagreements = sum((!is.na(m$pvalue_official) & m$pvalue_official < .05) !=
                                    (!is.na(m$pvalue_ess) & m$pvalue_ess < .05)),
    official_group_seconds = official_seconds, ess_seconds = ess_seconds)
  summary$pass <- with(summary, active_jaccard == 1 && max_abs_probability_difference <= 1e-15 &&
    pvalue_exact_fraction == 1 && significance_disagreements == 0)
  saveRDS(a, file.path(sub, "official_records.rds"))
  write.table(m, file.path(sub, "records_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary, file.path(sub, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(summary)
}
writeLines("completed; consult individual fidelity gates", file.path(out, "COMPLETED.txt"))
