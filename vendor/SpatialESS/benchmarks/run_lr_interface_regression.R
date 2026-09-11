#!/usr/bin/env Rscript
# One isolated R process per library/mode/case; no engine is patched in memory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: Rscript run_lr_interface_regression.R LIB CASE CELLS MODE OUT")
.libPaths(c(normalizePath(args[1]), .libPaths()))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SpatialESS))
case <- args[2]; cells <- as.integer(args[3]); mode <- args[4]; out <- args[5]
stopifnot(mode %in% c("old_manual", "new_explicit", "new_filter"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dataset_root <- Sys.getenv("LR_REGRESSION_COSMX", "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead")
cervical_root <- Sys.getenv("LR_REGRESSION_CERVICAL", "/home/dzf/cellchat_acceleration/external/cervical_xenium_20260910")
stratified_root <- Sys.getenv("LR_REGRESSION_STRATIFIED", "/data/dzf/SecAct_2026/results/spatialcellchat_v3_stratified_patched_ladder_20260909/cells010000")
if (case == "cervical") {
  path <- file.path(cervical_root, "results", sprintf("cells%06d", cells), "fixture.rds")
  f <- readRDS(path); expression <- f$expression; coordinates <- f$coordinates; group <- f$group
  prepared_lr <- f$lr
  parameters <- list(tol = 5, interaction.range = 30, contact.range = 10, nboot = 100, seed.use = 1)
} else {
  root <- if (case == "cosmx_stratified") stratified_root else dataset_root
  path <- file.path(root, "secact_cosmx_v3_common_fixture.rds")
  f <- readRDS(path); idx <- f$selection_order[seq_len(cells)]
  expression <- f$expression[, idx, drop = FALSE]
  coordinates <- f$coordinates_um[idx, , drop = FALSE]
  group <- droplevels(f$group[idx]); names(group) <- colnames(expression)
  prepared_lr <- f$lr_full
  parameters <- list(ratio = f$spatial_factors$ratio, tol = f$spatial_factors$tol,
    interaction.range = f$parameters$interaction_range_um, scale.distance = f$parameters$scale_distance,
    contact.range = f$parameters$contact_range_um, Kh = 0.5, n = 1, use.AGAN = TRUE,
    contact.dependent = TRUE, contact.dependent.forced = FALSE,
    min.percent = 0, min.cells.sr = 1, nboot = 20, seed.use = 20260803)
}
db_env <- new.env(); data("CellChatDB.human", package = "SpatialCellChat", envir = db_env)
db <- db_env$CellChatDB.human
genes <- rownames(expression)
# Freeze the old manual rule independently of the new helper.
parts <- function(name) {
  value <- if (name %in% rownames(db$complex)) {
    unlist(db$complex[name, grep("^subunit", names(db$complex)), drop = FALSE], use.names = FALSE)
  } else name
  value <- as.character(value); value[!is.na(value) & nzchar(value)]
}
ok <- vapply(seq_len(nrow(db$interaction)), function(i) {
  l <- parts(as.character(db$interaction$ligand[i])); r <- parts(as.character(db$interaction$receptor[i]))
  length(l) > 0L && length(r) > 0L && all(c(l, r) %in% genes)
}, logical(1))
manual_lr <- db$interaction[ok, , drop = FALSE]
stopifnot(identical(as.character(manual_lr$interaction_name), as.character(prepared_lr$interaction_name)))
base_args <- c(list(expression = expression, coordinates = coordinates, group = group,
                    complex = db$complex, cofactor = db$cofactor), parameters)
saveRDS(list(cell_ids = colnames(expression), genes = genes, group = group,
             coordinates = coordinates, lr = manual_lr, parameters = parameters), file.path(out, "input_contract.rds"))
writeLines(c(paste("fixture", normalizePath(path)), paste("library", find.package("SpatialESS")),
             paste("version", packageVersion("SpatialESS"))), file.path(out, "provenance.txt"))
gc()
start <- proc.time()[["elapsed"]]
if (mode == "old_manual") {
  result <- do.call(spatialess, c(base_args, list(lr = manual_lr)))
} else if (mode == "new_explicit") {
  selection <- filter_cellchat_lr(db$interaction, genes, db$complex, db$cofactor)
  stopifnot(identical(selection$lr, manual_lr))
  result <- do.call(spatialess, c(base_args, list(lr = selection$lr)))
} else {
  result <- do.call(spatialess, c(base_args, list(lr = db$interaction, lr_missing = "filter")))
  selection <- result$lr_filter
  stopifnot(identical(selection$lr, manual_lr))
}
seconds <- proc.time()[["elapsed"]] - start
saveRDS(result$records, file.path(out, "records.rds"), compress = FALSE)
saveRDS(.Random.seed, file.path(out, "rng.rds"))
saveRDS(list(lr = result$lr, parameters = result$parameters, diagnostics = result$diagnostics),
        file.path(out, "result_contract.rds"))
if (mode != "old_manual") {
  write.table(selection$audit, file.path(out, "lr_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(selection$missing_genes, file.path(out, "missing_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
row <- data.frame(case = case, cells = cells, mode = mode, version = as.character(packageVersion("SpatialESS")),
                  genes = length(genes), groups = nlevels(group), lr = nrow(manual_lr), nboot = parameters$nboot,
                  seed = parameters$seed.use, inference_and_filter_seconds = seconds,
                  graph_edges = length(result$graph$neighbors), active = nrow(result$records),
                  significant = sum(result$records$pvalue <= 0.05), status = "completed")
write.table(row, file.path(out, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(row)
