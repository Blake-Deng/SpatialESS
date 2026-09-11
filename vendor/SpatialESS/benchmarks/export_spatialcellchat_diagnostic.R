# Source this in the collaborator session, then call:
# export_collaborator_debug(data_1k, chat_1k, result_1k_ess, "cervical_debug")
export_collaborator_debug <- function(data_1k, chat_1k, ess_result = NULL,
                                     out = "cervical_debug") {
  if (dir.exists(out)) stop("Choose a new output directory; existing files are not overwritten.")
  dir.create(out, recursive = TRUE)
  saveRDS(data_1k, file.path(out, "data_1k.rds"))
  saveRDS(chat_1k, file.path(out, "chat_1k.rds"))
  if (!is.null(ess_result)) saveRDS(ess_result, file.path(out, "ess_result.rds"))
  net <- chat_1k@net
  prob <- net[["prob", exact = TRUE]]
  pval <- net[["pval", exact = TRUE]]
  cell_prob <- net[["prob.cell", exact = TRUE]]
  cells <- net[["tmp"]][["prob.cell"]]
  counts <- list(
    net_names = names(net),
    database_lr = nrow(chat_1k@DB$interaction),
    selected_lr = nrow(chat_1k@LR$LRsig),
    cell_prob_dim = dim(cell_prob),
    group_prob_dim = dim(prob),
    group_pval_dim = dim(pval),
    group_stage_present = !is.null(prob) && !is.null(pval),
    active_cell_links = if (length(cells)) sum(vapply(cells, function(x) sum(x@x > 0), numeric(1))) else NA_real_,
    active_group_records = if (!is.null(prob)) sum(prob > 0, na.rm = TRUE) else NA_real_,
    significant_group_records = if (!is.null(prob) && !is.null(pval)) sum(prob > 0 & pval < .05, na.rm = TRUE) else NA_real_,
    official_parameters = chat_1k@options$parameter,
    ess_parameters = if (!is.null(ess_result)) ess_result$parameters else NULL
  )
  dput(counts, file.path(out, "stage_audit.R"))
  function_names <- c("computeCommunProb", "computeCommunProb1", "computeCellDistance",
    "computeCellDistance1", "computeAvgCommunProb", "computeExpr_LR",
    "computeExpr_coreceptor", "computeExpr_agonist", "computeExpr_antagonist",
    "createPspatialFrom_dspatial", "createCellCellContactMatrixFrom_dspatial")
  origins <- lapply(function_names, function(n) {
    f <- get0(n, envir = .GlobalEnv, inherits = TRUE)
    if (is.function(f)) {
      writeLines(c(paste0(n, " <-"), deparse(f)), file.path(out, paste0(n, ".R")))
    }
    data.frame(function_name = n, found = is.function(f),
      environment = if (is.function(f)) environmentName(environment(f)) else NA_character_,
      search_path = paste(find(n), collapse = ";"))
  })
  write.table(do.call(rbind, origins), file.path(out, "function_origins.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  print(counts)
  invisible(out)
}
