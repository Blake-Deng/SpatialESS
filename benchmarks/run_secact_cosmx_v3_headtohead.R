#!/usr/bin/env Rscript

root <- Sys.getenv(
  "SECACT_V3_OUT",
  "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead"
)
method <- Sys.getenv("SECACT_V3_METHOD", "spatialess")
cells <- as.integer(Sys.getenv("SECACT_V3_CELLS", "250"))
panel <- Sys.getenv("SECACT_V3_PANEL", "small")
nperm <- as.integer(Sys.getenv("SECACT_V3_NPERM", "20"))
seed <- as.integer(Sys.getenv("SECACT_V3_SEED", "20260803"))
if (!method %in% c("spatialess", "spatialcellchat_v3")) stop("Unknown method: ", method)
if (!panel %in% c("small", "full")) stop("Unknown panel: ", panel)
if (is.na(cells) || cells < 20L || is.na(nperm) || nperm < 1L || is.na(seed)) {
  stop("Invalid cells, nperm or seed.")
}

suppressPackageStartupMessages(library(Matrix))
fixture <- readRDS(file.path(root, "secact_cosmx_v3_common_fixture.rds"))
if (cells > ncol(fixture$expression)) stop("Requested cells exceed fixture.")
idx <- fixture$selection_order[seq_len(cells)]
expression <- fixture$expression[, idx, drop = FALSE]
coordinates <- fixture$coordinates_um[idx, , drop = FALSE]
group <- droplevels(fixture$group[idx])
names(group) <- colnames(expression)
lr <- fixture[[paste0("lr_", panel)]]
db <- fixture$db

out_dir <- file.path(
  root,
  sprintf("cells%06d_%s_%s_perm%05d", cells, panel, method, nperm)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workflow_start <- proc.time()[["elapsed"]]

if (method == "spatialcellchat_v3") {
  # The OS-level benchmark limit controls memory. The future package's
  # unrelated 500 MiB globals guard otherwise aborts large sequential runs
  # before their true memory requirement can be measured.
  options(future.globals.maxSize = Inf)
  suppressPackageStartupMessages(library(SpatialCellChat))
  meta <- data.frame(labels = group, row.names = colnames(expression))
  prep_start <- proc.time()[["elapsed"]]
  chat <- createSpatialCellChat(
    object = expression, meta = meta, group.by = "labels",
    datatype = "spatial", coordinates = coordinates,
    spatial.factors = fixture$spatial_factors, do.sparse = TRUE
  )
  chat@DB <- db
  chat <- subsetData(chat)
  chat@LR$LRsig <- lr
  preparation_seconds <- proc.time()[["elapsed"]] - prep_start

  individual_start <- proc.time()[["elapsed"]]
  chat <- computeCommunProb(
    chat, LR.use = lr, raw.use = TRUE, Kh = 0.5, n = 1,
    distance.use = TRUE,
    interaction.range = fixture$parameters$interaction_range_um,
    scale.distance = fixture$parameters$scale_distance,
    contact.dependent = TRUE,
    contact.range = fixture$parameters$contact_range_um
  )
  individual_seconds <- proc.time()[["elapsed"]] - individual_start

  group_start <- proc.time()[["elapsed"]]
  chat <- computeAvgCommunProb(
    chat, avg.type = "avg", min.percent = 0, min.cells.sr = 1,
    do.permutation = TRUE, nboot = nperm, seed.use = seed,
    colocalization.use = FALSE
  )
  group_permutation_seconds <- proc.time()[["elapsed"]] - group_start

  prob <- chat@net$prob
  pval <- chat@net$pval
  active <- which(prob > 0, arr.ind = TRUE)
  records <- data.frame(
    sender_group_name = dimnames(prob)[[1L]][active[, 1L]],
    receiver_group_name = dimnames(prob)[[2L]][active[, 2L]],
    interaction_name = dimnames(prob)[[3L]][active[, 3L]],
    probability = prob[active],
    pvalue = pval[active],
    active_edges = NA_real_,
    stringsAsFactors = FALSE
  )
  stored_links <- sum(vapply(chat@net$tmp$prob.cell, function(x) length(x@x), numeric(1)))
  graph_edges <- NA_real_
  minimum_distance <- min(chat@images$result.computeCellDistance$d.spatial@x)
  rm(chat)
  gc()
} else {
  suppressPackageStartupMessages(library(SpatialESS))
  preparation_seconds <- 0
  individual_start <- proc.time()[["elapsed"]]
  result <- spatialess(
    expression = expression,
    coordinates = coordinates,
    group = group,
    lr = lr,
    complex = db$complex,
    cofactor = db$cofactor,
    ratio = fixture$spatial_factors$ratio,
    tol = fixture$spatial_factors$tol,
    interaction.range = fixture$parameters$interaction_range_um,
    scale.distance = fixture$parameters$scale_distance,
    contact.range = fixture$parameters$contact_range_um,
    Kh = 0.5, n = 1, use.AGAN = TRUE,
    contact.dependent = TRUE, contact.dependent.forced = FALSE,
    min.percent = 0, min.cells.sr = 1,
    nboot = nperm, seed.use = seed,
    max_edges = 5e8, max_output_records = 1e8,
    max_active_edges_per_lr = 5e8
  )
  individual_seconds <- proc.time()[["elapsed"]] - individual_start
  group_permutation_seconds <- NA_real_
  records <- result$records
  graph_edges <- length(result$graph$neighbors)
  stored_links <- sum(result$diagnostics$active_edges_per_lr)
  minimum_distance <- result$diagnostics$minimum_nonzero_distance
}

match_index <- match(records$interaction_name, lr$interaction_name)
records$pathway_name <- as.character(lr$pathway_name[match_index])
records$annotation <- as.character(lr$annotation[match_index])
records$method <- method
records$cells <- cells
records$panel <- panel
records$nperm <- nperm
saveRDS(records, file.path(out_dir, "records.rds"), compress = FALSE)
write.table(records, gzfile(file.path(out_dir, "records.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary <- data.frame(
  method = method,
  dataset = fixture$parameters$dataset,
  cells = cells,
  signaling_genes = nrow(expression),
  groups = nlevels(group),
  lr = nrow(lr),
  panel = panel,
  nperm = nperm,
  radius_um = fixture$parameters$interaction_range_um,
  tol_um = fixture$parameters$tol_um,
  contact_um = fixture$parameters$contact_range_um,
  scale_distance = fixture$parameters$scale_distance,
  minimum_nonzero_distance_um = minimum_distance,
  preparation_seconds = preparation_seconds,
  individual_or_stream_seconds = individual_seconds,
  group_permutation_seconds = group_permutation_seconds,
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start,
  graph_edges = graph_edges,
  stored_active_links = stored_links,
  active_group_lr_records = nrow(records),
  significant_p_0_05 = sum(records$pvalue <= 0.05)
)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
