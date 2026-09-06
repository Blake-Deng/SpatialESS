#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SpatialCellChat)
})

root <- "/home/dzf/cellchat_acceleration/SpatialESS"
prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
coords_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/Xenium_Prime_Human_Ovary_FF_cells.csv.gz"
out_root <- file.path(root, "benchmarks", "results", "xenium_spatial_new_20260731")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

n_cells <- as.integer(Sys.getenv("XENIUM_SPATIAL_NCELLS", "5000"))
nperm <- as.integer(Sys.getenv("XENIUM_SPATIAL_NPERM", "20"))
panel <- Sys.getenv("XENIUM_SPATIAL_PANEL", "full")
seed <- as.integer(Sys.getenv("XENIUM_SPATIAL_SEED", "20260731"))
scale_distance <- as.numeric(Sys.getenv("XENIUM_SPATIAL_SCALE_DISTANCE", "20"))
radius <- as.numeric(Sys.getenv("XENIUM_SPATIAL_RADIUS_UM", "30"))
tol <- as.numeric(Sys.getenv("XENIUM_SPATIAL_TOL_UM", "5"))
if (!panel %in% c("small", "full")) stop("Invalid XENIUM_SPATIAL_PANEL")
if (n_cells < 20L || nperm < 1L || radius <= 0 || tol < 0) {
  stop("Invalid XENIUM_SPATIAL_NCELLS, NPERM, RADIUS_UM or TOL_UM")
}

object <- readRDS(prepared_path)
expression_all <- object@data.signaling
cell_id_all <- colnames(expression_all)
group_all <- object@idents
names(group_all) <- cell_id_all
coords_all <- read.csv(
  gzfile(coords_path),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
matched <- match(cell_id_all, coords_all$cell_id)
stopifnot(!anyNA(matched))
coords_all <- as.matrix(coords_all[matched, c("x_centroid", "y_centroid")])
rownames(coords_all) <- cell_id_all

center <- colMeans(coords_all)
selection_order <- order((coords_all[, 1] - center[1])^2 +
                         (coords_all[, 2] - center[2])^2)
if (n_cells >= ncol(expression_all)) {
  idx <- seq_len(ncol(expression_all))
  run_label <- "full"
} else {
  idx <- selection_order[seq_len(n_cells)]
  run_label <- sprintf("n%06d", n_cells)
}
run_label <- paste0(run_label, "_", panel)
expression <- expression_all[, idx, drop = FALSE]
coords <- coords_all[idx, , drop = FALSE]
group <- droplevels(group_all[idx])
names(group) <- colnames(expression)
lr <- object@LR$LRsig
if (panel == "small") {
  pathway <- if ("pathway_name" %in% colnames(lr)) as.character(lr$pathway_name) else rep("", nrow(lr))
  focus <- which(pathway %in% c("CXCL", "TNF", "IL17", "TGFb", "MIF"))
  keep <- c(focus, setdiff(seq_len(nrow(lr)), focus))
  lr <- lr[keep[seq_len(min(50L, length(keep)))], , drop = FALSE]
}

out_dir <- file.path(out_root, paste0(run_label, "_spatialcellchat_v3"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workflow_start <- proc.time()[["elapsed"]]
meta <- data.frame(labels = group, row.names = colnames(expression))

prep_start <- proc.time()[["elapsed"]]
chat <- createSpatialCellChat(
  object = expression, meta = meta, group.by = "labels",
  datatype = "spatial", coordinates = coords,
  spatial.factors = list(ratio = 1, tol = tol), do.sparse = TRUE
)
chat@DB <- object@DB
chat <- subsetData(chat)
chat@LR$LRsig <- lr
preparation_seconds <- proc.time()[["elapsed"]] - prep_start

prob_start <- proc.time()[["elapsed"]]
chat <- computeCommunProb(
  chat, LR.use = lr, raw.use = TRUE, Kh = 0.5, n = 1,
  distance.use = TRUE, interaction.range = radius, scale.distance = scale_distance,
  contact.dependent = TRUE, contact.range = 10
)
probability_seconds <- proc.time()[["elapsed"]] - prob_start

perm_start <- proc.time()[["elapsed"]]
chat <- computeAvgCommunProb(
  chat, avg.type = "avg", min.percent = 0, min.cells.sr = 1,
  do.permutation = TRUE, nboot = nperm, seed.use = seed,
  colocalization.use = FALSE
)
permutation_seconds <- proc.time()[["elapsed"]] - perm_start

prob <- chat@net$prob
pval <- chat@net$pval
active <- which(prob > 0, arr.ind = TRUE)
if (nrow(active)) {
  records <- data.frame(
    sender_group_name = dimnames(prob)[[1L]][active[, 1L]],
    receiver_group_name = dimnames(prob)[[2L]][active[, 2L]],
    interaction_name = dimnames(prob)[[3L]][active[, 3L]],
    probability = prob[active], pvalue = pval[active],
    stringsAsFactors = FALSE
  )
  records$pathway_name <- as.character(lr$pathway_name[
    match(records$interaction_name, lr$interaction_name)
  ])
} else {
  records <- data.frame(
    sender_group_name = character(), receiver_group_name = character(),
    interaction_name = character(), probability = numeric(),
    pvalue = numeric(), pathway_name = character(),
    stringsAsFactors = FALSE
  )
}
saveRDS(records, file.path(out_dir, "records.rds"), compress = FALSE)
write.table(records, gzfile(file.path(out_dir, "records.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(
  method = "SpatialCellChat_v3",
  dataset = "Xenium Prime fresh-frozen human ovary",
  cells = ncol(expression), genes = nrow(expression),
  groups = nlevels(group), lr = nrow(lr), nperm = nperm,
  radius_um = radius, tol_um = tol,
  active_group_lr_records = nrow(records),
  significant_p_0_05 = sum(records$pvalue <= 0.05),
  preparation_seconds = preparation_seconds,
  probability_seconds = probability_seconds,
  permutation_seconds = permutation_seconds,
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start
)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
