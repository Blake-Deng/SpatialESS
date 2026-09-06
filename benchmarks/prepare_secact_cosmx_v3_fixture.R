#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))))
data_path <- Sys.getenv("SECACT_DATA", "/data/dzf/SecAct_2026/LIHC_CosMx_data.rda")
out_root <- Sys.getenv(
  "SECACT_V3_OUT",
  "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead"
)
fixture_path <- file.path(out_root, "secact_cosmx_v3_common_fixture.rds")
min_genes <- as.integer(Sys.getenv("SECACT_V3_MIN_GENES", "50"))
scale_factor <- as.numeric(Sys.getenv("SECACT_V3_LIBRARY_SCALE", "1000"))
interaction_range <- as.numeric(Sys.getenv("SECACT_V3_RADIUS_UM", "30"))
contact_range <- as.numeric(Sys.getenv("SECACT_V3_CONTACT_UM", "10"))
tol <- as.numeric(Sys.getenv("SECACT_V3_TOL_UM", "5"))

if (!file.exists(data_path)) stop("Missing prepared HCC CosMx input: ", data_path)
if (is.na(min_genes) || min_genes < 0L || !is.finite(scale_factor) ||
    scale_factor <= 0 || !is.finite(interaction_range) ||
    interaction_range <= 0 || !is.finite(contact_range) ||
    contact_range <= 0 || !is.finite(tol) || tol < 0) {
  stop("Invalid fixture parameters.")
}
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(SpatialCellChat)
  library(SpatialESS)
})

data_env <- new.env(parent = emptyenv())
loaded <- load(data_path, envir = data_env)
required <- c("counts", "spotCoordinates", "metaData")
if (!all(required %in% loaded)) {
  stop("Prepared RDA lacks: ", paste(setdiff(required, loaded), collapse = ", "))
}
counts <- methods::as(data_env$counts, "dgCMatrix")
coordinates <- as.matrix(data_env$spotCoordinates)
metadata <- as.data.frame(data_env$metaData, stringsAsFactors = FALSE)
rm(data_env)

cell_id <- colnames(counts)
if (!all(cell_id %in% rownames(coordinates)) ||
    !all(cell_id %in% rownames(metadata))) {
  stop("Counts, coordinates and metadata do not contain the same cells.")
}
coordinates <- coordinates[cell_id, 1:2, drop = FALSE]
metadata <- metadata[cell_id, , drop = FALSE]
storage.mode(coordinates) <- "double"
detected_genes <- Matrix::colSums(counts != 0)
library_size <- Matrix::colSums(counts)
keep <- detected_genes >= min_genes & library_size > 0 &
  !is.na(metadata$cellType) & nzchar(as.character(metadata$cellType))
counts <- counts[, keep, drop = FALSE]
coordinates <- coordinates[keep, , drop = FALSE]
metadata <- metadata[keep, , drop = FALSE]
library_size <- library_size[keep]

expression <- counts %*% Matrix::Diagonal(x = scale_factor / library_size)
dimnames(expression) <- dimnames(counts)
expression@x <- log1p(expression@x)
expression <- Matrix::drop0(methods::as(expression, "dgCMatrix"))
group <- droplevels(factor(as.character(metadata$cellType)))
names(group) <- colnames(expression)
rm(counts, metadata, library_size)
gc()

data(CellChatDB.human, package = "SpatialCellChat")
db <- CellChatDB.human
meta <- data.frame(labels = group, row.names = colnames(expression))
tmp <- createSpatialCellChat(
  object = expression, meta = meta, group.by = "labels",
  datatype = "spatial", coordinates = coordinates,
  spatial.factors = list(ratio = 1, tol = tol), do.sparse = TRUE
)
tmp@DB <- db
tmp <- subsetData(tmp)
expression_signaling <- methods::as(tmp@data.signaling, "dgCMatrix")
rm(tmp, expression)
gc()

component_genes <- function(value) {
  value <- as.character(value)
  if (!is.na(value) && value %in% rownames(db$complex)) {
    columns <- grep("^subunit", colnames(db$complex), value = TRUE)
    value <- unlist(db$complex[value, columns, drop = FALSE], use.names = FALSE)
  }
  value <- as.character(value)
  value[!is.na(value) & nzchar(value)]
}
lr <- db$interaction
valid <- vapply(seq_len(nrow(lr)), function(i) {
  required_genes <- c(
    component_genes(lr$ligand[[i]]),
    component_genes(lr$receptor[[i]])
  )
  length(required_genes) > 0L &&
    all(required_genes %in% rownames(expression_signaling))
}, logical(1))
lr <- lr[valid, , drop = FALSE]
lr <- lr[!duplicated(lr$interaction_name), , drop = FALSE]
rownames(lr) <- lr$interaction_name
db$interaction <- lr

priority_pathways <- c("CXCL", "TNF", "IL17", "TGFb", "MIF", "VEGF")
priority <- which(lr$pathway_name %in% priority_pathways)
fill <- setdiff(seq_len(nrow(lr)), priority)
small_index <- head(unique(c(priority, fill)), 50L)
lr_small <- lr[small_index, , drop = FALSE]

center <- apply(coordinates, 2L, stats::median)
distance_to_center <- sqrt(rowSums(sweep(coordinates, 2L, center)^2))
selection_order <- order(distance_to_center, rownames(coordinates), method = "radix")

message("Building one CSR graph to calibrate scale.distance...")
calibration_graph <- build_radius_graph_csr(
  coordinates, radius = interaction_range + tol,
  weight = "binary", store_distance = TRUE, max_edges = 5e8
)
positive_distance <- calibration_graph$distances[calibration_graph$distances > 0]
if (!length(positive_distance)) stop("No positive spatial distances inside the interaction range.")
minimum_nonzero_distance <- min(positive_distance)
scale_distance <- max(20, 1.5 / minimum_nonzero_distance)
calibration_edges <- length(calibration_graph$neighbors)
rm(calibration_graph, positive_distance)
gc()

fixture <- list(
  expression = expression_signaling,
  coordinates_um = coordinates,
  group = group,
  selection_order = selection_order,
  db = db,
  lr_small = lr_small,
  lr_full = lr,
  spatial_factors = list(ratio = 1, tol = tol),
  parameters = list(
    dataset = "Bruker/NanoString CosMx SMI CancerousLiver",
    source_article_doi = "10.1038/s41592-026-03172-0",
    source_data = normalizePath(data_path),
    selection = "nested Euclidean distance from median tissue coordinate",
    normalization = sprintf("log1p(count / cell library size * %g)", scale_factor),
    min_genes = min_genes,
    interaction_range_um = interaction_range,
    contact_range_um = contact_range,
    tol_um = tol,
    minimum_nonzero_distance_um = minimum_nonzero_distance,
    scale_distance = scale_distance,
    full_calibration_graph_edges = calibration_edges,
    seed = 20260803L,
    priority_pathways = priority_pathways
  )
)
saveRDS(fixture, fixture_path, compress = FALSE)

sizes <- c(250L, 1000L, 5000L, 10000L, 25000L, 50000L, 100000L, ncol(expression_signaling))
sizes <- unique(sizes[sizes <= ncol(expression_signaling)])
group_counts <- do.call(rbind, lapply(sizes, function(n) {
  idx <- selection_order[seq_len(n)]
  tab <- table(droplevels(group[idx]))
  data.frame(cells = n, group = names(tab), count = as.integer(tab))
}))
write.table(group_counts, file.path(out_root, "nested_group_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(
  cells = ncol(expression_signaling),
  signaling_genes = nrow(expression_signaling),
  groups = nlevels(group),
  valid_lr_full = nrow(lr),
  valid_lr_small = nrow(lr_small),
  minimum_nonzero_distance_um = minimum_nonzero_distance,
  scale_distance = scale_distance,
  interaction_range_um = interaction_range,
  tol_um = tol,
  contact_range_um = contact_range,
  calibration_graph_edges = calibration_edges
)
write.table(summary, file.path(out_root, "fixture_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_root, "prepare_sessionInfo.txt"))
print(summary)

