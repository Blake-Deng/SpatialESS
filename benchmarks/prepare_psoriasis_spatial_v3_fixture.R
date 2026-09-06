#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))))
reference_dir <- Sys.getenv(
  "SPATIALCELLCHAT_V3_SOURCE",
  "/home/dzf/cellchat_acceleration/SpatialESS_research_20260730/references/SpatialCellChat"
)
out_dir <- Sys.getenv(
  "SPATIALESS_PSORIASIS_OUT",
  file.path(project_dir, "benchmarks", "results", "spatialcellchat_v3_psoriasis")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SpatialCellChat)
})

data_file <- file.path(reference_dir, "tutorial", "visium_human_psoriasis.RData")
if (!file.exists(data_file)) stop("Missing official psoriasis fixture: ", data_file)
load(data_file)
if (!exists("seu") || !inherits(seu, "Seurat")) stop("Official RData lacks Seurat object 'seu'.")

DefaultAssay(seu) <- "Spatial"
seu@assays$Spatial@counts <- seu@assays$Spatial@data
seu <- NormalizeData(seu, verbose = FALSE)
expression <- GetAssayData(seu, assay = "Spatial", layer = "data")
expression <- methods::as(expression, "dgCMatrix")

decomposition <- t(as.matrix(seu@assays[["predictions"]]@data))
decomposition[decomposition < 0.1] <- 0

coordinates <- GetTissueCoordinates(
  seu, scale = NULL, cols = c("imagerow", "imagecol")
)
coordinates <- as.matrix(coordinates)
coordinates[, 2L] <- max(coordinates[, 2L]) - coordinates[, 2L]
tmp <- coordinates
coordinates[, 1L] <- tmp[, 2L]
coordinates[, 2L] <- tmp[, 1L]
coordinates[, 1L] <- max(coordinates[, 1L]) - coordinates[, 1L]
colnames(coordinates) <- c("x", "y")

scale_json <- file.path(
  reference_dir, "tutorial", "spatial_imaging_data_visium_psoriasisSkin",
  "scalefactors_json.json"
)
scale_factors <- jsonlite::fromJSON(scale_json)
spot_size_um <- 65
ratio <- spot_size_um / scale_factors$spot_diameter_fullres
coordinates_um <- coordinates * ratio

cell_id <- colnames(expression)
if (!identical(cell_id, rownames(coordinates))) {
  coordinates <- coordinates[cell_id, , drop = FALSE]
  coordinates_um <- coordinates_um[cell_id, , drop = FALSE]
}
if (!identical(cell_id, rownames(decomposition))) {
  decomposition <- decomposition[cell_id, , drop = FALSE]
}
group <- droplevels(factor(seu$final.clusters[cell_id]))
names(group) <- cell_id

data(CellChatDB.human, package = "SpatialCellChat")
db <- CellChatDB.human
lr <- db$interaction[
  db$interaction$annotation %in%
    c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
  , drop = FALSE
]

component_genes <- function(value) {
  value <- as.character(value)
  if (!is.na(value) && value %in% rownames(db$complex)) {
    cols <- grep("^subunit", colnames(db$complex))
    value <- unlist(db$complex[value, cols, drop = FALSE], use.names = FALSE)
  }
  value <- as.character(value)
  value[!is.na(value) & nzchar(value)]
}
valid <- vapply(seq_len(nrow(lr)), function(i) {
  all(c(component_genes(lr$ligand[[i]]), component_genes(lr$receptor[[i]])) %in%
        rownames(expression))
}, logical(1))
lr <- lr[valid, , drop = FALSE]
lr <- lr[!duplicated(lr$interaction_name), , drop = FALSE]
rownames(lr) <- lr$interaction_name
db$interaction <- lr

priority_pathways <- c("TNF", "CXCL", "IL17")
priority <- which(lr$pathway_name %in% priority_pathways)
fill <- setdiff(seq_len(nrow(lr)), priority)
small_index <- unique(c(priority, head(fill, max(0L, 50L - length(priority)))))
small_index <- head(small_index, 60L)
lr_small <- lr[small_index, , drop = FALSE]

center <- apply(coordinates_um, 2L, median)
distance_to_center <- sqrt(rowSums(sweep(coordinates_um, 2L, center)^2))
selection_order <- order(distance_to_center, rownames(coordinates_um), method = "radix")

fixture <- list(
  expression = expression,
  coordinates_pixel = coordinates,
  coordinates_um = coordinates_um,
  group = group,
  decomposition = decomposition,
  selection_order = selection_order,
  db = db,
  lr_small = lr_small,
  lr_full = lr,
  spatial_factors = list(ratio = ratio, tol = spot_size_um / 2),
  parameters = list(
    dataset = "SpatialCellChat official 10X Visium human psoriasis skin",
    source_file = normalizePath(data_file),
    source_commit = "be88f300464cf970b36a6165ea7d2e41a47a7511",
    normalization = "official tutorial: Spatial data copied to counts then Seurat NormalizeData",
    nesting = "ascending Euclidean distance from median tissue coordinate",
    interaction_range_um = 250,
    contact_range_um = 10,
    scale_distance = 0.2,
    priority_pathways = priority_pathways
  )
)
saveRDS(fixture, file.path(out_dir, "psoriasis_common_fixture.rds"), compress = FALSE)

sizes <- c(100L, 250L, 500L, ncol(expression))
sizes <- unique(sizes[sizes <= ncol(expression)])
group_counts <- do.call(rbind, lapply(sizes, function(n) {
  idx <- selection_order[seq_len(n)]
  tab <- table(droplevels(group[idx]))
  data.frame(spots = n, group = names(tab), count = as.integer(tab))
}))
write.table(group_counts, file.path(out_dir, "nested_group_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(
  spots = ncol(expression), genes = nrow(expression), groups = nlevels(group),
  valid_lr_full = nrow(lr), valid_lr_small = nrow(lr_small),
  priority_lr_small = sum(lr_small$pathway_name %in% priority_pathways),
  ratio_um_per_pixel = ratio, tol_um = spot_size_um / 2
)
write.table(summary, file.path(out_dir, "fixture_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "prepare_sessionInfo.txt"))
print(summary)

