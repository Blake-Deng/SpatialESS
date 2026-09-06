#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
coords_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/Xenium_Prime_Human_Ovary_FF_cells.csv.gz"
out_dir <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/csr_prepared1156091_xenium30"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading the exact prepared 50-grid CellChat object...")
object <- readRDS(prepared_path)
cell_id <- colnames(object@data.signaling)
group <- object@idents
stopifnot(length(cell_id) == length(group))
stopifnot(!anyDuplicated(cell_id))

message("Reading and matching full Xenium coordinates...")
cells <- read.csv(
  gzfile(coords_path),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
matched <- match(cell_id, cells$cell_id)
stopifnot(!anyNA(matched))
coords <- as.matrix(cells[matched, c("x_centroid", "y_centroid")])
rownames(coords) <- cell_id
rm(cells)
gc()

message("Building exact prepared-cell 30-um CSR graph...")
build_start <- proc.time()[["elapsed"]]
graph <- build_radius_graph_csr(
  coords, radius = 30, weight = "gaussian", scale = 15,
  store_distance = FALSE, max_edges = 2e8
)
build_seconds <- proc.time()[["elapsed"]] - build_start
graph$cells$group_code <- as.integer(group)
graph$cells$group_levels <- levels(group)
rm(coords, group)
gc()

message("Validating exact prepared-cell CSR graph...")
validation_start <- proc.time()[["elapsed"]]
checks <- validate_spatial_csr_graph(graph)
validation_seconds <- proc.time()[["elapsed"]] - validation_start
stopifnot(isTRUE(attr(checks, "valid")))

summary <- data.frame(
  cells = length(graph$cells$cell_id),
  signaling_genes = nrow(object@data.signaling),
  lr_records = nrow(object@LR$LRsig),
  groups = length(graph$cells$group_levels),
  directed_edges = length(graph$neighbors),
  radius_um = graph$parameters$radius,
  weight = graph$parameters$weight,
  mean_degree = mean(graph$degree),
  median_degree = median(graph$degree),
  max_degree = max(graph$degree),
  isolated_cells = sum(graph$degree == 0),
  graph_object_mib = as.numeric(object.size(graph)) / 1024^2,
  build_seconds = build_seconds,
  validation_seconds = validation_seconds,
  all_invariants_pass = isTRUE(attr(checks, "valid"))
)
write.table(summary, file.path(out_dir, "csr_prepared1156091_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(check = names(checks), pass = as.logical(checks)),
            file.path(out_dir, "csr_prepared1156091_checks.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saving exact prepared-cell compact graph...")
saveRDS(graph, file.path(out_dir, "xenium_prepared1156091_radius30_csr.rds"),
        compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
print(checks)
