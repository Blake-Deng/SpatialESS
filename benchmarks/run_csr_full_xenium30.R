#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SpatialESS))

input <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/Xenium_Prime_Human_Ovary_FF_cells.csv.gz"
out_dir <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/csr_full_xenium30"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading full Xenium coordinates...")
cells <- read.csv(
  gzfile(input),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
coords <- as.matrix(cells[, c("x_centroid", "y_centroid")])
rownames(coords) <- cells$cell_id
rm(cells)
gc()

message("Building 30-um compact CSR graph...")
build_start <- proc.time()[["elapsed"]]
graph <- build_radius_graph_csr(
  coords, radius = 30, weight = "gaussian", scale = 15,
  store_distance = FALSE, max_edges = 2e8
)
build_seconds <- proc.time()[["elapsed"]] - build_start
rm(coords)
gc()

message("Validating full CSR graph...")
validation_start <- proc.time()[["elapsed"]]
checks <- validate_spatial_csr_graph(graph)
validation_seconds <- proc.time()[["elapsed"]] - validation_start
stopifnot(isTRUE(attr(checks, "valid")))

degree <- graph$degree
summary <- data.frame(
  cells = length(graph$cells$cell_id),
  directed_edges = length(graph$neighbors),
  radius_um = graph$parameters$radius,
  weight = graph$parameters$weight,
  mean_degree = mean(degree),
  median_degree = median(degree),
  max_degree = max(degree),
  isolated_cells = sum(degree == 0),
  graph_object_mib = as.numeric(object.size(graph)) / 1024^2,
  build_seconds = build_seconds,
  validation_seconds = validation_seconds,
  all_invariants_pass = isTRUE(attr(checks, "valid"))
)
write.table(summary, file.path(out_dir, "csr_full_xenium30_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(check = names(checks), pass = as.logical(checks)),
            file.path(out_dir, "csr_full_xenium30_checks.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saving compact graph...")
saveRDS(graph, file.path(out_dir, "xenium_full_radius30_csr.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
print(checks)
