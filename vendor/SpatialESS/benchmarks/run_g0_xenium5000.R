#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialESS)
  library(data.table)
})

input <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/Xenium_Prime_Human_Ovary_FF_cells.csv.gz"
out_dir <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/g0_xenium5000"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- read.csv(gzfile(input), nrows = 5000L, check.names = FALSE)
cols <- c("cell_id", "x_centroid", "y_centroid")
cells <- cells[, cols, drop = FALSE]
#
coords <- as.matrix(cells[, c("x_centroid", "y_centroid")])
rownames(coords) <- cells$cell_id

run_one <- function(radius, weight) {
  gc()
  start <- proc.time()[["elapsed"]]
  graph <- build_radius_graph(coords, radius = radius, weight = weight,
                              scale = radius / 2)
  elapsed <- proc.time()[["elapsed"]] - start
  checks <- validate_spatial_graph(graph)
  stopifnot(isTRUE(attr(checks, "valid")))
  data.frame(
    cells = nrow(coords), radius_um = radius, weight = weight,
    directed_edges = nrow(graph$edges),
    mean_out_degree = nrow(graph$edges) / nrow(coords),
    isolated_cells = sum(tabulate(graph$edges$sender, nbins = nrow(coords)) == 0L),
    wall_seconds = elapsed,
    all_invariants_pass = isTRUE(attr(checks, "valid"))
  )
}

result <- do.call(rbind, list(
  run_one(15, "binary"),
  run_one(30, "gaussian"),
  run_one(60, "exponential")
))
write.table(result, file.path(out_dir, "g0_xenium5000.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(result)

