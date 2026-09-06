#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_n5000_grid4.rds"
coords_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/Xenium_Prime_Human_Ovary_FF_cells.csv.gz"
out_dir <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/f1_xenium4995"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(prepared_path)
expression <- object@data.signaling
cell_id <- colnames(expression)
coord_all <- read.csv(
  gzfile(coords_path),
  colClasses = c("character", "numeric", "numeric", rep("NULL", 11)),
  check.names = FALSE
)
matched <- match(cell_id, coord_all$cell_id)
stopifnot(!anyNA(matched))
coords <- as.matrix(coord_all[matched, c("x_centroid", "y_centroid")])
rownames(coords) <- cell_id
rm(coord_all)
gc()

graph <- build_radius_graph(coords, radius = 30, weight = "gaussian", scale = 15)
lr <- object@LR$LRsig
simple <- lr$ligand %in% rownames(expression) & lr$receptor %in% rownames(expression)
lr_simple <- lr[simple, , drop = FALSE]
support <- Matrix::rowSums(expression > 0)
priority <- support[lr_simple$ligand] * support[lr_simple$receptor]
lr_one <- lr_simple[which.max(priority), , drop = FALSE]

expression_max <- max(expression@x)
source <- as.numeric(expression[lr_one$ligand, , drop = TRUE]) / expression_max
receptor <- as.numeric(expression[lr_one$receptor, , drop = TRUE]) / expression_max
names(source) <- names(receptor) <- cell_id

start <- proc.time()[["elapsed"]]
field <- solve_ligand_field(
  graph, source = source, receptor = receptor,
  diffusion = 1, decay = 1, uptake = 0.25, production_rate = 1,
  tolerance = 1e-10, max_iterations = 2000
)
field_seconds <- proc.time()[["elapsed"]] - start
stopifnot(field$converged)
stopifnot(field$minimum_concentration >= -1e-10)

summary <- data.frame(
  cells = length(cell_id), graph_edges = nrow(graph$edges),
  ligand = lr_one$ligand, receptor = lr_one$receptor,
  interaction_name = lr_one$interaction_name,
  source_supported_cells = sum(source > 0),
  receptor_supported_cells = sum(receptor > 0),
  iterations = field$iterations, relative_residual = field$relative_residual,
  mass_balance_relative_error = field$mass_balance_relative_error,
  minimum_concentration = field$minimum_concentration,
  maximum_concentration = max(field$concentration),
  field_seconds = field_seconds, converged = field$converged
)
write.table(summary, file.path(out_dir, "f1_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
saveRDS(list(summary = summary, concentration = field$concentration,
             parameters = field$parameters),
        file.path(out_dir, "f1_field.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
