#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
graph_path <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/csr_prepared1156091_xenium30/xenium_prepared1156091_radius30_csr.rds"
out_dir <- "/home/dzf/cellchat_acceleration/SpatialESS/benchmarks/results/csr_field_full_sema4c_plxnb2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading exact prepared object and compact graph...")
object <- readRDS(prepared_path)
graph <- readRDS(graph_path)
expression <- object@data.signaling
stopifnot(identical(colnames(expression), graph$cells$cell_id))
stopifnot(all(c("SEMA4C", "PLXNB2") %in% rownames(expression)))

expression_max <- max(expression@x)
source <- as.numeric(expression["SEMA4C", , drop = TRUE]) / expression_max
receptor <- as.numeric(expression["PLXNB2", , drop = TRUE]) / expression_max
names(source) <- names(receptor) <- graph$cells$cell_id

message("Solving full graph reaction-diffusion field...")
start <- proc.time()[["elapsed"]]
field <- solve_ligand_field_csr(
  graph, source = source, receptor = receptor,
  diffusion = 1, decay = 1, uptake = 0.25, production_rate = 1,
  tolerance = 1e-8, max_iterations = 2000
)
field_seconds <- proc.time()[["elapsed"]] - start
stopifnot(field$converged)
stopifnot(field$minimum_concentration >= -1e-8)

summary <- data.frame(
  cells = length(source),
  graph_edges = length(graph$neighbors),
  ligand = "SEMA4C", receptor = "PLXNB2",
  source_supported_cells = sum(source > 0),
  receptor_supported_cells = sum(receptor > 0),
  diffusion = 1, decay = 1, uptake = 0.25,
  tolerance = 1e-8,
  iterations = field$iterations,
  relative_residual = field$relative_residual,
  mass_balance_relative_error = field$mass_balance_relative_error,
  minimum_concentration = field$minimum_concentration,
  maximum_concentration = max(field$concentration),
  mean_concentration = mean(field$concentration),
  field_seconds = field_seconds,
  converged = field$converged
)
write.table(summary, file.path(out_dir, "csr_field_full_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(summary = summary, concentration = field$concentration,
             parameters = field$parameters),
        file.path(out_dir, "csr_field_full_result.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
