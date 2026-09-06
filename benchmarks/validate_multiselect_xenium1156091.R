#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(CellChat)
  library(SpatialESS)
})

prepared_path <- "/data/dzf/Xenium_Prime_Human_Ovary_FF/benchmark/prepared_xenium_ovary_full_grid50.rds"
graph_path <- file.path(
  project_dir, "benchmarks", "results", "csr_prepared1156091_xenium30",
  "xenium_prepared1156091_radius30_csr.rds"
)
reference_dir <- file.path(
  project_dir, "benchmarks", "results",
  "cellchat_compat_all627_prepared1156091"
)
reference_path <- file.path(
  reference_dir, "group_trimean_534genes_2357groups.rds"
)
components_path <- file.path(reference_dir, "lr_components_all627.rds")
out_dir <- file.path(
  project_dir, "benchmarks", "results", "multiselect_xenium1156091"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading frozen Xenium expression, graph, LR components and sort reference...")
workflow_start <- proc.time()[["elapsed"]]
object <- readRDS(prepared_path)
graph <- readRDS(graph_path)
components <- readRDS(components_path)
reference <- readRDS(reference_path)
expression <- object@data.signaling
stopifnot(
  identical(colnames(expression), graph$cells$cell_id),
  ncol(expression) == 1156091L,
  length(graph$cells$group_levels) == 2357L,
  length(components$genes) == 534L,
  identical(rownames(reference$average), components$genes)
)

message("Recomputing 534-gene x 2,357-group exact triMeans with multi-select...")
kernel_start <- proc.time()[["elapsed"]]
candidate <- summarize_cellchat_trimean(
  expression, graph$cells$group_code, components, normalize = TRUE
)
kernel_seconds <- proc.time()[["elapsed"]] - kernel_start
dimnames(candidate$average) <- dimnames(reference$average)
candidate$group_levels <- reference$group_levels

difference <- abs(candidate$average - reference$average)
summary <- data.frame(
  cells = ncol(expression),
  signaling_genes = nrow(expression),
  referenced_genes = nrow(candidate$average),
  groups = ncol(candidate$average),
  values_compared = length(candidate$average),
  exact_identical = identical(candidate$average, reference$average),
  max_absolute_difference = max(difference),
  differing_values = sum(difference != 0),
  multiselect_kernel_seconds = kernel_seconds,
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start
)

saveRDS(candidate, file.path(out_dir, "group_trimean_multiselect.rds"),
        compress = FALSE)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
if (!isTRUE(summary$exact_identical) || summary$max_absolute_difference != 0) {
  stop("Multi-select triMean differs from the frozen sort reference.", call. = FALSE)
}

