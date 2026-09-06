#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(SpatialESS))

set.seed(20260731)
cell_counts <- c(100L, 1000L, 10000L, 100000L, 1000000L)
zero_fractions <- c(0, 0.25, 0.50, 0.75, 0.90, 0.99)
rows <- list()
index <- 1L
for (cell_count in cell_counts) {
  for (zero_fraction in zero_fractions) {
    positive_count <- max(1L, round(cell_count * (1 - zero_fraction)))
    positive <- rexp(positive_count)
    iterations <- max(3L, min(10000L, round(2e6 / positive_count)))
    result <- SpatialESS:::benchmark_sparse_type7_cpp(
      positive, cell_count, iterations
    )
    rows[[index]] <- data.frame(
      cells = cell_count,
      positive_values = positive_count,
      zero_fraction = zero_fraction,
      iterations = iterations,
      sort_seconds = result$sort_seconds,
      select_seconds = result$select_seconds,
      speedup_select_over_sort = result$sort_seconds / result$select_seconds,
      absolute_difference = abs(result$sort_value - result$select_value)
    )
    index <- index + 1L
  }
}
output <- do.call(rbind, rows)
out_dir <- file.path(project_dir, "benchmarks", "results", "sparse_type7_kernel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(output, file.path(out_dir, "sort_vs_multiselect.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(output)

