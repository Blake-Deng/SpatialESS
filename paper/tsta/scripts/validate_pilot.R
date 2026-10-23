args <- commandArgs(TRUE)
root <- normalizePath(args[1]); sid <- "TSTA_0001"
.libPaths(c(file.path(root, "lib"), .libPaths()))
suppressPackageStartupMessages({library(SpatialESS); library(data.table); library(jsonlite)})
cfg <- fromJSON(file.path(root, "config/design.json"))
out <- file.path(root, "results", sid)
result <- readRDS(file.path(out, "SpatialESS_result.rds"))
meta <- fread(cmd=paste("gzip -dc --", shQuote(file.path(out, "spot_metadata.tsv.gz"))))[retained == TRUE]
summary <- fread(file.path(out, "summary.tsv"))
xy <- cbind(meta$array_row * 0.5 * sqrt(3) * 100, meta$array_col * 0.5 * 100)
radius_sq <- cfg$parameters$interaction.range^2
epsilon <- max(1, radius_sq)*16*.Machine$double.eps
# Direct all-pair row checks independently verify the hashed CSR edge set.
graph <- result$graph
expected <- 0
for (i in seq_len(nrow(xy))) {
  d2 <- (xy[,1]-xy[i,1])^2 + (xy[,2]-xy[i,2])^2
  expected_neighbors <- which(seq_len(nrow(xy)) != i & d2 <= radius_sq+epsilon)
  begin <- graph$offsets[i]+1
  end <- graph$offsets[i+1L]
  actual <- if (end >= begin) graph$neighbors[seq.int(begin, end)] else integer()
  stopifnot(identical(sort(as.integer(actual)), as.integer(expected_neighbors)))
  expected <- expected+length(expected_neighbors)
}
stopifnot(expected == summary$E, length(result$graph$neighbors) == expected,
          result$parameters$nboot == 100L, result$parameters$seed.use == 1L,
          !anyNA(result$records$probability), !anyNA(result$records$pvalue),
          summary$unit == "Visium_spot", summary$contact_dependent == FALSE)
report <- list(status="PASS", pilot=sid, N=nrow(xy), E=expected,
               check="Every saved CSR neighbor set agrees with an independent direct all-pair distance oracle, including boundary tolerance.",
               nboot=result$parameters$nboot, active_records=nrow(result$records),
               checked=format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
write_json(report, file.path(root, "reports/pilot_validation.json"), pretty=TRUE, auto_unbox=TRUE)
cat("PILOT VALIDATION PASS", sid, "N", nrow(xy), "E", expected, "\n")
