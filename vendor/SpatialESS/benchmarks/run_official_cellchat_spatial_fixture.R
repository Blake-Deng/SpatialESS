#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
out_dir <- file.path(project_dir, "benchmarks", "results", "official_spatial_fixture")

suppressPackageStartupMessages(library(CellChat))
fixture <- readRDS(file.path(out_dir, "fixture.rds"))
data(CellChatDB.human, package = "CellChat")

workflow_start <- proc.time()[["elapsed"]]
object <- createCellChat(
  object = fixture$expression,
  meta = fixture$meta,
  group.by = "labels",
  datatype = "spatial",
  coordinates = fixture$coordinates_pixel,
  scale.factors = fixture$scale_factors,
  do.sparse = TRUE
)
object@DB <- CellChatDB.human
object <- subsetData(object)
object@LR$LRsig <- fixture$lr
probability_start <- proc.time()[["elapsed"]]
object <- computeCommunProb(
  object, type = "triMean", raw.use = TRUE,
  distance.use = TRUE,
  interaction.length = fixture$interaction_length_um,
  scale.distance = fixture$scale_distance,
  k.min = fixture$k_min,
  nboot = fixture$nboot, seed.use = fixture$seed,
  Kh = 0.5, n = 1
)
probability_seconds <- proc.time()[["elapsed"]] - probability_start
object <- filterCommunication(object, min.cells = 10)
object <- computeCommunProbPathway(object)
object <- aggregateNet(object)
summary <- data.frame(
  method = "official_CellChat_spatial_1.6.1",
  cells = ncol(fixture$expression),
  groups = nlevels(fixture$group),
  lr = nrow(fixture$lr),
  probability_seconds = probability_seconds,
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start,
  active_probabilities = sum(object@net$prob > 0),
  significant_p_0_05 = sum(object@net$pval <= 0.05 & object@net$prob > 0)
)
saveRDS(object, file.path(out_dir, "official_cellchat_spatial.rds"), compress = FALSE)
write.table(summary, file.path(out_dir, "official_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "official_sessionInfo.txt"))
print(summary)

