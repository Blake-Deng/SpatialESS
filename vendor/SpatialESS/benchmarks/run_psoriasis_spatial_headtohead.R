#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]))))
root <- Sys.getenv(
  "SPATIALESS_PSORIASIS_OUT",
  file.path(project_dir, "benchmarks", "results", "spatialcellchat_v3_psoriasis")
)
method <- Sys.getenv("SPATIAL_BENCH_METHOD", "spatialess")
spots <- as.integer(Sys.getenv("SPATIAL_BENCH_SPOTS", "100"))
panel <- Sys.getenv("SPATIAL_BENCH_PANEL", "small")
nperm <- as.integer(Sys.getenv("SPATIAL_BENCH_NPERM", "20"))
seed <- as.integer(Sys.getenv("SPATIAL_BENCH_SEED", "20260731"))
if (!method %in% c("spatialess", "spatialcellchat_v3")) stop("Unknown method: ", method)
if (!panel %in% c("small", "full")) stop("Unknown panel: ", panel)
if (is.na(spots) || spots < 20L || is.na(nperm) || nperm < 1L) stop("Invalid spots or nperm.")

suppressPackageStartupMessages(library(Matrix))
fixture <- readRDS(file.path(root, "psoriasis_common_fixture.rds"))
if (spots > ncol(fixture$expression)) stop("Requested spots exceed fixture.")
idx <- fixture$selection_order[seq_len(spots)]
expression <- fixture$expression[, idx, drop = FALSE]
coords_pixel <- fixture$coordinates_pixel[idx, , drop = FALSE]
coords_um <- fixture$coordinates_um[idx, , drop = FALSE]
group <- droplevels(fixture$group[idx])
names(group) <- colnames(expression)
lr <- fixture[[paste0("lr_", panel)]]
db <- fixture$db
db$interaction <- lr

out_dir <- file.path(root, sprintf("spots%04d_%s_%s_perm%05d", spots, panel, method, nperm))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
workflow_start <- proc.time()[["elapsed"]]

if (method == "spatialcellchat_v3") {
  suppressPackageStartupMessages(library(SpatialCellChat))
  meta <- data.frame(labels = group, row.names = colnames(expression))
  prep_start <- proc.time()[["elapsed"]]
  chat <- createSpatialCellChat(
    object = expression, meta = meta, group.by = "labels",
    datatype = "spatial", coordinates = coords_pixel,
    spatial.factors = fixture$spatial_factors, do.sparse = TRUE
  )
  chat@DB <- db
  chat <- subsetData(chat)
  chat@LR$LRsig <- lr
  preparation_seconds <- proc.time()[["elapsed"]] - prep_start

  individual_start <- proc.time()[["elapsed"]]
  chat <- computeCommunProb(
    chat, LR.use = lr, raw.use = TRUE, Kh = 0.5, n = 1,
    distance.use = TRUE,
    interaction.range = fixture$parameters$interaction_range_um,
    scale.distance = fixture$parameters$scale_distance,
    contact.dependent = TRUE,
    contact.range = fixture$parameters$contact_range_um
  )
  individual_seconds <- proc.time()[["elapsed"]] - individual_start

  group_start <- proc.time()[["elapsed"]]
  chat <- computeAvgCommunProb(
    chat, avg.type = "avg", min.percent = 0, min.cells.sr = 1,
    do.permutation = TRUE, nboot = nperm, seed.use = seed,
    colocalization.use = FALSE
  )
  group_permutation_seconds <- proc.time()[["elapsed"]] - group_start

  prob <- chat@net$prob
  pval <- chat@net$pval
  active <- which(prob > 0, arr.ind = TRUE)
  records <- data.frame(
    sender_group_name = dimnames(prob)[[1L]][active[, 1L]],
    receiver_group_name = dimnames(prob)[[2L]][active[, 2L]],
    interaction_name = dimnames(prob)[[3L]][active[, 3L]],
    probability = prob[active], pvalue = pval[active],
    stringsAsFactors = FALSE
  )
  link_count <- sum(vapply(chat@net$tmp$prob.cell, function(x) length(x@x), numeric(1)))
  graph_edges <- NA_real_
  exchangeable_fraction <- NA_real_
} else {
  suppressPackageStartupMessages(library(SpatialESS))
  preparation_start <- proc.time()[["elapsed"]]
  blocks <- build_spatial_blocks(
    coords_um, block_size = 500, sample_id = rep("psoriasis", spots),
    compartment = rep("tissue", spots)
  )
  preparation_seconds <- proc.time()[["elapsed"]] - preparation_start
  mechanism <- ifelse(lr$annotation == "Cell-Cell Contact", "contact", "diffusible")
  records_by_mechanism <- list()
  graph_edges <- 0
  individual_seconds <- 0
  group_permutation_seconds <- 0
  exchangeable <- numeric()
  for (mechanism_name in unique(mechanism)) {
    lr_part <- lr[mechanism == mechanism_name, , drop = FALSE]
    if (!nrow(lr_part)) next
    radius <- if (mechanism_name == "contact") {
      fixture$parameters$contact_range_um + fixture$spatial_factors$tol
    } else {
      fixture$parameters$interaction_range_um + fixture$spatial_factors$tol
    }
    part_start <- proc.time()[["elapsed"]]
    components <- prepare_cellchat_lr_components(
      lr_part, rownames(expression), db$complex, db$cofactor
    )
    graph <- build_radius_graph_csr(
      coords_um, radius = radius, sample_id = rep("psoriasis", spots),
      weight = "binary", store_distance = FALSE
    )
    support <- build_group_support_csr(graph, group)
    prepared <- prepare_sparse_trimean(expression, genes = components$genes, normalize = TRUE)
    individual_seconds <- individual_seconds + proc.time()[["elapsed"]] - part_start
    graph_edges <- graph_edges + length(graph$neighbors)

    perm_start <- proc.time()[["elapsed"]]
    result <- permutation_cellchat_group_support(
      prepared, group, support, components, blocks,
      nperm = nperm, seed = seed,
      finite_correction = FALSE, tail = "strict_greater", verbose = FALSE
    )
    group_permutation_seconds <-
      group_permutation_seconds + proc.time()[["elapsed"]] - perm_start
    part <- result$group_pairs
    if (nrow(part)) part$mechanism_graph <- mechanism_name
    records_by_mechanism[[mechanism_name]] <- part
    exchangeable <- c(exchangeable, result$diagnostics$exchangeable_cell_fraction)
  }
  records <- do.call(rbind, records_by_mechanism)
  rownames(records) <- NULL
  link_count <- graph_edges
  exchangeable_fraction <- min(exchangeable)
}

match_index <- match(records$interaction_name, lr$interaction_name)
records$pathway_name <- as.character(lr$pathway_name[match_index])
records$annotation <- as.character(lr$annotation[match_index])
records$method <- method
records$spots <- spots
records$panel <- panel
records$nperm <- nperm
saveRDS(records, file.path(out_dir, "records.rds"), compress = FALSE)
write.table(records, file.path(out_dir, "records.tsv.gz"), sep = "\t",
            quote = FALSE, row.names = FALSE)

summary <- data.frame(
  method = method, spots = spots, genes = nrow(expression),
  groups = nlevels(group), lr = nrow(lr), panel = panel, nperm = nperm,
  preparation_seconds = preparation_seconds,
  individual_or_graph_seconds = individual_seconds,
  group_permutation_seconds = group_permutation_seconds,
  end_to_end_seconds = proc.time()[["elapsed"]] - workflow_start,
  stored_links_or_edges = link_count,
  active_group_lr_records = nrow(records),
  significant_p_0_05 = sum(records$pvalue <= 0.05),
  exchangeable_cell_fraction = exchangeable_fraction
)
write.table(summary, file.path(out_dir, "summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)

