#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(SpatialESS)
})

env_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(sprintf("%s must be one positive integer.", name), call. = FALSE)
  }
  value
}

n_replicates <- env_integer("SPATIALESS_CALIBRATION_REPLICATES", 100L)
nperm <- env_integer("SPATIALESS_CALIBRATION_NPERM", 999L)
seed <- env_integer("SPATIALESS_CALIBRATION_SEED", 20260731L)
grid_side <- env_integer("SPATIALESS_CALIBRATION_GRID_SIDE", 20L)
if (grid_side %% 2L != 0L || grid_side < 10L) {
  stop("SPATIALESS_CALIBRATION_GRID_SIDE must be even and at least 10.", call. = FALSE)
}
out_dir <- Sys.getenv(
  "SPATIALESS_CALIBRATION_OUT",
  file.path(
    project_dir, "benchmarks", "results",
    sprintf("h1_semisimulated_rep%d_perm%d", n_replicates, nperm)
  )
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coordinates <- expand.grid(
  x = seq_len(grid_side) - 1L,
  y = seq_len(grid_side) - 1L
)
coordinates <- as.matrix(coordinates)
rownames(coordinates) <- paste0("cell", seq_len(nrow(coordinates)))
group <- factor(
  ifelse((coordinates[, "x"] + coordinates[, "y"]) %% 2L == 0L, "A", "B"),
  levels = c("A", "B")
)
sample_id <- rep("sample_1", nrow(coordinates))
compartment <- rep("tissue", nrow(coordinates))

graph <- build_radius_graph_csr(
  coordinates, radius = sqrt(2) + 1e-8,
  sample_id = sample_id, weight = "binary", store_distance = FALSE
)
graph$cells$group_code <- as.integer(group)
graph$cells$group_levels <- levels(group)
support <- build_group_support_csr(graph, group)
blocks <- build_spatial_blocks(
  coordinates, block_size = 5,
  sample_id = sample_id, compartment = compartment, origin = c(0, 0)
)

complex <- data.frame(subunit_1 = character(), row.names = character())
cofactor <- data.frame(cofactor1 = character(), row.names = character())
lr <- data.frame(interaction_name = "L_R", ligand = "L", receptor = "R")
components <- prepare_cellchat_lr_components(
  lr, c("L", "R"), complex, cofactor
)

simulate_expression <- function(replicate_seed, effect = 0) {
  set.seed(replicate_seed)
  block_baseline <- rgamma(length(blocks$block_levels), shape = 6, rate = 6)
  cell_depth <- rgamma(nrow(coordinates), shape = 10, rate = 10)
  baseline <- 3 * block_baseline[blocks$block_code] * cell_depth
  ligand_rate <- baseline + effect * (group == "A")
  receptor_rate <- baseline + effect * (group == "B")
  values <- rbind(
    L = rpois(nrow(coordinates), ligand_rate),
    R = rpois(nrow(coordinates), receptor_rate)
  )
  dropout <- matrix(
    runif(length(values)) < 0.15,
    nrow = nrow(values), ncol = ncol(values)
  )
  values[dropout] <- 0
  Matrix(
    values, sparse = TRUE,
    dimnames = list(c("L", "R"), rownames(coordinates))
  )
}

extract_pair_pvalues <- function(expression, permutation_seed) {
  prepared <- prepare_sparse_trimean(
    expression, genes = components$genes, normalize = FALSE
  )
  result <- permutation_cellchat_group_support(
    prepared, group, support, components, blocks,
    nperm = nperm, seed = permutation_seed,
    finite_correction = TRUE, tail = "greater_equal",
    retain_null_moments = FALSE, verbose = FALSE
  )
  keys <- paste(
    result$group_pairs$sender_group_name,
    result$group_pairs$receiver_group_name,
    sep = "->"
  )
  all_keys <- c("A->A", "A->B", "B->A", "B->B")
  pvalue <- rep(1, length(all_keys))
  names(pvalue) <- all_keys
  pvalue[keys] <- result$group_pairs$pvalue
  pvalue
}

run_one <- function(replicate_index) {
  data_seed <- seed + replicate_index * 2L
  permutation_seed <- seed + n_replicates * 2L +
    (replicate_index - 1L) * nperm * 2L
  if (as.double(permutation_seed) + 2 * nperm > .Machine$integer.max) {
    stop("Seed sequence exceeds the R integer range.", call. = FALSE)
  }
  null_p <- extract_pair_pvalues(
    simulate_expression(data_seed, effect = 0), permutation_seed
  )
  alternative_p <- extract_pair_pvalues(
    simulate_expression(data_seed + 1L, effect = 5),
    permutation_seed + nperm
  )
  null_q <- p.adjust(null_p, method = "BH")
  alternative_q <- p.adjust(alternative_p, method = "BH")
  alternative_discovery <- alternative_q <= 0.05
  false_discoveries <- sum(alternative_discovery[names(alternative_discovery) != "A->B"])
  discoveries <- sum(alternative_discovery)
  data.frame(
    replicate = replicate_index,
    null_target_p = null_p[["A->B"]],
    null_target_q = null_q[["A->B"]],
    null_any_discovery = any(null_q <= 0.05),
    alternative_target_p = alternative_p[["A->B"]],
    alternative_target_q = alternative_q[["A->B"]],
    alternative_target_discovery = alternative_discovery[["A->B"]],
    alternative_discoveries = discoveries,
    alternative_false_discoveries = false_discoveries,
    alternative_fdp = false_discoveries / max(1L, discoveries)
  )
}

message(sprintf(
  "Running %d null/alternative replicates with %d permutations each...",
  n_replicates, nperm
))
workflow_start <- proc.time()[["elapsed"]]
replicates <- lapply(seq_len(n_replicates), function(index) {
  if (index == 1L || index == n_replicates || index %% 10L == 0L) {
    message(sprintf("Calibration replicate %d/%d", index, n_replicates))
  }
  run_one(index)
})
replicates <- do.call(rbind, replicates)
wall_seconds <- proc.time()[["elapsed"]] - workflow_start

null_rejections <- sum(replicates$null_target_p <= 0.05)
type1_interval <- stats::binom.test(
  null_rejections, n_replicates, conf.level = 0.95
)$conf.int
summary <- data.frame(
  cells = nrow(coordinates),
  groups = nlevels(group),
  graph_edges = length(graph$neighbors),
  spatial_blocks = length(blocks$block_levels),
  replicates = n_replicates,
  permutations_per_condition = nperm,
  target_type1_rate = mean(replicates$null_target_p <= 0.05),
  target_type1_ci_lower = type1_interval[[1L]],
  target_type1_ci_upper = type1_interval[[2L]],
  null_familywise_rate_bh = mean(replicates$null_any_discovery),
  implanted_target_power_bh = mean(replicates$alternative_target_discovery),
  empirical_fdr = mean(replicates$alternative_fdp),
  median_null_target_p = stats::median(replicates$null_target_p),
  median_alternative_target_p = stats::median(replicates$alternative_target_p),
  wall_seconds = wall_seconds
)

write.table(replicates, file.path(out_dir, "replicates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(
  list(graph = graph, support = support, blocks = blocks,
       group = group, components = components),
  file.path(out_dir, "fixed_spatial_design.rds"), compress = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)

