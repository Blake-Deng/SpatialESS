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

n_replicates <- env_integer("SPATIALESS_ANALYTIC_REPLICATES", 200L)
nperm <- env_integer("SPATIALESS_ANALYTIC_NPERM", 999L)
seed <- env_integer("SPATIALESS_ANALYTIC_SEED", 20260731L)
grid_side <- env_integer("SPATIALESS_ANALYTIC_GRID_SIDE", 20L)
if (grid_side %% 2L != 0L || grid_side < 10L) {
  stop("SPATIALESS_ANALYTIC_GRID_SIDE must be even and at least 10.",
       call. = FALSE)
}
out_dir <- Sys.getenv(
  "SPATIALESS_ANALYTIC_OUT",
  file.path(
    project_dir, "benchmarks", "results",
    sprintf("h1_analytic_simple_rep%d_perm%d", n_replicates, nperm)
  )
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coordinates <- as.matrix(expand.grid(
  x = seq_len(grid_side) - 1L,
  y = seq_len(grid_side) - 1L
))
rownames(coordinates) <- paste0("cell", seq_len(nrow(coordinates)))
group <- factor(
  ifelse((coordinates[, "x"] + coordinates[, "y"]) %% 2L == 0L, "A", "B"),
  levels = c("A", "B")
)
graph <- build_radius_graph_csr(
  coordinates, radius = sqrt(2) + 1e-8,
  weight = "binary", store_distance = FALSE
)
support <- build_group_support_csr(graph, group)
blocks <- build_spatial_blocks(
  coordinates, block_size = grid_side + 1,
  origin = c(0, 0)
)
complex <- data.frame(subunit_1 = character(), row.names = character())
cofactor <- data.frame(cofactor1 = character(), row.names = character())
lr <- data.frame(interaction_name = "L_R", ligand = "L", receptor = "R")
components <- prepare_cellchat_lr_components(
  lr, c("L", "R"), complex, cofactor
)

simulate_expression <- function(replicate_seed, effect = 1) {
  set.seed(replicate_seed)
  ligand <- stats::rgamma(nrow(coordinates), shape = 5, rate = 2)
  receptor <- stats::rgamma(nrow(coordinates), shape = 6, rate = 2)
  ligand[group == "A"] <- ligand[group == "A"] * effect
  receptor[group == "B"] <- receptor[group == "B"] * effect
  Matrix(
    rbind(L = ligand, R = receptor), sparse = TRUE,
    dimnames = list(c("L", "R"), rownames(coordinates))
  )
}

extract_target <- function(expression, permutation_seed) {
  prepared <- prepare_sparse_trimean(
    expression, genes = components$genes, normalize = FALSE
  )
  analytic <- experimental_analytic_cellchat_group_support(
    prepared, group, support, components, blocks
  )
  permutation <- permutation_cellchat_group_support(
    prepared, group, support, components, blocks,
    nperm = nperm, seed = permutation_seed,
    finite_correction = TRUE, tail = "greater_equal",
    verbose = FALSE
  )
  analytic_target <- subset(
    analytic$group_pairs,
    sender_group_name == "A" & receiver_group_name == "B"
  )
  permutation_target <- subset(
    permutation$group_pairs,
    sender_group_name == "A" & receiver_group_name == "B"
  )
  if (nrow(analytic_target) != 1L || nrow(permutation_target) != 1L) {
    stop("A->B target record was not emitted exactly once.", call. = FALSE)
  }
  data.frame(
    analytic_p = analytic_target$pvalue,
    permutation_p = permutation_target$pvalue,
    inference_mode = analytic_target$inference_mode,
    fallback_reason = analytic_target$fallback_reason %||% NA_character_,
    analytic_z = analytic_target$analytic_z
  )
}

`%||%` <- function(x, y) if (length(x) && !is.na(x)) x else y

run_one <- function(replicate_index) {
  data_seed <- seed + replicate_index * 2L
  permutation_seed <- seed + n_replicates * 2L +
    (replicate_index - 1L) * nperm * 2L
  if (as.double(permutation_seed) + 2 * nperm > .Machine$integer.max) {
    stop("Seed sequence exceeds the R integer range.", call. = FALSE)
  }
  null <- extract_target(
    simulate_expression(data_seed, effect = 1), permutation_seed
  )
  alternative <- extract_target(
    simulate_expression(data_seed + 1L, effect = 1.35),
    permutation_seed + nperm
  )
  data.frame(
    replicate = replicate_index,
    null_analytic_p = null$analytic_p,
    null_permutation_p = null$permutation_p,
    null_analytic_z = null$analytic_z,
    null_mode = null$inference_mode,
    null_fallback = null$fallback_reason,
    alternative_analytic_p = alternative$analytic_p,
    alternative_permutation_p = alternative$permutation_p,
    alternative_analytic_z = alternative$analytic_z,
    alternative_mode = alternative$inference_mode,
    alternative_fallback = alternative$fallback_reason,
    stringsAsFactors = FALSE
  )
}

message(sprintf(
  "Running %d analytic/permutation null and alternative replicates with %d permutations...",
  n_replicates, nperm
))
workflow_start <- proc.time()[["elapsed"]]
replicates <- lapply(seq_len(n_replicates), function(index) {
  if (index == 1L || index == n_replicates || index %% 10L == 0L) {
    message(sprintf("Analytic calibration replicate %d/%d", index, n_replicates))
  }
  run_one(index)
})
replicates <- do.call(rbind, replicates)
wall_seconds <- proc.time()[["elapsed"]] - workflow_start

alpha <- 0.05
null_analytic_reject <- replicates$null_analytic_p <= alpha
null_permutation_reject <- replicates$null_permutation_p <= alpha
alternative_analytic_reject <- replicates$alternative_analytic_p <= alpha
alternative_permutation_reject <- replicates$alternative_permutation_p <= alpha
analytic_null_interval <- stats::binom.test(
  sum(null_analytic_reject), n_replicates
)$conf.int
permutation_null_interval <- stats::binom.test(
  sum(null_permutation_reject), n_replicates
)$conf.int

summary <- data.frame(
  cells = nrow(coordinates),
  groups = nlevels(group),
  graph_edges = length(graph$neighbors),
  spatial_blocks = length(blocks$block_levels),
  replicates = n_replicates,
  permutations_per_condition = nperm,
  analytic_completion_rate = mean(
    replicates$null_mode == "experimental_analytic" &
      replicates$alternative_mode == "experimental_analytic"
  ),
  analytic_type1_rate = mean(null_analytic_reject),
  analytic_type1_ci_lower = analytic_null_interval[[1L]],
  analytic_type1_ci_upper = analytic_null_interval[[2L]],
  permutation_type1_rate = mean(null_permutation_reject),
  permutation_type1_ci_lower = permutation_null_interval[[1L]],
  permutation_type1_ci_upper = permutation_null_interval[[2L]],
  null_threshold_agreement = mean(
    null_analytic_reject == null_permutation_reject
  ),
  null_mean_abs_p_difference = mean(abs(
    replicates$null_analytic_p - replicates$null_permutation_p
  )),
  null_max_abs_p_difference = max(abs(
    replicates$null_analytic_p - replicates$null_permutation_p
  )),
  analytic_power = mean(alternative_analytic_reject),
  permutation_power = mean(alternative_permutation_reject),
  alternative_threshold_agreement = mean(
    alternative_analytic_reject == alternative_permutation_reject
  ),
  alternative_mean_abs_p_difference = mean(abs(
    replicates$alternative_analytic_p - replicates$alternative_permutation_p
  )),
  alternative_max_abs_p_difference = max(abs(
    replicates$alternative_analytic_p - replicates$alternative_permutation_p
  )),
  null_p_spearman = suppressWarnings(stats::cor(
    replicates$null_analytic_p, replicates$null_permutation_p,
    method = "spearman"
  )),
  wall_seconds = wall_seconds
)

write.table(replicates, file.path(out_dir, "replicates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(
  list(graph = graph, support = support, blocks = blocks,
       group = group, components = components),
  file.path(out_dir, "fixed_design.rds"), compress = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
