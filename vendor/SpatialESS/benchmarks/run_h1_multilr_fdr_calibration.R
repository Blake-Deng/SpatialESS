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

env_double <- function(name, default, lower = -Inf, upper = Inf) {
  value <- suppressWarnings(as.numeric(
    Sys.getenv(name, unset = as.character(default))
  ))
  if (length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper) {
    stop(sprintf("%s must be one finite number in [%s, %s].",
                 name, lower, upper), call. = FALSE)
  }
  value
}

n_replicates <- env_integer("SPATIALESS_MULTILR_REPLICATES", 30L)
nperm <- env_integer("SPATIALESS_MULTILR_NPERM", 499L)
seed <- env_integer("SPATIALESS_MULTILR_SEED", 20260731L)
grid_side <- env_integer("SPATIALESS_MULTILR_GRID_SIDE", 20L)
block_side <- env_integer("SPATIALESS_MULTILR_BLOCK_SIDE", 5L)
zero_fraction <- env_double(
  "SPATIALESS_MULTILR_ZERO_FRACTION", 0, lower = 0, upper = 0.95
)
tie_step <- env_double(
  "SPATIALESS_MULTILR_TIE_STEP", 0, lower = 0, upper = Inf
)
lr_dependence <- env_double(
  "SPATIALESS_MULTILR_DEPENDENCE", 0, lower = 0, upper = 0.99
)
screening_margin <- env_double(
  "SPATIALESS_MULTILR_SCREENING_MARGIN", 0.01,
  lower = 0, upper = 0.949999
)
analytic_max_zero_fraction <- env_double(
  "SPATIALESS_MULTILR_MAX_ZERO_FRACTION", 0.10,
  lower = 0, upper = 0.95
)
lr_count <- 10L
signal_count <- 2L
out_dir <- Sys.getenv(
  "SPATIALESS_MULTILR_OUT",
  file.path(
    project_dir, "benchmarks", "results",
    sprintf("h1_multilr_conservative_rep%d_perm%d", n_replicates, nperm)
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
  coordinates, block_size = block_side, origin = c(0, 0)
)

simple_ligand <- paste0("L", seq_len(8L))
simple_receptor <- paste0("R", seq_len(8L))
all_genes <- c(
  simple_ligand, simple_receptor,
  "L9a", "L9b", "R9", "L10", "R10", "CO10"
)
lr <- data.frame(
  interaction_name = paste0("LR", seq_len(lr_count)),
  ligand = c(simple_ligand, "CX_L", "L10"),
  receptor = c(simple_receptor, "R9", "R10"),
  co_A_receptor = c(rep("", 9L), "CO10"),
  stringsAsFactors = FALSE
)
complex <- data.frame(
  subunit_1 = "L9a", subunit_2 = "L9b",
  row.names = "CX_L", stringsAsFactors = FALSE
)
cofactor <- data.frame(
  cofactor1 = "CO10", row.names = "CO10", stringsAsFactors = FALSE
)
components <- prepare_cellchat_lr_components(
  lr, all_genes, complex, cofactor
)
simple_lr_indices <- seq_len(8L)
control_lr_indices <- c(9L, 10L)
target_keys <- paste("A", "B", sep = "->")

simulate_expression <- function(replicate_seed, effect = 1) {
  set.seed(replicate_seed)
  cell_count <- nrow(coordinates)
  block_count <- length(blocks$block_levels)
  block_baseline <- matrix(
    stats::rgamma(block_count * length(all_genes), shape = 4, rate = 4),
    nrow = length(all_genes), ncol = block_count
  )
  expression <- matrix(0, nrow = length(all_genes), ncol = cell_count,
                       dimnames = list(all_genes, rownames(coordinates)))
  common_profile <- if (lr_dependence > 0) {
    stats::rgamma(cell_count, shape = 7, rate = 7)
  } else {
    NULL
  }
  for (gene in seq_along(all_genes)) {
    independent_profile <- stats::rgamma(cell_count, shape = 7, rate = 7)
    profile <- if (lr_dependence > 0) {
      (1 - lr_dependence) * independent_profile +
        lr_dependence * common_profile
    } else {
      independent_profile
    }
    expression[gene, ] <- block_baseline[gene, blocks$block_code] * profile
  }
  if (zero_fraction > 0) {
    expression[
      matrix(stats::runif(length(expression)) < zero_fraction,
             nrow = nrow(expression))
    ] <- 0
  }
  if (tie_step > 0) {
    expression <- round(expression / tie_step) * tie_step
  }
  if (effect != 1) {
    for (signal in seq_len(signal_count)) {
      ligand_index <- match(simple_ligand[[signal]], all_genes)
      receptor_index <- match(simple_receptor[[signal]], all_genes)
      expression[ligand_index, group == "A"] <-
        expression[ligand_index, group == "A"] * effect
      expression[receptor_index, group == "B"] <-
        expression[receptor_index, group == "B"] * effect
    }
  }
  Matrix(expression, sparse = TRUE)
}

extract_target <- function(expression, permutation_seed) {
  prepared <- prepare_sparse_trimean(
    expression, genes = components$genes, normalize = FALSE
  )
  hybrid <- experimental_hybrid_cellchat_group_support(
    prepared, group, support, components, blocks,
    nperm = nperm, seed = permutation_seed,
    alpha = 0.05, screening_margin = screening_margin,
    finite_correction = TRUE, tail = "greater_equal",
    min_group_size = 30L, min_shared_blocks = 5L,
    max_zero_fraction = analytic_max_zero_fraction
  )
  permutation <- permutation_cellchat_group_support(
    prepared, group, support, components, blocks,
    nperm = nperm, seed = permutation_seed,
    finite_correction = TRUE, tail = "greater_equal", verbose = FALSE
  )
  select_target <- function(x) {
    required <- c("sender_group_name", "receiver_group_name")
    if (!nrow(x) || !all(required %in% colnames(x))) {
      return(x[FALSE, , drop = FALSE])
    }
    x[
      x$sender_group_name == "A" & x$receiver_group_name == "B",
      , drop = FALSE
    ]
  }
  hybrid_target <- select_target(hybrid$group_pairs)
  permutation_target <- select_target(permutation$group_pairs)
  if (anyDuplicated(hybrid_target$lr_index) ||
      anyDuplicated(permutation_target$lr_index)) {
    stop("A->B LR records contain duplicate LR indices.", call. = FALSE)
  }
  lr_index <- seq_len(lr_count)
  hybrid_index <- match(lr_index, hybrid_target$lr_index)
  permutation_index <- match(lr_index, permutation_target$lr_index)
  hybrid_present <- !is.na(hybrid_index)
  permutation_present <- !is.na(permutation_index)

  hybrid_p <- rep(1, lr_count)
  analytic_p <- rep(1, lr_count)
  permutation_p <- rep(1, lr_count)
  inference_mode <- rep("inactive_zero_probability", lr_count)
  fallback_reason <- rep("inactive_zero_probability", lr_count)
  permutation_confirmed <- rep(FALSE, lr_count)

  hybrid_p[hybrid_present] <-
    hybrid_target$pvalue[hybrid_index[hybrid_present]]
  analytic_p[hybrid_present] <-
    hybrid_target$analytic_pvalue[hybrid_index[hybrid_present]]
  inference_mode[hybrid_present] <-
    hybrid_target$inference_mode[hybrid_index[hybrid_present]]
  fallback_reason[hybrid_present] <-
    hybrid_target$fallback_reason[hybrid_index[hybrid_present]]
  permutation_confirmed[hybrid_present] <-
    hybrid_target$permutation_confirmed[hybrid_index[hybrid_present]]
  permutation_p[permutation_present] <-
    permutation_target$pvalue[permutation_index[permutation_present]]

  screened <- inference_mode ==
    "experimental_hybrid_screened_nonsignificant"
  analytic_eligible <- inference_mode %in% c(
    "experimental_hybrid_screened_nonsignificant",
    "experimental_hybrid_permutation_candidate"
  )
  inactive <- inference_mode == "inactive_zero_probability"
  if (any(hybrid_p + 1e-15 < permutation_p)) {
    stop("Hybrid p-values must be coordinatewise no smaller than full permutation.")
  }
  list(
    analytic_p = analytic_p,
    permutation_p = permutation_p,
    hybrid_p = hybrid_p,
    analytic_eligible = analytic_eligible,
    fallback = !analytic_eligible & !inactive,
    inactive = inactive,
    screened = screened,
    permutation_confirmed = permutation_confirmed,
    inference_mode = inference_mode,
    fallback_reason = fallback_reason
  )
}

run_one <- function(replicate_index) {
  data_seed <- seed + replicate_index * 3L
  permutation_seed <- seed + n_replicates * 3L +
    (replicate_index - 1L) * nperm * 3L
  null <- extract_target(
    simulate_expression(data_seed, effect = 1), permutation_seed
  )
  alternative <- extract_target(
    simulate_expression(data_seed + 1L, effect = 1.8),
    permutation_seed + nperm
  )
  null_analytic_input <- null$analytic_p
  null_analytic_input[!null$analytic_eligible] <- 1
  null_analytic_q <- p.adjust(null_analytic_input, method = "BH")
  null_hybrid_q <- p.adjust(null$hybrid_p, method = "BH")
  null_permutation_q <- p.adjust(null$permutation_p, method = "BH")
  alternative_analytic_input <- alternative$analytic_p
  alternative_analytic_input[!alternative$analytic_eligible] <- 1
  alternative_analytic_q <- p.adjust(
    alternative_analytic_input, method = "BH"
  )
  alternative_hybrid_q <- p.adjust(alternative$hybrid_p, method = "BH")
  alternative_permutation_q <- p.adjust(
    alternative$permutation_p, method = "BH"
  )
  null_analytic_discovery <- sum(null_analytic_q <= 0.05)
  null_hybrid_discovery <- null_hybrid_q <= 0.05
  null_permutation_discovery <- null_permutation_q <= 0.05
  alternative_analytic_signal <- alternative_analytic_q <= 0.05
  alternative_hybrid_discovery <- alternative_hybrid_q <= 0.05
  alternative_permutation_discovery <- alternative_permutation_q <= 0.05
  alternative_analytic_signal <-
    alternative_analytic_signal[seq_len(signal_count)]
  alternative_hybrid_signal <-
    alternative_hybrid_discovery[seq_len(signal_count)]
  alternative_permutation_signal <-
    alternative_permutation_discovery[seq_len(signal_count)]
  alternative_hybrid_false <- sum(
    alternative_hybrid_discovery[-seq_len(signal_count)]
  )
  alternative_permutation_false <- sum(
    alternative_permutation_discovery[-seq_len(signal_count)]
  )
  data.frame(
    replicate = replicate_index,
    analytic_completion_rate = mean(
      null$analytic_eligible & alternative$analytic_eligible
    ),
    null_fallback_records = sum(null$fallback),
    null_inactive_records = sum(null$inactive),
    null_screened_records = sum(null$screened),
    null_permutation_records = sum(null$permutation_confirmed),
    null_analytic_fwer = null_analytic_discovery > 0,
    null_hybrid_fwer = any(null_hybrid_discovery),
    null_permutation_fwer = any(null_permutation_discovery),
    null_hybrid_discoveries = sum(null_hybrid_discovery),
    null_permutation_discoveries = sum(null_permutation_discovery),
    null_hybrid_threshold_agreement = mean(
      null_hybrid_discovery == null_permutation_discovery
    ),
    null_hybrid_exact_rejection_set = identical(
      null_hybrid_discovery, null_permutation_discovery
    ),
    null_all_discoveries_permutation_confirmed = all(
      null$permutation_confirmed[null_hybrid_discovery]
    ),
    null_simple_mean_abs_p_difference = mean(abs(
      null$analytic_p[null$analytic_eligible] -
        null$permutation_p[null$analytic_eligible]
    )),
    null_simple_p_spearman = suppressWarnings(stats::cor(
      null$analytic_p[null$analytic_eligible],
      null$permutation_p[null$analytic_eligible],
      method = "spearman"
    )),
    alternative_fallback_records = sum(alternative$fallback),
    alternative_inactive_records = sum(alternative$inactive),
    alternative_screened_records = sum(alternative$screened),
    alternative_permutation_records = sum(alternative$permutation_confirmed),
    alternative_analytic_signal_power = mean(alternative_analytic_signal),
    alternative_hybrid_signal_power = mean(alternative_hybrid_signal),
    alternative_permutation_signal_power = mean(alternative_permutation_signal),
    alternative_hybrid_fdp = alternative_hybrid_false /
      max(1L, sum(alternative_hybrid_discovery)),
    alternative_permutation_fdp = alternative_permutation_false /
      max(1L, sum(alternative_permutation_discovery)),
    alternative_hybrid_threshold_agreement = mean(
      alternative_hybrid_discovery == alternative_permutation_discovery
    ),
    alternative_hybrid_exact_rejection_set = identical(
      alternative_hybrid_discovery, alternative_permutation_discovery
    ),
    alternative_all_discoveries_permutation_confirmed = all(
      alternative$permutation_confirmed[alternative_hybrid_discovery]
    ),
    stringsAsFactors = FALSE
  )
}

message(sprintf(
  "Running %d multi-LR null/alternative replicates with %d permutations...",
  n_replicates, nperm
))
workflow_start <- proc.time()[["elapsed"]]
replicates <- lapply(seq_len(n_replicates), function(index) {
  if (index == 1L || index == n_replicates || index %% 5L == 0L) {
    message(sprintf("Multi-LR calibration replicate %d/%d", index, n_replicates))
  }
  run_one(index)
})
replicates <- do.call(rbind, replicates)
wall_seconds <- proc.time()[["elapsed"]] - workflow_start
summary <- data.frame(
  cells = nrow(coordinates), groups = nlevels(group),
  grid_side = grid_side, block_side = block_side,
  zero_fraction = zero_fraction, tie_step = tie_step,
  lr_dependence = lr_dependence,
  analytic_max_zero_fraction = analytic_max_zero_fraction,
  spatial_blocks = length(blocks$block_levels), lr_records = lr_count,
  simple_lr_records = length(simple_lr_indices),
  complex_or_cofactor_controls = length(control_lr_indices),
  implanted_signal_records = signal_count,
  replicates = n_replicates, permutations_per_condition = nperm,
  alpha = 0.05, screening_margin = screening_margin,
  analytic_completion_rate = mean(replicates$analytic_completion_rate),
  null_analytic_fwer = mean(replicates$null_analytic_fwer),
  null_hybrid_fwer = mean(replicates$null_hybrid_fwer),
  null_permutation_fwer = mean(replicates$null_permutation_fwer),
  null_hybrid_mean_discoveries = mean(replicates$null_hybrid_discoveries),
  null_mean_screened_records = mean(replicates$null_screened_records),
  null_mean_permutation_records = mean(replicates$null_permutation_records),
  null_hybrid_threshold_agreement = mean(
    replicates$null_hybrid_threshold_agreement
  ),
  null_exact_rejection_set_rate = mean(
    replicates$null_hybrid_exact_rejection_set
  ),
  null_all_discoveries_permutation_confirmed = all(
    replicates$null_all_discoveries_permutation_confirmed
  ),
  null_permutation_mean_discoveries = mean(
    replicates$null_permutation_discoveries
  ),
  null_simple_mean_abs_p_difference = mean(
    replicates$null_simple_mean_abs_p_difference
  ),
  null_simple_p_spearman = mean(replicates$null_simple_p_spearman),
  alternative_analytic_signal_power = mean(
    replicates$alternative_analytic_signal_power
  ),
  alternative_hybrid_signal_power = mean(
    replicates$alternative_hybrid_signal_power
  ),
  alternative_permutation_signal_power = mean(
    replicates$alternative_permutation_signal_power
  ),
  alternative_hybrid_mean_fdp = mean(replicates$alternative_hybrid_fdp),
  alternative_permutation_mean_fdp = mean(
    replicates$alternative_permutation_fdp
  ),
  alternative_hybrid_threshold_agreement = mean(
    replicates$alternative_hybrid_threshold_agreement
  ),
  alternative_mean_screened_records = mean(
    replicates$alternative_screened_records
  ),
  alternative_mean_permutation_records = mean(
    replicates$alternative_permutation_records
  ),
  alternative_exact_rejection_set_rate = mean(
    replicates$alternative_hybrid_exact_rejection_set
  ),
  alternative_all_discoveries_permutation_confirmed = all(
    replicates$alternative_all_discoveries_permutation_confirmed
  ),
  all_replicates_exact_rejection_sets = all(
    replicates$null_hybrid_exact_rejection_set &
      replicates$alternative_hybrid_exact_rejection_set
  ),
  wall_seconds = wall_seconds
)
write.table(replicates, file.path(out_dir, "replicates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(
  list(graph = graph, support = support, blocks = blocks,
       group = group, components = components, lr = lr),
  file.path(out_dir, "fixed_spatial_design.rds"), compress = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(summary)
