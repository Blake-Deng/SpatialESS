.spatialess_weighted_quantile <- function(values, weights, probabilities) {
  keep <- is.finite(weights) & weights > 0
  values <- values[keep]
  weights <- weights[keep]
  if (!length(values) || any(!is.finite(values)) || any(values < 0)) {
    return(rep(NA_real_, length(probabilities)))
  }
  if (max(weights) - min(weights) <= .Machine$double.eps * max(weights)) {
    return(as.numeric(stats::quantile(
      values, probabilities, names = FALSE, type = 7
    )))
  }
  order_index <- order(values, method = "radix")
  values <- values[order_index]
  weights <- weights[order_index]
  duplicate_value <- duplicated(values) | duplicated(values, fromLast = TRUE)
  if (any(duplicate_value)) {
    unique_values <- unique(values)
    weights <- as.numeric(rowsum(weights, match(values, unique_values), reorder = FALSE))
    values <- unique_values
  }
  weights <- weights / sum(weights)
  if (length(values) == 1L) return(rep(values, length(probabilities)))
  left_mass <- c(0, utils::head(cumsum(weights), -1L))
  plotting_position <- left_mass / pmax(.Machine$double.eps, 1 - weights)
  plotting_position[[1L]] <- 0
  plotting_position[[length(plotting_position)]] <- 1
  as.numeric(stats::approx(
    x = plotting_position, y = values, xout = probabilities,
    method = "linear", ties = "ordered", rule = 2
  )$y)
}

.spatialess_stratified_design <- function(block_code, group_code,
                                          block_count, group_count) {
  block_sizes <- tabulate(block_code, nbins = block_count)
  encoded <- block_code + (group_code - 1L) * block_count
  block_group_counts <- matrix(
    tabulate(encoded, nbins = block_count * group_count),
    nrow = block_count, ncol = group_count
  )
  list(
    block_sizes = block_sizes,
    block_group_counts = block_group_counts,
    group_sizes = colSums(block_group_counts)
  )
}

.spatialess_block_variance <- function(influence, group_index, design,
                                       block_members) {
  group_size <- design$group_sizes[[group_index]]
  total <- 0
  for (block in seq_along(block_members)) {
    members <- block_members[[block]]
    block_size <- design$block_sizes[[block]]
    selected <- design$block_group_counts[block, group_index]
    if (block_size <= 1L || selected <= 0L || selected >= block_size) next
    block_variance <- stats::var(influence[members])
    if (is.finite(block_variance)) {
      total <- total + selected * (1 - selected / block_size) * block_variance
    }
  }
  total / group_size^2
}

.spatialess_block_cross_covariance <- function(
    ligand_influence, receptor_influence, sender, receiver,
    design, block_members) {
  sender_size <- design$group_sizes[[sender]]
  receiver_size <- design$group_sizes[[receiver]]
  total <- 0
  for (block in seq_along(block_members)) {
    members <- block_members[[block]]
    block_size <- design$block_sizes[[block]]
    sender_count <- design$block_group_counts[block, sender]
    receiver_count <- design$block_group_counts[block, receiver]
    if (block_size <= 1L || sender_count <= 0L || receiver_count <= 0L) next
    block_covariance <- stats::cov(
      ligand_influence[members], receptor_influence[members]
    )
    if (is.finite(block_covariance)) {
      total <- total - sender_count * receiver_count / block_size *
        block_covariance
    }
  }
  total / (sender_size * receiver_size)
}

.spatialess_block_residual_correlation <- function(
    ligand_influence, receptor_influence, block_members) {
  ligand_residual <- ligand_influence
  receptor_residual <- receptor_influence
  for (members in block_members) {
    ligand_residual[members] <- ligand_residual[members] -
      mean(ligand_residual[members])
    receptor_residual[members] <- receptor_residual[members] -
      mean(receptor_residual[members])
  }
  suppressWarnings(stats::cor(ligand_residual, receptor_residual))
}

#' Experimental block-stratified analytic significance for simple LR records
#'
#' This Stage 1B prototype linearizes CellChat's joint Q1/median/Q3 triMean
#' under the same fixed-capacity spatial-block randomization used by
#' [permutation_cellchat_group_support()]. Null centers are group-specific
#' weighted empirical quantiles. Variances and ligand-receptor covariance are
#' sums of finite-population block contributions, with cell profiles never
#' exchanged across sample, compartment or spatial block.
#'
#' The method remains experimental. Unsupported records receive an explicit
#' permutation fallback reason rather than an approximate p-value.
#'
#' @param prepared A result from [prepare_sparse_trimean()].
#' @param group Fixed target-group labels.
#' @param group_support Spatially supported group pairs.
#' @param components Resolved CellChat LR components.
#' @param blocks Spatial exchangeability blocks.
#' @param Kh,n Positive Hill-function parameters used for the observed score.
#' @param min_group_size Minimum cells in sender and receiver groups.
#' @param min_shared_blocks Minimum blocks containing both groups.
#' @param density_window Probability half-width for quantile-density estimates.
#' @param max_tie_fraction Maximum weighted mass at a triMean quartile.
#' @param max_zero_fraction Maximum zero fraction allowed in either expression
#'   profile before requiring permutation fallback.
#' @param max_abs_influence_correlation Maximum block-residual ligand-receptor
#'   influence correlation.
#' @param min_log_variance Minimum delta-method log-product variance.
#' @param max_influence_entries Maximum cached cell-by-gene-by-group entries.
#' @export
experimental_stratified_analytic_cellchat_group_support <- function(
    prepared, group, group_support, components, blocks,
    Kh = 0.5, n = 1, min_group_size = 30L, min_shared_blocks = 5L,
    density_window = 0.05, max_tie_fraction = 0.10,
    max_zero_fraction = 0.10, max_abs_influence_correlation = 0.80,
    min_log_variance = 1e-12, max_influence_entries = 5e7) {
  if (!inherits(prepared, "SpatialESSPreparedTriMean") ||
      !inherits(group_support, "SpatialESSGroupSupport") ||
      !inherits(components, "SpatialESSLRComponents") ||
      !inherits(blocks, "SpatialESSBlocks")) {
    stop("Invalid prepared, group-support, LR-component or block object.",
         call. = FALSE)
  }
  if (!identical(prepared$cell_id, blocks$cell_id)) {
    stop("Prepared expression and spatial blocks have different cells/order.",
         call. = FALSE)
  }
  if (any(!components$genes %in% prepared$genes)) {
    stop("Prepared triMean kernel is missing referenced LR genes.",
         call. = FALSE)
  }
  if (length(Kh) != 1L || !is.finite(Kh) || Kh <= 0 ||
      length(n) != 1L || !is.finite(n) || n <= 0) {
    stop("Kh and n must be positive finite scalars.", call. = FALSE)
  }
  min_group_size <- as.integer(min_group_size)
  min_shared_blocks <- as.integer(min_shared_blocks)
  if (length(min_group_size) != 1L || is.na(min_group_size) ||
      min_group_size < 2L || length(min_shared_blocks) != 1L ||
      is.na(min_shared_blocks) || min_shared_blocks < 1L) {
    stop("Group and shared-block minima must be positive valid integers.",
         call. = FALSE)
  }
  validate_probability <- function(value, name, lower_open = FALSE) {
    lower_valid <- if (lower_open) value > 0 else value >= 0
    if (length(value) != 1L || !is.finite(value) || !lower_valid || value >= 1) {
      stop(sprintf("%s is outside its valid probability range.", name),
           call. = FALSE)
    }
    as.double(value)
  }
  density_window <- validate_probability(
    density_window, "density_window", lower_open = TRUE
  )
  if (density_window >= 0.25) {
    stop("density_window must be smaller than 0.25.", call. = FALSE)
  }
  max_tie_fraction <- validate_probability(
    max_tie_fraction, "max_tie_fraction"
  )
  max_zero_fraction <- validate_probability(
    max_zero_fraction, "max_zero_fraction"
  )
  max_abs_influence_correlation <- validate_probability(
    max_abs_influence_correlation, "max_abs_influence_correlation",
    lower_open = TRUE
  )
  if (length(min_log_variance) != 1L || !is.finite(min_log_variance) ||
      min_log_variance <= 0 || length(max_influence_entries) != 1L ||
      !is.finite(max_influence_entries) || max_influence_entries < 1) {
    stop("Variance and workspace limits must be positive finite scalars.",
         call. = FALSE)
  }

  info <- .spatialess_group_info(group)
  cell_count <- length(prepared$cell_id)
  group_levels <- group_support$group_levels
  if (length(info$code) != cell_count || info$count != length(group_levels) ||
      !identical(as.character(info$levels[seq_len(info$count)]),
                 as.character(group_levels))) {
    stop("Target groups do not align with prepared cells/group support.",
         call. = FALSE)
  }
  block_count <- length(blocks$block_levels)
  design <- .spatialess_stratified_design(
    blocks$block_code, info$code, block_count, info$count
  )
  block_members <- split(seq_len(cell_count), blocks$block_code)

  observed_expression <- summarize_prepared_trimean(
    prepared, info$code, group_levels = group_levels
  )
  observed <- .spatialess_group_score_matrix(
    observed_expression$average, group_support, components, Kh, n
  )
  active <- which(observed > 0, arr.ind = TRUE)
  support <- group_support$group_pairs
  lr <- components$lr
  interaction <- if ("interaction_name" %in% colnames(lr)) {
    as.character(lr$interaction_name)
  } else {
    paste(lr$ligand, lr$receptor, sep = "_")
  }
  if (!nrow(active)) {
    return(structure(
      list(
        group_pairs = data.frame(),
        parameters = list(calibration_status = "experimental_unvalidated"),
        diagnostics = list(active_records = 0L, analytic_records = 0L,
                           fallback_records = 0L)
      ),
      class = "SpatialESSExperimentalStratifiedAnalytic"
    ))
  }

  simple_lr <- lengths(components$ligand) == 1L &
    lengths(components$receptor) == 1L &
    lengths(components$co_a) == 0L & lengths(components$co_i) == 0L &
    lengths(components$agonist) == 0L &
    lengths(components$antagonist) == 0L &
    !components$has_agonist & !components$has_antagonist
  output <- data.frame(
    lr_index = active[, "col"],
    interaction_name = interaction[active[, "col"]],
    ligand = as.character(lr$ligand)[active[, "col"]],
    receptor = as.character(lr$receptor)[active[, "col"]],
    support_pair_index = active[, "row"],
    sender_group = support$sender_group[active[, "row"]],
    receiver_group = support$receiver_group[active[, "row"]],
    sender_group_name = group_levels[support$sender_group[active[, "row"]]],
    receiver_group_name = group_levels[support$receiver_group[active[, "row"]]],
    supported_edges = support$supported_edges[active[, "row"]],
    probability = observed[active],
    pvalue = NA_real_, analytic_z = NA_real_, analytic_log_sd = NA_real_,
    influence_correlation = NA_real_, shared_blocks = 0L,
    inference_mode = "permutation_fallback_required",
    fallback_reason = NA_character_, stringsAsFactors = FALSE
  )

  candidate_keys <- character()
  for (row in seq_len(nrow(output))) {
    lr_index <- output$lr_index[[row]]
    if (!simple_lr[[lr_index]]) next
    candidate_keys <- c(
      candidate_keys,
      paste(components$ligand[[lr_index]][[1L]],
            output$sender_group[[row]], sep = "\r"),
      paste(components$receptor[[lr_index]][[1L]],
            output$receiver_group[[row]], sep = "\r")
    )
  }
  candidate_keys <- unique(candidate_keys)
  workspace_entries <- as.double(cell_count) * length(candidate_keys)
  workspace_available <- workspace_entries <= max_influence_entries
  influence_cache <- new.env(parent = emptyenv())

  build_influence <- function(gene, group_index) {
    key <- paste(gene, group_index, sep = "\r")
    if (exists(key, envir = influence_cache, inherits = FALSE)) {
      return(get(key, envir = influence_cache, inherits = FALSE))
    }
    values <- as.numeric(
      prepared$cell_by_gene[, match(gene, prepared$genes)]
    )
    group_size <- design$group_sizes[[group_index]]
    cell_weights <- design$block_group_counts[blocks$block_code, group_index] /
      design$block_sizes[blocks$block_code] / group_size
    probabilities <- c(0.25, 0.50, 0.75)
    trimean_weights <- c(0.25, 0.50, 0.25)
    quartiles <- .spatialess_weighted_quantile(
      values, cell_weights, probabilities
    )
    zero_fraction <- sum(cell_weights[values == 0]) / sum(cell_weights)
    reason <- if (zero_fraction > max_zero_fraction) {
      "excessive_zero_fraction"
    } else {
      NA_character_
    }
    if (is.na(reason) &&
        (any(!is.finite(quartiles)) || any(quartiles <= 0))) {
      reason <- "zero_or_invalid_quartile"
    }
    tie_fraction <- vapply(quartiles, function(value) {
      sum(cell_weights[values == value]) / sum(cell_weights)
    }, numeric(1))
    if (is.na(reason) && any(tie_fraction > max_tie_fraction)) {
      reason <- "excessive_quartile_ties"
    }
    lower <- .spatialess_weighted_quantile(
      values, cell_weights, probabilities - density_window
    )
    upper <- .spatialess_weighted_quantile(
      values, cell_weights, probabilities + density_window
    )
    density <- 2 * density_window / (upper - lower)
    if (is.na(reason) &&
        any(!is.finite(density) | density <= sqrt(.Machine$double.eps))) {
      reason <- "unstable_quantile_density"
    }
    if (is.na(reason)) {
      influence <- numeric(cell_count)
      for (j in seq_along(probabilities)) {
        influence <- influence + trimean_weights[[j]] *
          (probabilities[[j]] - as.double(values <= quartiles[[j]])) /
          density[[j]]
      }
    } else {
      influence <- numeric()
    }
    result <- list(
      center = sum(trimean_weights * quartiles), influence = influence,
      quartiles = quartiles, density = density,
      tie_fraction = tie_fraction, zero_fraction = zero_fraction,
      reason = reason
    )
    assign(key, result, envir = influence_cache)
    result
  }

  for (row in seq_len(nrow(output))) {
    lr_index <- output$lr_index[[row]]
    sender <- output$sender_group[[row]]
    receiver <- output$receiver_group[[row]]
    reason <- NA_character_
    if (!workspace_available) reason <- "influence_workspace_limit"
    if (is.na(reason) && !simple_lr[[lr_index]]) reason <- "non_simple_lr"
    if (is.na(reason) && sender == receiver) reason <- "autocrine_pair"
    if (is.na(reason) &&
        (design$group_sizes[[sender]] < min_group_size ||
         design$group_sizes[[receiver]] < min_group_size)) {
      reason <- "small_group"
    }
    shared <- sum(
      design$block_group_counts[, sender] > 0L &
        design$block_group_counts[, receiver] > 0L &
        design$block_sizes > 1L
    )
    output$shared_blocks[[row]] <- shared
    if (is.na(reason) && shared < min_shared_blocks) {
      reason <- "insufficient_shared_blocks"
    }
    if (!is.na(reason)) {
      output$fallback_reason[[row]] <- reason
      next
    }

    ligand_gene <- components$ligand[[lr_index]][[1L]]
    receptor_gene <- components$receptor[[lr_index]][[1L]]
    ligand <- build_influence(ligand_gene, sender)
    receptor <- build_influence(receptor_gene, receiver)
    if (!is.na(ligand$reason)) {
      output$fallback_reason[[row]] <- paste0("ligand_", ligand$reason)
      next
    }
    if (!is.na(receptor$reason)) {
      output$fallback_reason[[row]] <- paste0("receptor_", receptor$reason)
      next
    }

    influence_correlation <- .spatialess_block_residual_correlation(
      ligand$influence, receptor$influence, block_members
    )
    output$influence_correlation[[row]] <- influence_correlation
    if (!is.finite(influence_correlation) ||
        abs(influence_correlation) > max_abs_influence_correlation) {
      output$fallback_reason[[row]] <-
        "strong_block_residual_ligand_receptor_dependence"
      next
    }
    ligand_variance <- .spatialess_block_variance(
      ligand$influence, sender, design, block_members
    )
    receptor_variance <- .spatialess_block_variance(
      receptor$influence, receiver, design, block_members
    )
    cross_covariance <- .spatialess_block_cross_covariance(
      ligand$influence, receptor$influence, sender, receiver,
      design, block_members
    )
    log_variance <- ligand_variance / ligand$center^2 +
      receptor_variance / receptor$center^2 +
      2 * cross_covariance / (ligand$center * receptor$center)
    if (!is.finite(log_variance) || log_variance <= min_log_variance) {
      output$fallback_reason[[row]] <- "degenerate_delta_variance"
      next
    }
    observed_ligand <- observed_expression$average[ligand_gene, sender]
    observed_receptor <- observed_expression$average[receptor_gene, receiver]
    if (!(observed_ligand > 0 && observed_receptor > 0)) {
      output$fallback_reason[[row]] <- "zero_observed_component"
      next
    }
    z <- (log(observed_ligand) + log(observed_receptor) -
          log(ligand$center) - log(receptor$center)) / sqrt(log_variance)
    output$pvalue[[row]] <- stats::pnorm(z, lower.tail = FALSE)
    output$analytic_z[[row]] <- z
    output$analytic_log_sd[[row]] <- sqrt(log_variance)
    output$inference_mode[[row]] <- "experimental_stratified_analytic"
    output$fallback_reason[[row]] <- NA_character_
  }

  structure(
    list(
      group_pairs = output,
      parameters = list(
        calibration_status = "experimental_unvalidated",
        null = paste(
          "fixed-capacity profile permutation within sample x compartment x",
          "spatial block"
        ),
        method = paste(
          "group-weighted joint triMean influence functions with blockwise",
          "finite-population covariance and log-product delta method"
        ),
        Kh = Kh, n = n, min_group_size = min_group_size,
        min_shared_blocks = min_shared_blocks,
        density_window = density_window,
        max_tie_fraction = max_tie_fraction,
        max_zero_fraction = max_zero_fraction,
        max_abs_influence_correlation = max_abs_influence_correlation,
        min_log_variance = min_log_variance,
        max_influence_entries = max_influence_entries
      ),
      diagnostics = list(
        cells = cell_count, groups = info$count, blocks = block_count,
        active_records = nrow(output),
        analytic_records = sum(
          output$inference_mode == "experimental_stratified_analytic"
        ),
        fallback_records = sum(
          output$inference_mode == "permutation_fallback_required"
        ),
        fallback_reasons = sort(table(output$fallback_reason), decreasing = TRUE),
        candidate_influence_keys = length(candidate_keys),
        influence_workspace_entries = workspace_entries
      )
    ),
    class = "SpatialESSExperimentalStratifiedAnalytic"
  )
}
