#' Experimental analytic significance for simple CellChat LR records
#'
#' This prototype uses a joint influence-function approximation for the
#' ligand and receptor triMeans under one homogeneous finite-population label
#' permutation null. The covariance calculation retains both within-profile
#' ligand-receptor dependence and the negative dependence induced when two
#' different groups partition the same cells. Records outside the validated
#' simple-LR domain are returned with an explicit fallback reason and no
#' analytic p-value.
#'
#' This function is intentionally conservative and is not yet the default
#' SpatialESS inference mode. In particular, spatially stratified blocks,
#' autocrine records, complexes, cofactors, strong dependence, small groups,
#' zero quartiles and heavily tied quartiles require permutation fallback.
#'
#' @param prepared A result from [prepare_sparse_trimean()].
#' @param group Fixed target-group labels.
#' @param group_support Spatially supported group pairs.
#' @param components Resolved CellChat LR components.
#' @param blocks Spatial exchangeability blocks. Exactly one block is required
#'   by this first homogeneous-null prototype.
#' @param Kh,n Positive Hill-function parameters used for the observed score.
#' @param min_group_size Minimum cells in both sender and receiver groups.
#' @param density_window Probability half-width used for quantile-density
#'   estimation.
#' @param max_tie_fraction Maximum mass allowed at any triMean quartile.
#' @param max_zero_fraction Maximum zero fraction allowed in either expression
#'   profile before requiring permutation fallback.
#' @param max_abs_influence_correlation Dependence guard for ligand and
#'   receptor influence functions.
#' @param min_log_variance Minimum delta-method variance.
#' @param max_influence_entries Maximum cells times candidate genes retained in
#'   the influence-function cache.
#' @export
experimental_analytic_cellchat_group_support <- function(
    prepared, group, group_support, components, blocks,
    Kh = 0.5, n = 1, min_group_size = 30L,
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
  if (length(min_group_size) != 1L || is.na(min_group_size) ||
      min_group_size < 2L) {
    stop("min_group_size must be one integer of at least two.", call. = FALSE)
  }
  scalar_probability <- function(value, name, lower_open = FALSE) {
    valid_lower <- if (lower_open) value > 0 else value >= 0
    if (length(value) != 1L || !is.finite(value) || !valid_lower || value >= 1) {
      stop(sprintf("%s must be a finite scalar in %s.", name,
                   if (lower_open) "(0, 1)" else "[0, 1)"), call. = FALSE)
    }
    as.double(value)
  }
  density_window <- scalar_probability(
    density_window, "density_window", lower_open = TRUE
  )
  if (density_window >= 0.25) {
    stop("density_window must be smaller than 0.25.", call. = FALSE)
  }
  max_tie_fraction <- scalar_probability(max_tie_fraction, "max_tie_fraction")
  max_zero_fraction <- scalar_probability(
    max_zero_fraction, "max_zero_fraction"
  )
  max_abs_influence_correlation <- scalar_probability(
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
  if (length(info$code) != cell_count ||
      info$count != length(group_support$group_levels)) {
    stop("Target groups do not align with prepared cells/group support.",
         call. = FALSE)
  }
  group_levels <- group_support$group_levels
  group_sizes <- tabulate(info$code, nbins = info$count)
  if (!identical(as.character(info$levels[seq_len(info$count)]),
                 as.character(group_levels))) {
    stop("Target group levels differ from spatial group support.", call. = FALSE)
  }

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
        parameters = list(
          calibration_status = "experimental_unvalidated",
          null = "one homogeneous finite-population label permutation"
        ),
        diagnostics = list(active_records = 0L, analytic_records = 0L,
                           fallback_records = 0L)
      ),
      class = "SpatialESSExperimentalAnalytic"
    ))
  }

  simple_lr <- lengths(components$ligand) == 1L &
    lengths(components$receptor) == 1L &
    lengths(components$co_a) == 0L & lengths(components$co_i) == 0L &
    lengths(components$agonist) == 0L &
    lengths(components$antagonist) == 0L &
    !components$has_agonist & !components$has_antagonist
  candidate_genes <- unique(unlist(
    c(components$ligand[simple_lr], components$receptor[simple_lr]),
    use.names = FALSE
  ))
  workspace_entries <- as.double(cell_count) * length(candidate_genes)
  workspace_available <- workspace_entries <= max_influence_entries

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
    pvalue = NA_real_,
    analytic_z = NA_real_,
    analytic_log_sd = NA_real_,
    influence_correlation = NA_real_,
    inference_mode = "permutation_fallback_required",
    fallback_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  global_reason <- if (length(blocks$block_levels) != 1L) {
    "spatially_stratified_null"
  } else if (!workspace_available) {
    "influence_workspace_limit"
  } else {
    NA_character_
  }
  influence_cache <- new.env(parent = emptyenv())

  build_influence <- function(gene) {
    if (exists(gene, envir = influence_cache, inherits = FALSE)) {
      return(get(gene, envir = influence_cache, inherits = FALSE))
    }
    gene_index <- match(gene, prepared$genes)
    values <- as.numeric(prepared$cell_by_gene[, gene_index])
    probabilities <- c(0.25, 0.50, 0.75)
    weights <- c(0.25, 0.50, 0.25)
    quartiles <- as.numeric(stats::quantile(
      values, probabilities, names = FALSE, type = 7
    ))
    zero_fraction <- mean(values == 0)
    reason <- if (zero_fraction > max_zero_fraction) {
      "excessive_zero_fraction"
    } else {
      NA_character_
    }
    if (is.na(reason) &&
        (any(!is.finite(quartiles)) || any(quartiles <= 0))) {
      reason <- "zero_or_invalid_quartile"
    }
    tie_fraction <- vapply(
      quartiles, function(value) mean(values == value), numeric(1)
    )
    if (is.na(reason) && any(tie_fraction > max_tie_fraction)) {
      reason <- "excessive_quartile_ties"
    }
    lower_probability <- probabilities - density_window
    upper_probability <- probabilities + density_window
    lower_quantile <- as.numeric(stats::quantile(
      values, lower_probability, names = FALSE, type = 7
    ))
    upper_quantile <- as.numeric(stats::quantile(
      values, upper_probability, names = FALSE, type = 7
    ))
    quantile_span <- upper_quantile - lower_quantile
    density <- 2 * density_window / quantile_span
    if (is.na(reason) &&
        any(!is.finite(density) | density <= sqrt(.Machine$double.eps))) {
      reason <- "unstable_quantile_density"
    }
    if (is.na(reason)) {
      influence <- numeric(cell_count)
      for (j in seq_along(probabilities)) {
        influence <- influence + weights[j] *
          (probabilities[j] - as.double(values <= quartiles[j])) / density[j]
      }
      influence <- influence - mean(influence)
    } else {
      influence <- numeric()
    }
    result <- list(
      center = sum(weights * quartiles), influence = influence,
      density = density, quartiles = quartiles,
      tie_fraction = tie_fraction, zero_fraction = zero_fraction,
      reason = reason
    )
    assign(gene, result, envir = influence_cache)
    result
  }

  for (row in seq_len(nrow(output))) {
    lr_index <- output$lr_index[row]
    sender <- output$sender_group[row]
    receiver <- output$receiver_group[row]
    reason <- global_reason
    if (is.na(reason) && !simple_lr[lr_index]) reason <- "non_simple_lr"
    if (is.na(reason) && sender == receiver) reason <- "autocrine_pair"
    if (is.na(reason) &&
        (group_sizes[sender] < min_group_size ||
         group_sizes[receiver] < min_group_size)) {
      reason <- "small_group"
    }
    if (!is.na(reason)) {
      output$fallback_reason[row] <- reason
      next
    }

    ligand_gene <- components$ligand[[lr_index]][[1L]]
    receptor_gene <- components$receptor[[lr_index]][[1L]]
    ligand <- build_influence(ligand_gene)
    receptor <- build_influence(receptor_gene)
    if (!is.na(ligand$reason)) {
      output$fallback_reason[row] <- paste0("ligand_", ligand$reason)
      next
    }
    if (!is.na(receptor$reason)) {
      output$fallback_reason[row] <- paste0("receptor_", receptor$reason)
      next
    }

    influence_correlation <- suppressWarnings(stats::cor(
      ligand$influence, receptor$influence
    ))
    output$influence_correlation[row] <- influence_correlation
    if (!is.finite(influence_correlation) ||
        abs(influence_correlation) > max_abs_influence_correlation) {
      output$fallback_reason[row] <- "strong_ligand_receptor_dependence"
      next
    }

    ligand_variance_population <- stats::var(ligand$influence)
    receptor_variance_population <- stats::var(receptor$influence)
    cross_covariance_population <- stats::cov(
      ligand$influence, receptor$influence
    )
    ligand_variance <-
      (1 / group_sizes[sender] - 1 / cell_count) *
      ligand_variance_population
    receptor_variance <-
      (1 / group_sizes[receiver] - 1 / cell_count) *
      receptor_variance_population
    cross_covariance <- -cross_covariance_population / cell_count
    log_variance <- ligand_variance / ligand$center^2 +
      receptor_variance / receptor$center^2 +
      2 * cross_covariance / (ligand$center * receptor$center)
    if (!is.finite(log_variance) || log_variance <= min_log_variance) {
      output$fallback_reason[row] <- "degenerate_delta_variance"
      next
    }

    observed_ligand <- observed_expression$average[ligand_gene, sender]
    observed_receptor <- observed_expression$average[receptor_gene, receiver]
    if (!(observed_ligand > 0 && observed_receptor > 0)) {
      output$fallback_reason[row] <- "zero_observed_component"
      next
    }
    z <- (log(observed_ligand) + log(observed_receptor) -
          log(ligand$center) - log(receptor$center)) / sqrt(log_variance)
    output$pvalue[row] <- stats::pnorm(z, lower.tail = FALSE)
    output$analytic_z[row] <- z
    output$analytic_log_sd[row] <- sqrt(log_variance)
    output$inference_mode[row] <- "experimental_analytic"
    output$fallback_reason[row] <- NA_character_
  }

  structure(
    list(
      group_pairs = output,
      parameters = list(
        calibration_status = "experimental_unvalidated",
        null = "one homogeneous finite-population label permutation",
        method = "joint triMean influence functions plus log-product delta method",
        Kh = Kh, n = n, min_group_size = min_group_size,
        density_window = density_window,
        max_tie_fraction = max_tie_fraction,
        max_zero_fraction = max_zero_fraction,
        max_abs_influence_correlation = max_abs_influence_correlation,
        min_log_variance = min_log_variance,
        max_influence_entries = max_influence_entries
      ),
      diagnostics = list(
        active_records = nrow(output),
        analytic_records = sum(output$inference_mode == "experimental_analytic"),
        fallback_records = sum(
          output$inference_mode == "permutation_fallback_required"
        ),
        fallback_reasons = sort(table(output$fallback_reason), decreasing = TRUE),
        candidate_genes = length(candidate_genes),
        influence_workspace_entries = workspace_entries
      )
    ),
    class = "SpatialESSExperimentalAnalytic"
  )
}
