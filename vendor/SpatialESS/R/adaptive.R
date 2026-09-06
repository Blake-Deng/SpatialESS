#' Experimental adaptive block-permutation inference
#'
#' This threshold-level prototype runs block-preserving permutations in batches
#' and stops each active record when a Bonferroni-adjusted Clopper-Pearson
#' interval lies wholly above or at/below `alpha`. It is intended for single
#' threshold screening. Complete fixed-permutation inference remains required
#' for BH/FDR, pathway/network aggregation and publication p-values.
#'
#' Each batch currently evaluates the complete supported score matrix because
#' the lower-level CellChat-compatible scorer has not yet been record-sliced.
#' The API nevertheless records per-record stopping counts and provides the
#' correct place for future active-record batching.
#'
#' @param prepared A prepared sparse triMean kernel.
#' @param group Fixed target-group labels.
#' @param group_support Spatially supported group pairs.
#' @param components Resolved CellChat LR components.
#' @param blocks Spatial exchangeability blocks.
#' @param min_permutations Minimum permutations before checking intervals.
#' @param batch_size Number of permutations added at each check.
#' @param max_permutations Hard upper bound on permutations.
#' @param seed First non-negative permutation seed.
#' @param alpha Single decision threshold.
#' @param confidence_level Simultaneous confidence target over batch checks.
#' @param Kh,n Positive Hill-function parameters.
#' @param tail Empirical upper-tail convention.
#' @param max_score_entries Maximum dense score-workspace entries.
#' @export
experimental_adaptive_permutation_cellchat_group_support <- function(
    prepared, group, group_support, components, blocks,
    min_permutations = 100L, batch_size = 100L,
    max_permutations = 1000L, seed = 1L, alpha = 0.05,
    confidence_level = 0.99, Kh = 0.5, n = 1,
    tail = c("greater_equal", "strict_greater"),
    max_score_entries = 1e8) {
  min_permutations <- as.integer(min_permutations)
  batch_size <- as.integer(batch_size)
  max_permutations <- as.integer(max_permutations)
  seed <- as.integer(seed)
  tail <- match.arg(tail)
  if (length(min_permutations) != 1L || is.na(min_permutations) ||
      min_permutations < 1L || length(batch_size) != 1L ||
      is.na(batch_size) || batch_size < 1L ||
      length(max_permutations) != 1L || is.na(max_permutations) ||
      max_permutations < min_permutations || length(seed) != 1L ||
      is.na(seed) || as.double(seed) + max_permutations - 1 >
        .Machine$integer.max) {
    stop("Permutation limits and seed must define a valid integer sequence.",
         call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(confidence_level) != 1L || !is.finite(confidence_level) ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("alpha and confidence_level must lie strictly between zero and one.",
         call. = FALSE)
  }
  max_batches <- ceiling(max_permutations / batch_size) + 1L
  batch_confidence <- 1 - (1 - confidence_level) / max_batches
  reject <- NULL
  keys <- NULL
  output_template <- NULL
  active <- NULL
  permutations_used <- NULL
  stop_reason <- NULL
  lower <- NULL
  upper <- NULL
  permutations_done <- 0L
  batch_index <- 0L
  record_permutation_evaluations <- 0

  while (permutations_done < max_permutations &&
         (is.null(active) || any(active))) {
    batch_index <- batch_index + 1L
    remaining <- max_permutations - permutations_done
    batch_n <- min(batch_size, remaining)
    if (is.null(keys)) {
      result <- permutation_cellchat_group_support(
        prepared, group, group_support, components, blocks,
        nperm = batch_n, seed = seed + permutations_done,
        Kh = Kh, n = n, finite_correction = FALSE, tail = tail,
        max_score_entries = max_score_entries, verbose = FALSE
      )
    } else {
      selected_records <- output_template[
        active, c("support_pair_index", "lr_index"), drop = FALSE
      ]
      result <- permutation_cellchat_records(
        prepared, group, group_support, components, blocks, selected_records,
        nperm = batch_n, seed = seed + permutations_done,
        Kh = Kh, n = n, finite_correction = FALSE, tail = tail,
        max_records = max_score_entries, verbose = FALSE
      )
    }
    current <- result$group_pairs
    record_permutation_evaluations <- record_permutation_evaluations +
      as.double(nrow(current)) * batch_n
    current_keys <- paste(current$lr_index, current$support_pair_index, sep = "\r")
    if (is.null(keys)) {
      keys <- current_keys
      output_template <- current
      output_template$pvalue <- NA_real_
      output_template$n_reject <- NA_integer_
      output_template$permutations_used <- NA_integer_
      output_template$pvalue_lower <- NA_real_
      output_template$pvalue_upper <- NA_real_
      output_template$decision <- NA_character_
      output_template$stop_reason <- NA_character_
      reject <- integer(length(keys))
      active <- rep(TRUE, length(keys))
      permutations_used <- integer(length(keys))
      stop_reason <- rep(NA_character_, length(keys))
      lower <- rep(NA_real_, length(keys))
      upper <- rep(NA_real_, length(keys))
    }
    active_index <- which(active)
    match_index <- match(keys[active_index], current_keys)
    if (anyNA(match_index) || nrow(current) != length(active_index)) {
      stop("Permutation batch did not return every requested active record.",
           call. = FALSE)
    }
    reject[active_index] <- reject[active_index] +
      current$n_reject[match_index]
    permutations_done <- permutations_done + batch_n
    inspect <- permutations_done >= min_permutations
    if (inspect) {
      check_index <- which(active)
      for (index in check_index) {
        interval <- stats::binom.test(
          reject[[index]], permutations_done,
          conf.level = batch_confidence
        )$conf.int
        lower[[index]] <- (interval[[1L]] * permutations_done + 1) /
          (permutations_done + 1)
        upper[[index]] <- (interval[[2L]] * permutations_done + 1) /
          (permutations_done + 1)
        if (upper[[index]] <= alpha) {
          active[[index]] <- FALSE
          stop_reason[[index]] <- "pvalue_upper_at_or_below_alpha"
          permutations_used[[index]] <- permutations_done
        } else if (lower[[index]] > alpha) {
          active[[index]] <- FALSE
          stop_reason[[index]] <- "pvalue_lower_above_alpha"
          permutations_used[[index]] <- permutations_done
        }
      }
    }
    if (permutations_done >= max_permutations) {
      still_active <- which(active)
      stop_reason[still_active] <- "max_permutations"
      permutations_used[still_active] <- permutations_done
      active[still_active] <- FALSE
    }
  }
  if (is.null(output_template)) {
    stop("No active records were emitted by the permutation reference.",
         call. = FALSE)
  }
  used <- pmax(1L, permutations_used)
  output_template$n_reject <- reject
  output_template$permutations_used <- used
  output_template$pvalue <- (reject + 1) / (used + 1)
  output_template$pvalue_lower <- lower
  output_template$pvalue_upper <- upper
  output_template$decision <- ifelse(
    upper <= alpha, "significant", "not_significant"
  )
  output_template$stop_reason <- stop_reason
  structure(
    list(
      group_pairs = output_template,
      parameters = list(
        mode = "experimental_adaptive_permutation",
        alpha = alpha, confidence_level = confidence_level,
        batch_confidence = batch_confidence,
        min_permutations = min_permutations, batch_size = batch_size,
        max_permutations = max_permutations, seed = seed,
        Kh = Kh, n = n, tail = tail,
        multiple_testing = "not calibrated for BH/FDR; use fixed permutation"
      ),
      diagnostics = list(
        batches = batch_index, total_permutations = permutations_done,
        stopped_early = permutations_done < max_permutations,
        active_at_maximum = sum(stop_reason == "max_permutations"),
        record_permutation_evaluations = record_permutation_evaluations,
        full_matrix_record_permutation_evaluations =
          as.double(length(keys)) * permutations_done,
        pvalue_upper_mean = mean(upper), pvalue_lower_mean = mean(lower)
      )
    ),
    class = "SpatialESSExperimentalAdaptivePermutation"
  )
}
