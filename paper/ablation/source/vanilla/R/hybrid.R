#' Experimental analytic/permutation hybrid dispatcher
#'
#' The analytic approximation is used only as a one-way non-significance
#' screen. Eligible records with analytic p-values above
#' `alpha + screening_margin` are assigned p = 1. Every other active record is
#' evaluated by the block-preserving fixed-permutation reference. Consequently,
#' every record that can be declared significant has a permutation-confirmed
#' p-value. This function is still experimental: screening can reduce power if
#' an analytic p-value incorrectly excludes a permutation-significant record.
#'
#' @param prepared A result from [prepare_sparse_trimean()].
#' @param group Fixed target-group labels.
#' @param group_support Spatially supported group pairs.
#' @param components Resolved CellChat LR components.
#' @param blocks Spatial exchangeability blocks.
#' @param nperm Number of fallback permutations.
#' @param seed First fallback permutation seed.
#' @param Kh,n Positive Hill-function parameters.
#' @param alpha Significance threshold. Analytic screening can only exclude
#'   records above `alpha + screening_margin`.
#' @param screening_margin Non-negative margin added to `alpha` for the
#'   conservative analytic non-significance screen.
#' @param finite_correction Whether fallback uses the plus-one correction.
#' @param tail Empirical upper-tail convention for fallback.
#' @param ... Conservative analytic guard parameters passed to the stratified
#'   analytic engine.
#' @export
experimental_hybrid_cellchat_group_support <- function(
    prepared, group, group_support, components, blocks,
    nperm = 999L, seed = 1L, Kh = 0.5, n = 1,
    alpha = 0.05, screening_margin = 0.01,
    finite_correction = TRUE,
    tail = c("greater_equal", "strict_greater"), ...) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(screening_margin) != 1L || !is.finite(screening_margin) ||
      screening_margin < 0 || alpha + screening_margin >= 1) {
    stop(
      "alpha must lie in (0, 1), and alpha + screening_margin must lie in (0, 1).",
      call. = FALSE
    )
  }
  tail <- match.arg(tail)
  analytic <- experimental_stratified_analytic_cellchat_group_support(
    prepared, group, group_support, components, blocks, Kh = Kh, n = n, ...
  )
  output <- analytic$group_pairs
  output$analytic_pvalue <- output$pvalue
  output$permutation_confirmed <- rep(FALSE, nrow(output))
  output$n_reject <- rep(NA_integer_, nrow(output))
  if (!nrow(output)) {
    analytic$parameters$mode <- "experimental_hybrid"
    analytic$parameters$permutation_called <- FALSE
    analytic$parameters$screening_margin <- screening_margin
    analytic$parameters$screening_threshold <- alpha + screening_margin
    analytic$parameters$screening_rule <-
      "analytic_only_excludes; every_p_at_or_below_alpha_is_permutation_confirmed"
    analytic$parameters$alpha <- alpha
    analytic$diagnostics$hybrid_screened_nonsignificant_records <- 0L
    analytic$diagnostics$hybrid_permutation_records <- 0L
    class(analytic) <- "SpatialESSExperimentalHybrid"
    return(analytic)
  }
  analytic_eligible <- output$inference_mode ==
    "experimental_stratified_analytic"
  analytic_screened <- analytic_eligible &
    output$analytic_pvalue > alpha + screening_margin
  needs_permutation <- !analytic_screened

  output$pvalue[analytic_screened] <- 1
  output$inference_mode[analytic_screened] <-
    "experimental_hybrid_screened_nonsignificant"
  if (any(needs_permutation)) {
    selected_records <- output[
      needs_permutation, c("support_pair_index", "lr_index"), drop = FALSE
    ]
    permutation <- permutation_cellchat_records(
      prepared, group, group_support, components, blocks, selected_records,
      nperm = nperm, seed = seed, Kh = Kh, n = n,
      finite_correction = finite_correction, tail = tail,
      verbose = FALSE
    )
    permutation_key <- paste(
      permutation$group_pairs$lr_index,
      permutation$group_pairs$support_pair_index, sep = "\r"
    )
    output_key <- paste(
      output$lr_index[needs_permutation],
      output$support_pair_index[needs_permutation], sep = "\r"
    )
    match_index <- match(output_key, permutation_key)
    if (anyNA(match_index)) {
      stop("Permutation fallback did not return every selected record.",
           call. = FALSE)
    }
    fallback_index <- match_index
    output$pvalue[needs_permutation] <-
      permutation$group_pairs$pvalue[fallback_index]
    output$n_reject[needs_permutation] <-
      permutation$group_pairs$n_reject[fallback_index]
    output$permutation_confirmed[needs_permutation] <- TRUE
    output$inference_mode[needs_permutation] <- ifelse(
      analytic_eligible[needs_permutation],
      "experimental_hybrid_permutation_candidate",
      "experimental_hybrid_permutation_fallback"
    )
  }
  analytic$group_pairs <- output
  analytic$parameters$mode <- "experimental_hybrid"
  analytic$parameters$permutation_called <- any(needs_permutation)
  analytic$parameters$permutation_scope <- if (any(needs_permutation)) {
    "selected_candidate_and_fallback_records"
  } else {
    "none"
  }
  analytic$parameters$nperm <- if (any(needs_permutation)) nperm else 0L
  analytic$parameters$seed <- if (any(needs_permutation)) seed else NA_integer_
  analytic$parameters$finite_correction <- isTRUE(finite_correction)
  analytic$parameters$tail <- tail
  analytic$parameters$alpha <- alpha
  analytic$parameters$screening_margin <- screening_margin
  analytic$parameters$screening_threshold <- alpha + screening_margin
  analytic$parameters$screening_rule <-
    "analytic_only_excludes; every_p_at_or_below_alpha_is_permutation_confirmed"
  analytic$diagnostics$hybrid_screened_nonsignificant_records <-
    sum(analytic_screened)
  analytic$diagnostics$hybrid_permutation_records <- sum(needs_permutation)
  if (any(output$pvalue <= alpha & !output$permutation_confirmed)) {
    stop(
      "Internal error: a significant hybrid result lacks permutation confirmation.",
      call. = FALSE
    )
  }
  class(analytic) <- "SpatialESSExperimentalHybrid"
  analytic
}
