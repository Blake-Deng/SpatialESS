#' Guarded hybrid spatial significance
#'
#' Uses the validated one-way analytic screen only to exclude clearly
#' non-significant records. Every record at or below `alpha` is confirmed by
#' fixed, block-preserving permutation. Direct analytic and adaptive p-values
#' are intentionally not exposed by this production wrapper.
#'
#' @inheritParams experimental_hybrid_cellchat_group_support
#' @export
guarded_hybrid_cellchat_group_support <- function(
    prepared, group, group_support, components, blocks,
    nperm = 999L, seed = 1L, Kh = 0.5, n = 1,
    alpha = 0.05, screening_margin = 0.01,
    finite_correction = TRUE,
    tail = c("greater_equal", "strict_greater"), ...) {
  result <- experimental_hybrid_cellchat_group_support(
    prepared = prepared, group = group, group_support = group_support,
    components = components, blocks = blocks,
    nperm = nperm, seed = seed, Kh = Kh, n = n,
    alpha = alpha, screening_margin = screening_margin,
    finite_correction = finite_correction, tail = tail, ...
  )
  result$parameters$mode <- "guarded_hybrid"
  result$parameters$calibration_status <-
    "validated_one_way_screen_with_fixed_permutation_confirmation"
  result$parameters$publication_boundary <- paste(
    "Every discovery is fixed-permutation confirmed; direct analytic and",
    "adaptive p-values are not publication inference."
  )
  class(result) <- c("SpatialESSGuardedHybrid", class(result))
  result
}

.spatialess_validate_integrated_inputs <- function(
    expression, coordinates, group, sample_id, compartment) {
  if (!inherits(expression, "Matrix")) {
    expression <- Matrix::Matrix(expression, sparse = TRUE)
  }
  expression <- methods::as(expression, "dgCMatrix")
  if (is.null(rownames(expression)) || is.null(colnames(expression)) ||
      anyNA(rownames(expression)) || anyNA(colnames(expression)) ||
      anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression))) {
    stop("expression requires unique gene and cell names.", call. = FALSE)
  }
  if (any(!is.finite(expression@x)) || any(expression@x < 0)) {
    stop("expression must contain finite non-negative values.", call. = FALSE)
  }
  coordinates <- as.matrix(coordinates)
  storage.mode(coordinates) <- "double"
  if (nrow(coordinates) != ncol(expression) ||
      !ncol(coordinates) %in% c(2L, 3L) || any(!is.finite(coordinates))) {
    stop("coordinates must provide two or three finite columns per cell.",
         call. = FALSE)
  }
  if (is.null(rownames(coordinates))) {
    rownames(coordinates) <- colnames(expression)
  }
  if (!identical(rownames(coordinates), colnames(expression))) {
    stop("coordinate and expression cell order differs.", call. = FALSE)
  }
  info <- .spatialess_group_info(group)
  if (length(info$code) != ncol(expression)) {
    stop("group must provide one label per cell.", call. = FALSE)
  }
  normalize_label <- function(value, default, name) {
    if (is.null(value)) value <- rep(default, ncol(expression))
    value <- as.character(value)
    if (length(value) != ncol(expression) || anyNA(value) || any(!nzchar(value))) {
      stop(sprintf("%s must provide one non-empty label per cell.", name),
           call. = FALSE)
    }
    value
  }
  list(
    expression = expression, coordinates = coordinates,
    group = factor(info$code, levels = seq_len(info$count),
                   labels = info$levels[seq_len(info$count)]),
    sample_id = normalize_label(sample_id, "sample_1", "sample_id"),
    compartment = normalize_label(compartment, "all", "compartment")
  )
}

.spatialess_run_mechanism <- function(
    mechanism, lr_index, graph, expression, group, components, blocks,
    inference, nperm, seed, Kh, n, alpha, screening_margin,
    finite_correction, tail, max_group_pairs, max_score_entries,
    max_records, analytic_guard, verbose) {
  if (!length(lr_index)) return(NULL)
  mechanism_components <- components
  component_names <- c("ligand", "receptor", "co_a", "co_i",
                       "agonist", "antagonist")
  for (name in component_names) {
    mechanism_components[[name]] <- components[[name]][lr_index]
  }
  mechanism_components$has_agonist <- components$has_agonist[lr_index]
  mechanism_components$has_antagonist <- components$has_antagonist[lr_index]
  mechanism_components$lr <- components$lr[lr_index, , drop = FALSE]
  mechanism_components$genes <- unique(unlist(
    mechanism_components[component_names], recursive = TRUE,
    use.names = FALSE
  ))
  class(mechanism_components) <- "SpatialESSLRComponents"

  support <- build_group_support_csr(
    graph, group, max_group_pairs = max_group_pairs
  )
  prepared <- prepare_sparse_trimean(
    expression, genes = mechanism_components$genes, normalize = TRUE
  )
  if (inference == "permutation") {
    result <- permutation_cellchat_group_support(
      prepared, group, support, mechanism_components, blocks,
      nperm = nperm, seed = seed, Kh = Kh, n = n,
      finite_correction = finite_correction,
      retain_null_moments = FALSE, tail = tail,
      max_score_entries = max_score_entries, verbose = verbose
    )
    result$group_pairs$inference_mode <- "fixed_spatial_permutation"
    result$group_pairs$permutation_confirmed <- TRUE
    result$parameters$mode <- "fixed_spatial_permutation"
  } else {
    arguments <- c(list(
      prepared = prepared, group = group, group_support = support,
      components = mechanism_components, blocks = blocks,
      nperm = nperm, seed = seed, Kh = Kh, n = n,
      alpha = alpha, screening_margin = screening_margin,
      finite_correction = finite_correction, tail = tail
    ), analytic_guard)
    result <- do.call(guarded_hybrid_cellchat_group_support, arguments)
  }
  records <- result$group_pairs
  records$global_lr_index <- lr_index[records$lr_index]
  records$mechanism <- mechanism
  records$graph_radius <- graph$parameters$radius
  records$graph_edges <- length(graph$neighbors)
  if (nrow(records) > max_records) {
    stop("Mechanism result exceeds max_records.", call. = FALSE)
  }
  list(
    records = records, graph = graph, support = support,
    components = mechanism_components, inference = result,
    diagnostics = data.frame(
      mechanism = mechanism, lr = length(lr_index),
      referenced_genes = length(mechanism_components$genes),
      graph_radius = graph$parameters$radius,
      graph_edges = length(graph$neighbors),
      supported_group_pairs = nrow(support$group_pairs),
      active_records = nrow(records),
      significant_records = sum(records$pvalue <= alpha),
      permutation_confirmed_significant = sum(
        records$pvalue <= alpha & records$permutation_confirmed
      ), stringsAsFactors = FALSE
    )
  )
}

#' Integrated single-sample SpatialESS workflow
#'
#' Builds separate compact CSR graphs for contact and diffusible signaling,
#' resolves only LR-referenced genes, computes exact sparse type-7 triMeans and
#' CellChat-compatible molecular scores, and performs either fixed spatial
#' block permutation or the validated guarded hybrid significance procedure.
#'
#' @param expression Non-negative gene-by-cell matrix.
#' @param coordinates Cell-by-2/3 spatial coordinates aligned to expression.
#' @param group Cell state/type labels.
#' @param lr CellChat-style ligand-receptor table.
#' @param complex,cofactor CellChat component tables.
#' @param sample_id,compartment Optional cell-aligned exchangeability labels.
#' @param contact_radius,diffusion_radius Positive mechanism-specific radii.
#' @param block_size Positive spatial permutation block size.
#' @param inference Validated significance mode.
#' @param nperm,seed Fixed permutation count and first seed.
#' @param Kh,n CellChat Hill parameters.
#' @param alpha,screening_margin Guarded hybrid decision boundary.
#' @param finite_correction,tail Empirical p-value convention.
#' @param contact_annotation LR annotation identifying contact signaling.
#' @param same_compartment Whether graph edges must stay in compartments.
#' @param graph_weight,graph_scale CSR edge-weight configuration.
#' @param max_edges,max_group_pairs,max_score_entries,max_records Memory guards.
#' @param analytic_guard Named list of calibrated analytic guard overrides.
#' @param verbose Whether to report permutation progress.
#' @export
spatialess <- function(
    expression, coordinates, group, lr, complex, cofactor,
    sample_id = NULL, compartment = NULL,
    contact_radius = 15, diffusion_radius = 35, block_size = 100,
    inference = c("permutation", "guarded_hybrid"),
    nperm = 999L, seed = 1L, Kh = 0.5, n = 1,
    alpha = 0.05, screening_margin = 0.01,
    finite_correction = TRUE,
    tail = c("greater_equal", "strict_greater"),
    contact_annotation = "Cell-Cell Contact",
    same_compartment = FALSE,
    graph_weight = c("binary", "gaussian", "exponential"),
    graph_scale = diffusion_radius / 2,
    max_edges = 5e8, max_group_pairs = 1e7,
    max_score_entries = 1e8, max_records = 1e8,
    analytic_guard = list(), verbose = interactive()) {
  inference <- match.arg(inference)
  tail <- match.arg(tail)
  graph_weight <- match.arg(graph_weight)
  if (!is.list(analytic_guard) ||
      (length(analytic_guard) && is.null(names(analytic_guard)))) {
    stop("analytic_guard must be a named list.", call. = FALSE)
  }
  radius <- c(contact = contact_radius, diffusion = diffusion_radius)
  if (length(contact_radius) != 1L || length(diffusion_radius) != 1L ||
      length(block_size) != 1L || any(!is.finite(c(radius, block_size))) ||
      any(radius <= 0) || block_size <= 0) {
    stop("Radii and block_size must be positive finite scalars.", call. = FALSE)
  }
  input <- .spatialess_validate_integrated_inputs(
    expression, coordinates, group, sample_id, compartment
  )
  lr <- as.data.frame(lr, stringsAsFactors = FALSE)
  if (!all(c("ligand", "receptor") %in% colnames(lr))) {
    stop("lr is missing ligand or receptor columns.", call. = FALSE)
  }
  annotation <- if ("annotation" %in% colnames(lr)) {
    as.character(lr$annotation)
  } else {
    rep("Secreted Signaling", nrow(lr))
  }
  components <- prepare_cellchat_lr_components(
    lr, rownames(input$expression), complex, cofactor
  )
  blocks <- build_spatial_blocks(
    input$coordinates, block_size = block_size,
    sample_id = input$sample_id, compartment = input$compartment
  )
  graph_arguments <- list(
    coords = input$coordinates, sample_id = input$sample_id,
    compartment = input$compartment,
    same_compartment = isTRUE(same_compartment), weight = graph_weight,
    store_distance = FALSE, max_edges = max_edges
  )
  contact_graph <- do.call(build_radius_graph_csr, c(
    graph_arguments,
    list(radius = contact_radius,
         scale = min(graph_scale, contact_radius))
  ))
  diffusion_graph <- if (identical(contact_radius, diffusion_radius)) {
    contact_graph
  } else {
    do.call(build_radius_graph_csr, c(
      graph_arguments,
      list(radius = diffusion_radius, scale = graph_scale)
    ))
  }
  contact_index <- which(annotation == contact_annotation)
  diffusion_index <- which(annotation != contact_annotation)
  run <- list(
    contact = .spatialess_run_mechanism(
      "contact", contact_index, contact_graph, input$expression, input$group,
      components, blocks, inference, nperm, seed, Kh, n, alpha,
      screening_margin, finite_correction, tail, max_group_pairs,
      max_score_entries, max_records, analytic_guard, verbose
    ),
    diffusion = .spatialess_run_mechanism(
      "diffusion", diffusion_index, diffusion_graph, input$expression,
      input$group, components, blocks, inference, nperm, seed, Kh, n, alpha,
      screening_margin, finite_correction, tail, max_group_pairs,
      max_score_entries, max_records, analytic_guard, verbose
    )
  )
  run <- run[!vapply(run, is.null, logical(1))]
  records <- do.call(rbind, lapply(run, `[[`, "records"))
  rownames(records) <- NULL
  lr_columns <- setdiff(c("interaction_name", "pathway_name", "annotation"),
                        colnames(records))
  for (name in lr_columns) {
    records[[name]] <- lr[[name]][records$global_lr_index]
  }
  if (any(records$pvalue <= alpha & !records$permutation_confirmed)) {
    stop("Internal error: a significant result lacks permutation confirmation.",
         call. = FALSE)
  }
  diagnostics <- do.call(rbind, lapply(run, `[[`, "diagnostics"))
  structure(
    list(
      records = records, diagnostics = diagnostics,
      graphs = lapply(run, `[[`, "graph"),
      group_support = lapply(run, `[[`, "support"),
      inference = lapply(run, `[[`, "inference"),
      blocks = blocks, lr = lr,
      parameters = list(
        engine = "integrated_spatialess_csr",
        inference = inference, contact_radius = contact_radius,
        diffusion_radius = diffusion_radius, block_size = block_size,
        nperm = as.integer(nperm), seed = as.integer(seed),
        Kh = Kh, n = n, alpha = alpha,
        screening_margin = screening_margin,
        finite_correction = isTRUE(finite_correction), tail = tail,
        null = paste(
          "complete cell profiles permuted within sample x compartment x",
          "spatial block"
        ),
        significance_boundary =
          "every p <= alpha is fixed-permutation confirmed"
      )
    ),
    class = "SpatialESSIntegrated"
  )
}
