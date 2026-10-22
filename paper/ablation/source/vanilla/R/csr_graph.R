#' Build a compact CSR radius graph
#'
#' Uses a two-pass uniform spatial hash. The first pass counts neighbors and the
#' second fills compact row offsets, neighbor indices and edge weights without
#' materializing a sender column or an edge data frame.
#'
#' @param coords A finite cell-by-coordinate matrix.
#' @param radius Positive neighborhood radius.
#' @param sample_id,compartment Optional cell-aligned boundary labels.
#' @param same_compartment Whether edges must remain within compartments.
#' @param weight Spatial edge weighting rule.
#' @param scale Positive decay scale for weighted graphs.
#' @param store_distance Whether to retain edge distances.
#' @param max_edges Maximum directed edges.
#' @export
build_radius_graph_csr <- function(coords, radius, sample_id = NULL,
                                   compartment = NULL,
                                   same_compartment = FALSE,
                                   weight = c("binary", "gaussian", "exponential"),
                                   scale = radius / 2,
                                   store_distance = FALSE,
                                   max_edges = 5e8) {
  weight <- match.arg(weight)
  if (length(radius) != 1L || !is.finite(radius) || radius <= 0) {
    stop("radius must be one positive finite number.", call. = FALSE)
  }
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("scale must be one positive finite number.", call. = FALSE)
  }
  max_edges <- as.double(max_edges)
  if (!is.finite(max_edges) || max_edges < 0) {
    stop("max_edges must be finite and non-negative.", call. = FALSE)
  }
  input <- .validate_graph_inputs(coords, sample_id, compartment)
  result <- radius_graph_csr_cpp(
    input$coords, input$sample_int, input$compartment_int,
    radius, same_compartment, weight, scale,
    isTRUE(store_distance), max_edges
  )
  structure(
    list(
      offsets = result$offsets,
      neighbors = result$neighbors,
      weights = result$weights,
      distances = result$distances,
      degree = result$degree,
      weighted_degree = result$weighted_degree,
      cells = list(
        cell_id = input$cell_id,
        sample_code = input$sample_int,
        sample_levels = unique(input$sample_id),
        compartment_code = input$compartment_int,
        compartment_levels = unique(input$compartment)
      ),
      graph_type = "radius",
      symmetric = TRUE,
      parameters = list(
        radius = radius, weight = weight, scale = scale,
        same_compartment = same_compartment,
        store_distance = isTRUE(store_distance),
        engine = "two_pass_spatial_hash_csr_cpp"
      )
    ),
    class = "SpatialESSCSRGraph"
  )
}

#' Validate a compact SpatialESS CSR graph
#'
#' @param graph A compact SpatialESS CSR graph.
#' @export
validate_spatial_csr_graph <- function(graph) {
  if (!inherits(graph, "SpatialESSCSRGraph")) {
    stop("graph must be a SpatialESSCSRGraph.", call. = FALSE)
  }
  result <- validate_spatial_csr_cpp(
    graph$offsets, graph$neighbors, graph$weights,
    graph$cells$sample_code, graph$cells$compartment_code,
    isTRUE(graph$parameters$same_compartment)
  )
  checks <- unlist(result, use.names = TRUE)
  structure(checks, valid = all(checks))
}

#' Materialize a compact graph for small-data inspection
#'
#' @param graph A compact SpatialESS CSR graph.
#' @param max_edges Maximum edges allowed for materialization.
#' @export
materialize_csr_edges <- function(graph, max_edges = 1e6) {
  if (!inherits(graph, "SpatialESSCSRGraph")) {
    stop("graph must be a SpatialESSCSRGraph.", call. = FALSE)
  }
  if (length(graph$neighbors) > max_edges) {
    stop("CSR graph exceeds max_edges; keep it compact.", call. = FALSE)
  }
  result <- as.data.frame(materialize_csr_edges_cpp(
    graph$offsets, graph$neighbors, graph$weights, graph$distances
  ), stringsAsFactors = FALSE)
  result$sender_id <- graph$cells$cell_id[result$sender]
  result$receiver_id <- graph$cells$cell_id[result$receiver]
  result
}
