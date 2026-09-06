.validate_graph_inputs <- function(coords, sample_id = NULL, compartment = NULL) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  if (!nrow(coords) || ncol(coords) < 2L || ncol(coords) > 3L) {
    stop("coords must be a non-empty numeric matrix with two or three columns.", call. = FALSE)
  }
  if (any(!is.finite(coords))) {
    stop("coords contains non-finite values.", call. = FALSE)
  }
  n <- nrow(coords)
  if (is.null(rownames(coords))) rownames(coords) <- sprintf("cell_%d", seq_len(n))
  if (anyNA(rownames(coords)) || anyDuplicated(rownames(coords))) {
    stop("coords row names must be unique, non-missing cell identifiers.", call. = FALSE)
  }

  if (is.null(sample_id)) sample_id <- rep.int("sample_1", n)
  if (length(sample_id) != n || anyNA(sample_id)) {
    stop("sample_id must have one non-missing value per cell.", call. = FALSE)
  }
  sample_int <- as.integer(factor(sample_id, levels = unique(as.character(sample_id))))

  if (is.null(compartment)) compartment <- rep.int("all", n)
  if (length(compartment) != n || anyNA(compartment)) {
    stop("compartment must have one non-missing value per cell.", call. = FALSE)
  }
  compartment_int <- as.integer(factor(compartment, levels = unique(as.character(compartment))))

  list(
    coords = coords,
    cell_id = rownames(coords),
    sample_id = as.character(sample_id),
    sample_int = sample_int,
    compartment = as.character(compartment),
    compartment_int = compartment_int
  )
}

.edge_weights <- function(distance, weight, scale) {
  switch(
    weight,
    binary = rep.int(1, length(distance)),
    gaussian = exp(-(distance^2) / (2 * scale^2)),
    exponential = exp(-distance / scale),
    stop("Unsupported weight function.", call. = FALSE)
  )
}

.new_spatial_graph <- function(edges, input, graph_type, parameters) {
  if (nrow(edges)) {
    edges$sender_id <- input$cell_id[edges$sender]
    edges$receiver_id <- input$cell_id[edges$receiver]
    edges$sample_id <- input$sample_id[edges$sender]
    edges <- edges[, c("sender", "receiver", "sender_id", "receiver_id",
                       "sample_id", "distance", "weight"), drop = FALSE]
  } else {
    edges <- data.frame(
      sender = integer(), receiver = integer(), sender_id = character(),
      receiver_id = character(), sample_id = character(), distance = double(),
      weight = double(), stringsAsFactors = FALSE
    )
  }
  structure(
    list(
      edges = edges,
      cells = data.frame(
        cell_index = seq_along(input$cell_id),
        cell_id = input$cell_id,
        sample_id = input$sample_id,
        compartment = input$compartment,
        stringsAsFactors = FALSE
      ),
      graph_type = graph_type,
      directed = TRUE,
      parameters = parameters
    ),
    class = "SpatialESSGraph"
  )
}

#' Build a deterministic sparse radius graph
#'
#' Candidate pairs are generated with a uniform spatial hash. Pairs are never
#' formed across samples, and an optional compartment barrier can be enforced.
#'
#' @param coords A finite cell-by-coordinate matrix.
#' @param radius Positive neighborhood radius.
#' @param sample_id,compartment Optional cell-aligned boundary labels.
#' @param same_compartment Whether edges must remain within compartments.
#' @param weight Spatial edge weighting rule.
#' @param scale Positive decay scale for weighted graphs.
#' @param include_self Whether to retain self edges.
#' @export
build_radius_graph <- function(coords, radius, sample_id = NULL,
                               compartment = NULL,
                               same_compartment = FALSE,
                               weight = c("binary", "gaussian", "exponential"),
                               scale = radius / 2,
                               include_self = FALSE) {
  weight <- match.arg(weight)
  if (length(radius) != 1L || !is.finite(radius) || radius <= 0) {
    stop("radius must be one positive finite number.", call. = FALSE)
  }
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("scale must be one positive finite number.", call. = FALSE)
  }
  input <- .validate_graph_inputs(coords, sample_id, compartment)
  edges <- as.data.frame(radius_graph_cpp(
    input$coords, input$sample_int, input$compartment_int,
    radius, same_compartment, include_self, weight, scale
  ), stringsAsFactors = FALSE)
  edges$sender <- as.integer(edges$sender)
  edges$receiver <- as.integer(edges$receiver)
  .new_spatial_graph(
    edges, input, "radius",
    list(radius = radius, weight = weight, scale = scale,
         same_compartment = same_compartment, include_self = include_self,
         engine = "spatial_hash_cpp")
  )
}

#' Brute-force radius graph for correctness testing
#'
#' @param coords A finite cell-by-coordinate matrix.
#' @param radius Positive neighborhood radius.
#' @param sample_id,compartment Optional cell-aligned boundary labels.
#' @param same_compartment Whether edges must remain within compartments.
#' @param weight Spatial edge weighting rule.
#' @param scale Positive decay scale for weighted graphs.
#' @param include_self Whether to retain self edges.
#' @param max_cells Maximum cells accepted by the quadratic oracle.
#' @export
build_radius_graph_bruteforce <- function(coords, radius, sample_id = NULL,
                                          compartment = NULL,
                                          same_compartment = FALSE,
                                          weight = c("binary", "gaussian", "exponential"),
                                          scale = radius / 2,
                                          include_self = FALSE,
                                          max_cells = 5000L) {
  weight <- match.arg(weight)
  input <- .validate_graph_inputs(coords, sample_id, compartment)
  n <- nrow(input$coords)
  if (n > max_cells) {
    stop("The brute-force oracle is limited to max_cells.", call. = FALSE)
  }
  d <- as.matrix(stats::dist(input$coords, upper = TRUE, diag = TRUE))
  eligible <- outer(input$sample_int, input$sample_int, `==`) & d <= radius
  if (same_compartment) {
    eligible <- eligible & outer(input$compartment_int, input$compartment_int, `==`)
  }
  if (!include_self) diag(eligible) <- FALSE
  idx <- which(eligible, arr.ind = TRUE)
  if (nrow(idx)) {
    ord <- order(idx[, 1L], idx[, 2L])
    idx <- idx[ord, , drop = FALSE]
    distance <- d[idx]
    edges <- data.frame(
      sender = idx[, 1L], receiver = idx[, 2L], distance = distance,
      weight = .edge_weights(distance, weight, scale)
    )
  } else {
    edges <- data.frame(sender = integer(), receiver = integer(),
                        distance = double(), weight = double())
  }
  .new_spatial_graph(
    edges, input, "radius",
    list(radius = radius, weight = weight, scale = scale,
         same_compartment = same_compartment, include_self = include_self,
         engine = "bruteforce_R")
  )
}

#' Build a directed k-nearest-neighbor graph within each sample
#'
#' @param coords A finite cell-by-coordinate matrix.
#' @param k Positive number of directed nearest neighbors.
#' @param sample_id,compartment Optional cell-aligned boundary labels.
#' @param same_compartment Whether edges must remain within compartments.
#' @param max_radius Optional maximum neighbor distance.
#' @param weight Spatial edge weighting rule.
#' @param scale Positive decay scale for weighted graphs.
#' @export
build_knn_graph <- function(coords, k = 10L, sample_id = NULL,
                            compartment = NULL,
                            same_compartment = FALSE,
                            max_radius = Inf,
                            weight = c("binary", "gaussian", "exponential"),
                            scale = NULL) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for kNN graph construction.", call. = FALSE)
  }
  weight <- match.arg(weight)
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) stop("k must be positive.", call. = FALSE)
  input <- .validate_graph_inputs(coords, sample_id, compartment)
  edge_parts <- vector("list", max(input$sample_int))
  for (s in seq_along(edge_parts)) {
    cells <- which(input$sample_int == s)
    if (length(cells) <= 1L) next
    k_use <- min(k + 1L, length(cells))
    nn <- RANN::nn2(input$coords[cells, , drop = FALSE], k = k_use)
    from <- rep(cells, each = k_use)
    to <- cells[as.vector(t(nn$nn.idx))]
    distance <- as.vector(t(nn$nn.dists))
    keep <- from != to & distance <= max_radius
    if (same_compartment) keep <- keep & input$compartment_int[from] == input$compartment_int[to]
    edge_parts[[s]] <- data.frame(sender = from[keep], receiver = to[keep],
                                  distance = distance[keep])
  }
  edges <- do.call(rbind, edge_parts)
  if (is.null(edges)) {
    edges <- data.frame(sender = integer(), receiver = integer(), distance = double())
  } else {
    edges <- unique(edges)
    edges <- edges[order(edges$sender, edges$receiver), , drop = FALSE]
    rownames(edges) <- NULL
  }
  if (is.null(scale)) {
    positive <- edges$distance[edges$distance > 0]
    scale <- if (length(positive)) stats::median(positive) else 1
  }
  if (!is.finite(scale) || scale <= 0) stop("scale must be positive.", call. = FALSE)
  edges$weight <- .edge_weights(edges$distance, weight, scale)
  .new_spatial_graph(
    edges, input, "knn",
    list(k = k, max_radius = max_radius, weight = weight, scale = scale,
         same_compartment = same_compartment, engine = "RANN")
  )
}

#' Validate structural invariants of a SpatialESS graph
#'
#' @param graph A SpatialESS graph.
#' @param tolerance Non-negative numerical comparison tolerance.
#' @export
validate_spatial_graph <- function(graph, tolerance = 1e-12) {
  stopifnot(inherits(graph, "SpatialESSGraph"))
  e <- graph$edges
  n <- nrow(graph$cells)
  checks <- c(
    indices_in_range = !nrow(e) || all(e$sender >= 1L & e$sender <= n &
                                      e$receiver >= 1L & e$receiver <= n),
    no_duplicate_edges = !anyDuplicated(paste(e$sender, e$receiver, sep = ":")),
    no_cross_sample_edges = !nrow(e) || all(graph$cells$sample_id[e$sender] ==
                                            graph$cells$sample_id[e$receiver]),
    finite_nonnegative_distance = !nrow(e) || all(is.finite(e$distance) & e$distance >= 0),
    finite_nonnegative_weight = !nrow(e) || all(is.finite(e$weight) & e$weight >= -tolerance),
    canonical_order = !nrow(e) || !is.unsorted(order(e$sender, e$receiver))
  )
  structure(checks, valid = all(checks))
}

