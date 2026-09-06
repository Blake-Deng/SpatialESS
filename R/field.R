#' Solve a steady-state ligand field on a sparse spatial graph
#'
#' Solves (decay I + diffusion L + uptake diag(receptor)) concentration =
#' production_rate source with a matrix-free Jacobi-preconditioned conjugate
#' gradient method. The current prototype requires a symmetric radius graph.
#'
#' @param graph A symmetric SpatialESS graph.
#' @param source,receptor Non-negative cell-aligned vectors.
#' @param diffusion,decay,uptake,production_rate Non-negative field parameters; `decay` must be positive.
#' @param tolerance Positive solver tolerance.
#' @param max_iterations Positive maximum PCG iterations.
#' @export
solve_ligand_field <- function(graph, source, receptor = NULL,
                               diffusion = 1, decay = 1, uptake = 0,
                               production_rate = 1,
                               tolerance = 1e-8, max_iterations = 1000L) {
  if (!inherits(graph, "SpatialESSGraph")) {
    stop("graph must be a SpatialESSGraph.", call. = FALSE)
  }
  if (!identical(graph$graph_type, "radius")) {
    stop("The F1 solver currently requires a symmetric radius graph.", call. = FALSE)
  }
  cell_id <- graph$cells$cell_id
  align_vector <- function(x, label) {
    if (is.null(x)) return(rep.int(0, length(cell_id)))
    x_names <- names(x)
    x <- as.numeric(x)
    if (length(x) != length(cell_id)) {
      stop(label, " must have one value per graph cell.", call. = FALSE)
    }
    if (!is.null(x_names)) {
      if (!setequal(x_names, cell_id)) {
        stop(label, " names do not match graph cell identifiers.", call. = FALSE)
      }
      x <- x[match(cell_id, x_names)]
    }
    if (any(!is.finite(x)) || any(x < 0)) {
      stop(label, " must contain finite non-negative values.", call. = FALSE)
    }
    x
  }
  source <- align_vector(source, "source")
  receptor <- align_vector(receptor, "receptor")
  parameters <- c(diffusion = diffusion, decay = decay, uptake = uptake,
                  production_rate = production_rate, tolerance = tolerance)
  if (any(!is.finite(parameters)) || diffusion < 0 || decay <= 0 ||
      uptake < 0 || production_rate < 0 || tolerance <= 0) {
    stop("Invalid reaction-diffusion parameter.", call. = FALSE)
  }
  max_iterations <- as.integer(max_iterations)
  if (is.na(max_iterations) || max_iterations < 1L) {
    stop("max_iterations must be a positive integer.", call. = FALSE)
  }
  e <- graph$edges
  result <- graph_field_cg_cpp(
    e$sender, e$receiver, e$weight, length(cell_id),
    source, receptor, diffusion, decay, uptake, production_rate,
    tolerance, max_iterations
  )
  names(result$concentration) <- cell_id
  result$parameters <- as.list(parameters)
  result$parameters$max_iterations <- max_iterations
  result$graph_type <- graph$graph_type
  class(result) <- "SpatialESSField"
  result
}
