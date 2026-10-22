#' Score simple ligand-receptor interactions on supported spatial edges
#'
#' This is the first edge-local molecular primitive for SpatialESS. It accepts
#' single-gene ligand-receptor pairs only. Complexes and cofactors are added in
#' a later, separately validated kernel.
#'
#' @param expression A non-negative gene-by-cell matrix.
#' @param graph A SpatialESS graph.
#' @param lr A ligand-receptor interaction data frame.
#' @param ligand_col,receptor_col,interaction_col LR column names.
#' @param Kh,n Positive Hill-function parameters.
#' @param normalize Whether to divide by the global expression maximum.
#' @param min_expression Non-negative expression-support threshold.
#' @param max_output_edges Maximum LR-edge output records.
#' @export
score_simple_lr_edges <- function(expression, graph, lr,
                                  ligand_col = "ligand",
                                  receptor_col = "receptor",
                                  interaction_col = "interaction_name",
                                  Kh = 0.5, n = 1,
                                  normalize = TRUE,
                                  min_expression = 0,
                                  max_output_edges = 1e7) {
  if (!inherits(graph, "SpatialESSGraph")) {
    stop("graph must be a SpatialESSGraph.", call. = FALSE)
  }
  if (!inherits(expression, "Matrix")) expression <- Matrix::Matrix(expression, sparse = TRUE)
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression must have gene row names and cell column names.", call. = FALSE)
  }
  cell_id <- graph$cells$cell_id
  if (!identical(colnames(expression), cell_id)) {
    if (length(colnames(expression)) != length(cell_id) ||
        !setequal(colnames(expression), cell_id)) {
      stop("Expression columns and graph cell identifiers do not match.", call. = FALSE)
    }
    expression <- expression[, cell_id, drop = FALSE]
  }
  required <- c(ligand_col, receptor_col)
  if (!all(required %in% colnames(lr))) {
    stop("lr is missing ligand or receptor columns.", call. = FALSE)
  }
  if (!interaction_col %in% colnames(lr)) {
    lr[[interaction_col]] <- paste(lr[[ligand_col]], lr[[receptor_col]], sep = "_")
  }
  if (any(!is.finite(expression@x)) || any(expression@x < 0)) {
    stop("expression must contain finite non-negative values.", call. = FALSE)
  }
  if (!is.finite(Kh) || Kh <= 0 || !is.finite(n) || n <= 0) {
    stop("Kh and n must be positive finite values.", call. = FALSE)
  }
  if (!is.finite(min_expression) || min_expression < 0) {
    stop("min_expression must be finite and non-negative.", call. = FALSE)
  }
  max_output_edges <- as.double(max_output_edges)
  if (!is.finite(max_output_edges) || max_output_edges < 0) {
    stop("max_output_edges must be finite and non-negative.", call. = FALSE)
  }

  genes <- rownames(expression)
  simple <- !is.na(lr[[ligand_col]]) & !is.na(lr[[receptor_col]]) &
    lr[[ligand_col]] %in% genes & lr[[receptor_col]] %in% genes
  lr_use <- lr[simple, , drop = FALSE]
  if (!nrow(lr_use)) {
    return(structure(list(edges = data.frame(), lr = lr_use,
                          diagnostics = data.frame()),
                     class = "SpatialESSEdgeScores"))
  }

  global_max <- if (length(expression@x)) max(expression@x) else 0
  if (normalize && global_max <= 0) stop("expression has no positive values.", call. = FALSE)
  divisor <- if (normalize) global_max else 1
  e <- graph$edges
  results <- vector("list", nrow(lr_use))
  diagnostics <- vector("list", nrow(lr_use))
  emitted <- 0

  for (q in seq_len(nrow(lr_use))) {
    ligand <- as.character(lr_use[[ligand_col]][q])
    receptor <- as.character(lr_use[[receptor_col]][q])
    ligand_expression <- as.numeric(expression[ligand, , drop = TRUE]) / divisor
    receptor_expression <- as.numeric(expression[receptor, , drop = TRUE]) / divisor
    remaining <- max_output_edges - emitted
    scored <- score_simple_lr_edge_cpp(
      e$sender, e$receiver, e$distance, e$weight,
      ligand_expression, receptor_expression,
      Kh, n, min_expression, remaining
    )
    out <- as.data.frame(scored$edges, stringsAsFactors = FALSE)
    if (nrow(out)) {
      out$lr_index <- q
      out$interaction_name <- as.character(lr_use[[interaction_col]][q])
      out$ligand <- ligand
      out$receptor <- receptor
      out$sender_id <- cell_id[out$sender]
      out$receiver_id <- cell_id[out$receiver]
      out <- out[, c("lr_index", "interaction_name", "ligand", "receptor",
                     "sender", "receiver", "sender_id", "receiver_id",
                     "distance", "spatial_weight", "ligand_expression",
                     "receptor_expression", "molecular_product", "hill",
                     "spatial_score"), drop = FALSE]
    }
    results[[q]] <- out
    emitted <- emitted + nrow(out)
    diagnostics[[q]] <- data.frame(
      lr_index = q,
      interaction_name = as.character(lr_use[[interaction_col]][q]),
      ligand = ligand,
      receptor = receptor,
      ligand_supported_cells = sum(ligand_expression > min_expression),
      receptor_supported_cells = sum(receptor_expression > min_expression),
      graph_edges = nrow(e),
      supported_edges = nrow(out),
      supported_fraction = if (nrow(e)) nrow(out) / nrow(e) else 0
    )
  }
  nonempty <- vapply(results, nrow, integer(1)) > 0L
  edge_result <- if (any(nonempty)) do.call(rbind, results[nonempty]) else data.frame()
  rownames(edge_result) <- NULL
  structure(
    list(edges = edge_result, lr = lr_use,
         diagnostics = do.call(rbind, diagnostics),
         parameters = list(Kh = Kh, n = n, normalize = normalize,
                           expression_max = global_max,
                           min_expression = min_expression)),
    class = "SpatialESSEdgeScores"
  )
}
