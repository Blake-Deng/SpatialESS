#' Stream simple LR scores over a compact CSR graph
#'
#' Supported cell edges are accumulated directly into sparse sender-group by
#' receiver-group records. No LR-by-edge table or dense group square is formed.
#'
#' @param expression A non-negative gene-by-cell matrix.
#' @param graph A compact SpatialESS CSR graph.
#' @param lr A ligand-receptor interaction data frame.
#' @param group Optional cell-aligned group labels.
#' @param ligand_col,receptor_col,interaction_col LR column names.
#' @param Kh,n Positive Hill-function parameters.
#' @param normalize Whether to divide by the global expression maximum.
#' @param min_expression Non-negative expression-support threshold.
#' @param max_group_pairs Maximum combined sparse output records.
#' @export
aggregate_simple_lr_csr <- function(expression, graph, lr, group = NULL,
                                    ligand_col = "ligand",
                                    receptor_col = "receptor",
                                    interaction_col = "interaction_name",
                                    Kh = 0.5, n = 1,
                                    normalize = TRUE,
                                    min_expression = 0,
                                    max_group_pairs = 1e7) {
  if (!inherits(graph, "SpatialESSCSRGraph")) {
    stop("graph must be a SpatialESSCSRGraph.", call. = FALSE)
  }
  if (!inherits(expression, "Matrix")) expression <- Matrix::Matrix(expression, sparse = TRUE)
  cell_id <- graph$cells$cell_id
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression must have gene and cell names.", call. = FALSE)
  }
  if (!identical(colnames(expression), cell_id)) {
    if (length(colnames(expression)) != length(cell_id) ||
        !setequal(colnames(expression), cell_id)) {
      stop("Expression columns and graph cell identifiers do not match.", call. = FALSE)
    }
    expression <- expression[, cell_id, drop = FALSE]
  }
  if (is.null(group)) group <- graph$cells$group_code
  if (is.null(group) || length(group) != length(cell_id)) {
    stop("group must provide one label per graph cell.", call. = FALSE)
  }
  if (is.factor(group)) {
    group_levels <- levels(group)
    group_code <- as.integer(group)
  } else if (is.numeric(group) && !is.null(graph$cells$group_levels) &&
             all(group == as.integer(group), na.rm = TRUE)) {
    group_code <- as.integer(group)
    group_levels <- graph$cells$group_levels
  } else {
    group_factor <- factor(group, levels = unique(as.character(group)))
    group_code <- as.integer(group_factor)
    group_levels <- levels(group_factor)
  }
  if (anyNA(group_code) || any(group_code < 1L)) {
    stop("group contains missing or invalid labels.", call. = FALSE)
  }
  K <- max(group_code)
  if (length(group_levels) < K) group_levels <- as.character(seq_len(K))

  required <- c(ligand_col, receptor_col)
  if (!all(required %in% colnames(lr))) stop("lr is missing ligand or receptor columns.", call. = FALSE)
  if (!interaction_col %in% colnames(lr)) {
    lr[[interaction_col]] <- paste(lr[[ligand_col]], lr[[receptor_col]], sep = "_")
  }
  if (any(!is.finite(expression@x)) || any(expression@x < 0)) {
    stop("expression must contain finite non-negative values.", call. = FALSE)
  }
  if (!is.finite(Kh) || Kh <= 0 || !is.finite(n) || n <= 0 ||
      !is.finite(min_expression) || min_expression < 0) {
    stop("Invalid scoring parameter.", call. = FALSE)
  }

  genes <- rownames(expression)
  simple <- !is.na(lr[[ligand_col]]) & !is.na(lr[[receptor_col]]) &
    lr[[ligand_col]] %in% genes & lr[[receptor_col]] %in% genes
  lr_use <- lr[simple, , drop = FALSE]
  if (!nrow(lr_use)) {
    return(structure(list(group_pairs = data.frame(), diagnostics = data.frame(),
                          lr = lr_use), class = "SpatialESSAggregates"))
  }
  expression_max <- if (length(expression@x)) max(expression@x) else 0
  if (normalize && expression_max <= 0) stop("expression has no positive values.", call. = FALSE)
  divisor <- if (normalize) expression_max else 1

  outputs <- vector("list", nrow(lr_use))
  diagnostics <- vector("list", nrow(lr_use))
  emitted <- 0
  for (q in seq_len(nrow(lr_use))) {
    ligand <- as.character(lr_use[[ligand_col]][q])
    receptor <- as.character(lr_use[[receptor_col]][q])
    ligand_expression <- as.numeric(expression[ligand, , drop = TRUE]) / divisor
    receptor_expression <- as.numeric(expression[receptor, , drop = TRUE]) / divisor
    result <- aggregate_simple_lr_csr_cpp(
      graph$offsets, graph$neighbors, graph$weights,
      ligand_expression, receptor_expression, group_code, K,
      Kh, n, min_expression, max_group_pairs
    )
    out <- as.data.frame(result$group_pairs, stringsAsFactors = FALSE)
    emitted <- emitted + nrow(out)
    if (emitted > max_group_pairs) {
      stop("Combined sparse group-pair output exceeded max_group_pairs.", call. = FALSE)
    }
    if (nrow(out)) {
      out$lr_index <- q
      out$interaction_name <- as.character(lr_use[[interaction_col]][q])
      out$ligand <- ligand
      out$receptor <- receptor
      out$sender_group_name <- group_levels[out$sender_group]
      out$receiver_group_name <- group_levels[out$receiver_group]
      out <- out[, c("lr_index", "interaction_name", "ligand", "receptor",
                     "sender_group", "receiver_group",
                     "sender_group_name", "receiver_group_name",
                     "supported_edges", "sum_hill", "sum_spatial_score",
                     "mean_spatial_score"), drop = FALSE]
    }
    outputs[[q]] <- out
    diagnostics[[q]] <- data.frame(
      lr_index = q,
      interaction_name = as.character(lr_use[[interaction_col]][q]),
      ligand = ligand, receptor = receptor,
      ligand_supported_cells = result$ligand_supported_cells,
      receptor_supported_cells = result$receptor_supported_cells,
      visited_out_edges = result$visited_out_edges,
      supported_edges = result$supported_edges,
      sparse_group_pairs = nrow(out),
      total_graph_edges = length(graph$neighbors)
    )
  }
  nonempty <- vapply(outputs, nrow, integer(1)) > 0L
  pairs <- if (any(nonempty)) do.call(rbind, outputs[nonempty]) else data.frame()
  rownames(pairs) <- NULL
  structure(
    list(group_pairs = pairs, diagnostics = do.call(rbind, diagnostics),
         lr = lr_use,
         parameters = list(Kh = Kh, n = n, normalize = normalize,
                           expression_max = expression_max,
                           min_expression = min_expression)),
    class = "SpatialESSAggregates"
  )
}
