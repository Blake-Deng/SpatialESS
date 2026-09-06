#' SpatialCellChat v3-compatible streamed spatial communication
#'
#' Reproduces the individual-cell probability and group-average ordering used
#' by SpatialCellChat v3 while streaming over a compact CSR radius graph.
#'
#' @param expression Non-negative gene-by-cell matrix.
#' @param coordinates Cell-by-2/3 coordinate matrix.
#' @param group Cell-aligned factor or labels.
#' @param lr CellChat ligand-receptor interaction table.
#' @param complex,cofactor CellChat database component tables.
#' @param ratio Coordinate-unit to micrometre conversion.
#' @param tol Distance tolerance in micrometres.
#' @param interaction.range,contact.range SpatialCellChat distance parameters.
#' @param scale.distance SpatialCellChat inverse-distance scaling parameter.
#' @param Kh,n Hill-function parameters.
#' @param use.AGAN Whether to include agonist and antagonist factors.
#' @param contact.dependent,contact.dependent.forced Contact-mode controls.
#' @param min.percent,min.cells.sr Group-average filters.
#' @param nboot,seed.use Global label-permutation settings.
#' @param max_edges,max_output_records,max_active_edges_per_lr Memory guards.
#' @export
spatial_v3_compatible <- function(
    expression, coordinates, group, lr, complex, cofactor,
    ratio = 1, tol = 5, interaction.range = 250,
    scale.distance = 0.01, contact.range = 10,
    Kh = 0.5, n = 1, use.AGAN = TRUE,
    contact.dependent = TRUE, contact.dependent.forced = FALSE,
    min.percent = 0.1, min.cells.sr = 5L,
    nboot = 0L, seed.use = 1L,
    max_edges = 5e8, max_output_records = 1e8,
    max_active_edges_per_lr = 5e8) {
  if (!inherits(expression, "Matrix")) {
    expression <- Matrix::Matrix(expression, sparse = TRUE)
  }
  expression <- methods::as(expression, "dgCMatrix")
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression requires gene and cell names.", call. = FALSE)
  }
  coordinates <- as.matrix(coordinates)
  if (nrow(coordinates) != ncol(expression) ||
      !ncol(coordinates) %in% c(2L, 3L)) {
    stop("coordinates must provide two or three columns per cell.", call. = FALSE)
  }
  if (!is.null(rownames(coordinates)) &&
      !identical(rownames(coordinates), colnames(expression))) {
    stop("coordinate and expression cell order differs.", call. = FALSE)
  }
  info <- .spatialess_group_info(group)
  if (length(info$code) != ncol(expression)) {
    stop("group must provide one label per cell.", call. = FALSE)
  }
  nboot <- as.integer(nboot)
  seed.use <- as.integer(seed.use)
  if (is.na(nboot) || nboot < 0L || is.na(seed.use)) {
    stop("nboot and seed.use are invalid.", call. = FALSE)
  }

  lr <- as.data.frame(lr, stringsAsFactors = FALSE)
  if ("annotation" %in% colnames(lr) &&
      length(unique(lr$annotation)) > 1L) {
    annotation <- factor(
      lr$annotation,
      levels = c("Secreted Signaling", "ECM-Receptor",
                 "Non-protein Signaling", "Cell-Cell Contact")
    )
    lr <- lr[order(annotation), , drop = FALSE]
  }
  components <- prepare_cellchat_lr_components(
    lr, rownames(expression), complex, cofactor
  )
  expression_max <- if (length(expression@x)) max(expression@x) else 0
  if (!is.finite(expression_max) || expression_max <= 0) {
    stop("expression has no positive finite value.", call. = FALSE)
  }
  selected <- expression[components$genes, , drop = FALSE]
  selected@x <- selected@x / expression_max
  cell_by_gene <- methods::as(Matrix::t(selected), "dgCMatrix")
  indices <- lapply(components[c(
    "ligand", "receptor", "co_a", "co_i", "agonist", "antagonist"
  )], .spatialess_index_matrix, genes = components$genes)

  coordinates_um <- coordinates * ratio
  rownames(coordinates_um) <- colnames(expression)
  graph <- build_radius_graph_csr(
    coordinates_um, radius = interaction.range + tol,
    weight = "binary", store_distance = TRUE, max_edges = max_edges
  )
  annotation <- if ("annotation" %in% colnames(lr)) {
    as.character(lr$annotation)
  } else {
    rep("Secreted Signaling", nrow(lr))
  }
  contact_lr <- if (isTRUE(contact.dependent.forced)) {
    rep(TRUE, nrow(lr))
  } else if (isTRUE(contact.dependent)) {
    annotation == "Cell-Cell Contact"
  } else {
    rep(FALSE, nrow(lr))
  }
  permutations <- matrix(integer(), nrow = ncol(expression), ncol = 0L)
  if (nboot > 0L) {
    set.seed(seed.use)
    if (requireNamespace("future.apply", quietly = TRUE)) {
      invisible(future.apply::future_lapply(1L, identity, future.seed = TRUE))
    }
    permutations <- replicate(
      nboot, sample.int(ncol(expression), size = ncol(expression))
    )
    if (nboot == 1L) {
      permutations <- matrix(permutations, nrow = ncol(expression), ncol = 1L)
    }
    storage.mode(permutations) <- "integer"
  }

  result <- spatial_v3_compatible_stream_cpp(
    cell_by_gene, graph$offsets, graph$neighbors, graph$distances,
    info$code, info$count,
    indices$ligand, indices$receptor, indices$co_a, indices$co_i,
    indices$agonist, indices$antagonist,
    components$has_agonist, components$has_antagonist,
    contact_lr, permutations, Kh, n, scale.distance,
    contact.range + tol, min.percent, as.integer(min.cells.sr),
    isTRUE(use.AGAN), as.double(max_output_records),
    as.double(max_active_edges_per_lr)
  )
  records <- as.data.frame(result$records, stringsAsFactors = FALSE)
  if (nrow(records)) {
    interaction <- if ("interaction_name" %in% colnames(lr)) {
      as.character(lr$interaction_name)
    } else {
      paste(lr$ligand, lr$receptor, sep = "_")
    }
    records$interaction_name <- interaction[records$lr_index]
    records$ligand <- as.character(lr$ligand)[records$lr_index]
    records$receptor <- as.character(lr$receptor)[records$lr_index]
    records$sender_group_name <- info$levels[records$sender_group]
    records$receiver_group_name <- info$levels[records$receiver_group]
  }
  structure(
    list(
      records = records,
      lr = lr,
      graph = graph,
      diagnostics = list(
        active_edges_per_lr = result$active_edges_per_lr,
        minimum_nonzero_distance = result$minimum_nonzero_distance,
        self_spatial_weight = result$self_spatial_weight
      ),
      parameters = list(
        engine = "SpatialCellChat_v3_compatible_CSR_stream",
        ratio = ratio, tol = tol, interaction.range = interaction.range,
        scale.distance = scale.distance, contact.range = contact.range,
        Kh = Kh, n = n, use.AGAN = isTRUE(use.AGAN),
        min.percent = min.percent, min.cells.sr = min.cells.sr,
        nboot = nboot, seed.use = seed.use
      )
    ),
    class = "SpatialESSV3Compatible"
  )
}
