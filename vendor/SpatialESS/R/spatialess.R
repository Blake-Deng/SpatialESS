#' SpatialESS streamed spatial communication
#'
#' Reproduces the individual-cell probability and group-average ordering used
#' by SpatialCellChat V3 while streaming over a compact CSR radius graph.
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
#' @param lr_missing How to handle incomplete mandatory LR genes: `"error"`
#'   (default) preserves strict checking; `"filter"` explicitly excludes
#'   incomplete LR with a warning and stores the full audit in `$lr_filter`.
#' @param aggregation Group aggregation arithmetic. `"reference"` (default)
#'   uses the same full Matrix multiplication and crossproduct as the official
#'   SpatialCellChat reference, preserving its floating-point reduction order.
#'   `"stream"` retains the compact native sum but can differ at exact permutation
#'   boundaries and must not be described as bitwise reference equivalent.
#' @param max_reference_entries Maximum entries in each dense reference
#'   workspace matrix (`cells * groups` or `groups * groups`). Several matrices
#'   are allocated, so this is not a total process-memory limit. The default is
#'   `1e8`. Exceeding it stops rather than silently changing aggregation.
#' @return A SpatialESSResult with records, LR table, CSR graph, diagnostics
#'   and parameters. When `lr_missing = "filter"`, `lr_filter` contains the
#'   original-order coverage audit, excluded rows and missing-gene table.
#' @export
spatialess <- function(
    expression, coordinates, group, lr, complex, cofactor,
    ratio = 1, tol = 5, interaction.range = 250,
    scale.distance = 0.01, contact.range = 10,
    Kh = 0.5, n = 1, use.AGAN = TRUE,
    contact.dependent = TRUE, contact.dependent.forced = FALSE,
    min.percent = 0.1, min.cells.sr = 5L,
    nboot = 0L, seed.use = 1L,
    max_edges = 5e8, max_output_records = 1e8,
    max_active_edges_per_lr = 5e8,
    lr_missing = c("error", "filter"),
    aggregation = c("reference", "stream"),
    max_reference_entries = 1e8) {
  lr_missing <- match.arg(lr_missing)
  aggregation <- match.arg(aggregation)
  if (!is.numeric(ratio) || length(ratio) != 1L ||
      !is.finite(ratio) || ratio <= 0) {
    stop("ratio must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(max_reference_entries) || length(max_reference_entries) != 1L ||
      !is.finite(max_reference_entries) || max_reference_entries < 1) {
    stop("max_reference_entries must be a positive finite entry limit.", call. = FALSE)
  }
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
  reference_cell_group_entries <- as.double(ncol(expression)) * info$count
  reference_group_entries <- as.double(info$count)^2
  if (aggregation == "reference" &&
      max(reference_cell_group_entries, reference_group_entries) > max_reference_entries) {
    stop(paste0("Reference aggregation exceeds max_reference_entries. Increase ",
      "the explicit workspace limit, or explicitly choose aggregation='stream' ",
      "without claiming bitwise reference equivalence."), call. = FALSE)
  }
  nboot <- as.integer(nboot)
  seed.use <- as.integer(seed.use)
  if (is.na(nboot) || nboot < 0L || is.na(seed.use)) {
    stop("nboot and seed.use are invalid.", call. = FALSE)
  }

  lr <- as.data.frame(lr, stringsAsFactors = FALSE)
  lr_filter <- NULL
  if (lr_missing == "filter") {
    lr_filter <- filter_cellchat_lr(lr, rownames(expression), complex, cofactor)
    lr <- lr_filter$lr
  }
  if (nrow(lr) == 0L) {
    stop("No complete LR records remain. Inspect filter_cellchat_lr()$audit and $missing_genes; check the panel and gene identifiers.", call. = FALSE)
  }
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

  # The reference computes distances in the original coordinate units and only
  # then multiplies retained distances by ratio. Pre-scaling the coordinates can
  # change last-bit rounding and membership at the distance threshold.
  radius_um <- interaction.range + tol
  radius_coordinate <- radius_um / ratio
  if (length(radius_coordinate) != 1L || !is.finite(radius_coordinate) ||
      radius_coordinate <= 0) {
    stop("The converted interaction radius must be positive and finite.", call. = FALSE)
  }
  rownames(coordinates) <- colnames(expression)
  graph <- build_radius_graph_csr(
    coordinates, radius = radius_coordinate,
    weight = "binary", store_distance = TRUE, max_edges = max_edges
  )
  # The general graph builder admits a numerical candidate tolerance. Match the
  # reference's strict distance mask here, including exclusion of distinct
  # cells at identical coordinates. Self interactions are added separately.
  keep_distance <- graph$distances > 0 & graph$distances <= radius_coordinate
  if (!all(keep_distance)) {
    removed_rows <- findInterval(which(!keep_distance) - 1, graph$offsets)
    graph$degree <- graph$degree - tabulate(removed_rows, nrow(coordinates))
    graph$weighted_degree <- graph$degree
    graph$offsets <- c(0, cumsum(graph$degree))
    graph$neighbors <- graph$neighbors[keep_distance]
    graph$weights <- graph$weights[keep_distance]
    graph$distances <- graph$distances[keep_distance]
  }
  graph$distances <- graph$distances * ratio
  graph$parameters$radius <- radius_um
  graph$parameters$scale <- radius_um / 2
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
    # Match the reference V3 future RNG advance regardless of optional-package
    # discovery. future.apply is a declared runtime dependency.
    invisible(future.apply::future_lapply(1L, identity, future.seed = TRUE))
    permutations <- replicate(
      nboot, sample.int(ncol(expression), size = ncol(expression))
    )
    if (nboot == 1L) {
      permutations <- matrix(permutations, nrow = ncol(expression), ncol = 1L)
    }
    storage.mode(permutations) <- "integer"
  }

  result <- spatialess_stream_cpp(
    cell_by_gene, graph$offsets, graph$neighbors, graph$distances,
    info$code, info$count,
    indices$ligand, indices$receptor, indices$co_a, indices$co_i,
    indices$agonist, indices$antagonist,
    components$has_agonist, components$has_antagonist,
    contact_lr, permutations, Kh, n, scale.distance,
    contact.range + tol, min.percent, as.integer(min.cells.sr),
    isTRUE(use.AGAN), as.double(max_output_records),
    as.double(max_active_edges_per_lr), aggregation,
    as.double(max_reference_entries)
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
      lr_filter = lr_filter,
      graph = graph,
      diagnostics = list(
        active_edges_per_lr = result$active_edges_per_lr,
        minimum_nonzero_distance = result$minimum_nonzero_distance,
        self_spatial_weight = result$self_spatial_weight,
        aggregation = aggregation,
        reference_cell_group_entries = if (aggregation == "reference") {
          reference_cell_group_entries
        } else 0,
        reference_group_matrix_entries = if (aggregation == "reference") {
          reference_group_entries
        } else 0,
        reference_workspace_note = if (aggregation == "reference") {
          "Full Matrix/BLAS reference order; several dense cell-by-group/group-by-group workspaces, not a cell-by-cell dense matrix."
        } else {
          "Native flat sums; exact permutation-boundary agreement with the reference is not guaranteed."
        }
      ),
      parameters = list(
        engine = "SpatialESS_CSR_stream",
        percentage_filter = "SpatialCellChat_V3_format_digits1",
        ratio = ratio, tol = tol, interaction.range = interaction.range,
        scale.distance = scale.distance, contact.range = contact.range,
        Kh = Kh, n = n, use.AGAN = isTRUE(use.AGAN),
        min.percent = min.percent, min.cells.sr = min.cells.sr,
        nboot = nboot, seed.use = seed.use,
        aggregation = aggregation, max_reference_entries = max_reference_entries
      )
    ),
    class = "SpatialESSResult"
  )
}

# V3 formats each aggregated expression-fraction column before comparing it.
# Keep the data.frame dispatch and character comparison, including R options.
.spatialess_v3_percent_mask <- function(positive, sizes, min.percent) {
  fractions <- data.frame(
    ligand = positive[, 1L] / sizes,
    receptor = positive[, 2L] / sizes
  )
  unname(format(fractions, digits = 1) >= min.percent)
}
