.spatialess_component_genes <- function(name, table, columns, genes,
                                        require_all = FALSE) {
  name <- as.character(name)
  if (length(name) != 1L || is.na(name) || !nzchar(name)) return(character())
  values <- if (!is.null(table) && name %in% rownames(table)) {
    unlist(table[name, columns, drop = FALSE], use.names = FALSE)
  } else {
    name
  }
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (require_all && any(!values %in% genes)) return(character())
  values[values %in% genes]
}

.spatialess_index_matrix <- function(sets, genes) {
  positions <- seq_along(genes)
  names(positions) <- genes
  indices <- lapply(sets, function(x) unname(positions[x]))
  width <- max(1L, max(lengths(indices)))
  output <- matrix(0L, nrow = length(indices), ncol = width)
  for (i in seq_along(indices)) {
    if (length(indices[[i]])) output[i, seq_along(indices[[i]])] <- indices[[i]]
  }
  output
}

.spatialess_group_info <- function(group, graph = NULL) {
  if (is.null(group) && !is.null(graph)) group <- graph$cells$group_code
  if (is.null(group)) stop("group is required.", call. = FALSE)
  if (is.factor(group)) {
    code <- as.integer(group)
    level <- levels(group)
  } else if (is.numeric(group) &&
             all(group == as.integer(group), na.rm = TRUE)) {
    code <- as.integer(group)
    level <- if (!is.null(graph) && !is.null(graph$cells$group_levels)) {
      graph$cells$group_levels
    } else {
      as.character(seq_len(max(code, na.rm = TRUE)))
    }
  } else {
    value <- as.character(group)
    level <- unique(value)
    code <- match(value, level)
  }
  if (anyNA(code) || any(code < 1L)) stop("group contains invalid labels.", call. = FALSE)
  count <- max(code)
  if (!identical(sort(unique(code)), seq_len(count))) {
    stop("group codes must not contain unused levels.", call. = FALSE)
  }
  if (length(level) < count) level <- as.character(seq_len(count))
  list(code = code, levels = level, count = count)
}

#' Resolve CellChat ligand, receptor, complex and cofactor definitions
#'
#' Complexes retain all listed subunits. Cofactor sets follow CellChat by using
#' only genes present in the expression matrix.
#'
#' @param lr A ligand-receptor interaction data frame.
#' @param genes Unique genes available in the expression matrix.
#' @param complex,cofactor CellChat complex and cofactor definition tables.
#' @param ligand_col,receptor_col Ligand and receptor column names.
#' @export
prepare_cellchat_lr_components <- function(lr, genes, complex, cofactor,
                                           ligand_col = "ligand",
                                           receptor_col = "receptor") {
  required <- c(ligand_col, receptor_col)
  if (!all(required %in% colnames(lr))) {
    stop("lr is missing ligand or receptor columns.", call. = FALSE)
  }
  genes <- as.character(genes)
  if (!length(genes) || anyNA(genes) || anyDuplicated(genes)) {
    stop("genes must be unique non-missing identifiers.", call. = FALSE)
  }
  if (is.null(rownames(complex)) || is.null(rownames(cofactor))) {
    stop("complex and cofactor tables require row names.", call. = FALSE)
  }
  complex_columns <- grep("^subunit", colnames(complex))
  cofactor_columns <- grep("cofactor", colnames(cofactor))
  if (!length(complex_columns) || !length(cofactor_columns)) {
    stop("Cannot identify complex subunit or cofactor columns.", call. = FALSE)
  }

  resolve_complex <- function(x) .spatialess_component_genes(
    x, complex, complex_columns, genes, require_all = TRUE
  )
  resolve_cofactor <- function(x) .spatialess_component_genes(
    x, cofactor, cofactor_columns, genes, require_all = FALSE
  )
  ligand <- lapply(lr[[ligand_col]], resolve_complex)
  receptor <- lapply(lr[[receptor_col]], resolve_complex)
  valid <- lengths(ligand) > 0L & lengths(receptor) > 0L
  if (any(!valid)) {
    bad <- which(!valid)
    stop(sprintf("%d LR records have unresolved ligand/receptor components (first row: %d).",
                 length(bad), bad[[1L]]), call. = FALSE)
  }

  get_column <- function(name) {
    if (name %in% colnames(lr)) lr[[name]] else rep("", nrow(lr))
  }
  agonist_name <- as.character(get_column("agonist"))
  antagonist_name <- as.character(get_column("antagonist"))
  components <- list(
    ligand = ligand,
    receptor = receptor,
    co_a = lapply(get_column("co_A_receptor"), resolve_cofactor),
    co_i = lapply(get_column("co_I_receptor"), resolve_cofactor),
    agonist = lapply(agonist_name, resolve_cofactor),
    antagonist = lapply(antagonist_name, resolve_cofactor)
  )
  used_genes <- unique(unlist(components, recursive = TRUE, use.names = FALSE))
  structure(
    c(components, list(
      has_agonist = !is.na(agonist_name) & nzchar(agonist_name),
      has_antagonist = !is.na(antagonist_name) & nzchar(antagonist_name),
      genes = genes[genes %in% used_genes],
      lr = lr
    )),
    class = "SpatialESSLRComponents"
  )
}

#' Compute CellChat-compatible sparse group triMeans
#'
#' Only genes referenced by the supplied LR components are selected. Sparse
#' zeros remain implicit while exact type-7 quartiles are computed.
#'
#' @param expression A non-negative gene-by-cell matrix.
#' @param group One group label per cell.
#' @param components Optional resolved LR components.
#' @param genes Optional genes to summarize.
#' @param normalize Whether to divide by the global expression maximum.
#' @export
summarize_cellchat_trimean <- function(expression, group, components = NULL,
                                       genes = NULL, normalize = TRUE) {
  if (!inherits(expression, "Matrix")) {
    expression <- Matrix::Matrix(expression, sparse = TRUE)
  }
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression must have gene and cell names.", call. = FALSE)
  }
  if (any(!is.finite(expression@x)) || any(expression@x < 0)) {
    stop("expression must contain finite non-negative values.", call. = FALSE)
  }
  info <- .spatialess_group_info(group)
  if (length(info$code) != ncol(expression)) {
    stop("group must provide one label per expression column.", call. = FALSE)
  }
  if (!is.null(components)) {
    if (!inherits(components, "SpatialESSLRComponents")) {
      stop("components must come from prepare_cellchat_lr_components().", call. = FALSE)
    }
    genes <- components$genes
  }
  if (is.null(genes)) genes <- rownames(expression)
  genes <- as.character(genes)
  if (!length(genes) || any(!genes %in% rownames(expression))) {
    stop("Requested genes are absent from expression.", call. = FALSE)
  }
  expression_max <- if (length(expression@x)) max(expression@x) else 0
  if (normalize && expression_max <= 0) {
    stop("expression has no positive values.", call. = FALSE)
  }
  selected <- expression[genes, , drop = FALSE]
  if (normalize) selected <- selected / expression_max
  cell_by_gene <- methods::as(Matrix::t(selected), "dgCMatrix")
  average <- group_tri_mean_dgc_cpp(cell_by_gene, info$code, info$count)
  rownames(average) <- genes
  colnames(average) <- info$levels[seq_len(info$count)]
  structure(
    list(average = average, genes = genes, group_code = info$code,
         group_levels = info$levels[seq_len(info$count)],
         expression_max = expression_max, normalize = normalize),
    class = "SpatialESSGroupExpression"
  )
}

#' Collapse a cell CSR graph to sparse group-pair support
#'
#' @param graph A compact SpatialESS CSR graph.
#' @param group Optional group labels aligned to graph cells.
#' @param max_group_pairs Maximum sparse group-pair records.
#' @export
build_group_support_csr <- function(graph, group = NULL, max_group_pairs = 1e7) {
  if (!inherits(graph, "SpatialESSCSRGraph")) {
    stop("graph must be a SpatialESSCSRGraph.", call. = FALSE)
  }
  info <- .spatialess_group_info(group, graph)
  if (length(info$code) != length(graph$cells$cell_id)) {
    stop("group must provide one label per graph cell.", call. = FALSE)
  }
  pairs <- as.data.frame(aggregate_group_support_csr_cpp(
    graph$offsets, graph$neighbors, graph$weights, info$code, info$count,
    as.double(max_group_pairs)
  ), stringsAsFactors = FALSE)
  pairs$sender_group_name <- info$levels[pairs$sender_group]
  pairs$receiver_group_name <- info$levels[pairs$receiver_group]
  structure(
    list(group_pairs = pairs, group_code = info$code,
         group_levels = info$levels[seq_len(info$count)],
         cell_edges = length(graph$neighbors)),
    class = "SpatialESSGroupSupport"
  )
}

#' Score CellChat molecular probabilities on sparse spatial group support
#'
#' The molecular probability follows CellChat's ordering exactly: group
#' triMean, complex geometric mean, coreceptor adjustment, Hill activation,
#' then sender/receiver agonist and antagonist factors. Spatial weighting is
#' reported separately and does not alter the molecular column.
#'
#' @param group_expression A prepared group-expression summary.
#' @param group_support Spatially supported group pairs.
#' @param components Resolved CellChat LR components.
#' @param Kh,n Positive Hill-function parameters.
#' @param max_output_records Maximum LR-group output records.
#' @export
score_cellchat_group_support <- function(group_expression, group_support,
                                         components, Kh = 0.5, n = 1,
                                         max_output_records = 1e8) {
  if (!inherits(group_expression, "SpatialESSGroupExpression") ||
      !inherits(group_support, "SpatialESSGroupSupport") ||
      !inherits(components, "SpatialESSLRComponents")) {
    stop("Invalid group expression, group support or LR components object.", call. = FALSE)
  }
  average <- group_expression$average
  if (!identical(colnames(average), group_support$group_levels)) {
    stop("Group expression and graph support levels differ.", call. = FALSE)
  }
  if (any(!components$genes %in% rownames(average))) {
    stop("Group expression is missing referenced LR genes.", call. = FALSE)
  }
  indices <- lapply(components[c(
    "ligand", "receptor", "co_a", "co_i", "agonist", "antagonist"
  )], .spatialess_index_matrix, genes = rownames(average))
  support <- group_support$group_pairs
  result <- score_cellchat_group_support_cpp(
    average, support$sender_group, support$receiver_group,
    support$supported_edges, support$sum_spatial_weight,
    support$mean_spatial_weight,
    indices$ligand, indices$receptor, indices$co_a, indices$co_i,
    indices$agonist, indices$antagonist,
    components$has_agonist, components$has_antagonist,
    Kh, n, as.double(max_output_records)
  )
  pairs <- as.data.frame(result$group_pairs, stringsAsFactors = FALSE)
  lr <- components$lr
  interaction <- if ("interaction_name" %in% colnames(lr)) {
    as.character(lr$interaction_name)
  } else {
    paste(lr$ligand, lr$receptor, sep = "_")
  }
  if (nrow(pairs)) {
    pairs$interaction_name <- interaction[pairs$lr_index]
    pairs$ligand <- as.character(lr$ligand)[pairs$lr_index]
    pairs$receptor <- as.character(lr$receptor)[pairs$lr_index]
    pairs$sender_group_name <- group_support$group_levels[pairs$sender_group]
    pairs$receiver_group_name <- group_support$group_levels[pairs$receiver_group]
  }
  diagnostics <- data.frame(
    lr_index = seq_len(nrow(lr)), interaction_name = interaction,
    ligand = as.character(lr$ligand), receptor = as.character(lr$receptor),
    emitted_group_pairs = as.numeric(result$emitted_group_pairs),
    spatially_supported_group_pairs = nrow(support)
  )
  structure(
    list(group_pairs = pairs, diagnostics = diagnostics, lr = lr,
         parameters = list(Kh = Kh, n = n,
                           spatial_columns = c("mean_weighted_probability",
                                               "sum_weighted_probability"))),
    class = "SpatialESSCellChatCompatibility"
  )
}
