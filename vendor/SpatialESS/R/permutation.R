#' Construct sample- and compartment-aware rectangular spatial blocks
#'
#' Cells are exchangeable only when sample, compartment and all spatial-bin
#' coordinates agree. The returned block identifiers are deterministic for a
#' fixed coordinate origin and block size.
#'
#' @param coords A finite cell-by-coordinate matrix.
#' @param block_size Positive scalar or per-dimension block sizes.
#' @param sample_id,compartment Optional cell-aligned exchangeability labels.
#' @param origin Optional scalar or coordinate origin vector.
#' @export
build_spatial_blocks <- function(coords, block_size, sample_id = NULL,
                                 compartment = NULL, origin = NULL) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  if (!is.numeric(coords) || nrow(coords) < 1L ||
      !ncol(coords) %in% c(2L, 3L) || any(!is.finite(coords))) {
    stop("coords must be a finite 2D or 3D numeric matrix.", call. = FALSE)
  }
  cell_id <- rownames(coords)
  if (is.null(cell_id)) cell_id <- as.character(seq_len(nrow(coords)))
  if (anyNA(cell_id) || anyDuplicated(cell_id)) {
    stop("Coordinate cell identifiers must be unique and non-missing.", call. = FALSE)
  }
  dimension <- ncol(coords)
  if (!length(block_size) %in% c(1L, dimension)) {
    stop("block_size must be scalar or match coordinate dimensions.", call. = FALSE)
  }
  block_size <- rep(as.double(block_size), length.out = dimension)
  if (length(block_size) != dimension || any(!is.finite(block_size)) ||
      any(block_size <= 0)) {
    stop("block_size must contain positive finite values.", call. = FALSE)
  }
  if (is.null(origin)) {
    origin <- floor(apply(coords, 2L, min) / block_size) * block_size
  }
  if (!length(origin) %in% c(1L, dimension)) {
    stop("origin must be scalar or match coordinate dimensions.", call. = FALSE)
  }
  origin <- rep(as.double(origin), length.out = dimension)
  if (length(origin) != dimension || any(!is.finite(origin))) {
    stop("origin must contain one finite value per coordinate dimension.",
         call. = FALSE)
  }
  if (is.null(sample_id)) sample_id <- rep("sample_1", nrow(coords))
  if (is.null(compartment)) compartment <- rep("all", nrow(coords))
  sample_id <- as.character(sample_id)
  compartment <- as.character(compartment)
  if (length(sample_id) != nrow(coords) || length(compartment) != nrow(coords) ||
      anyNA(sample_id) || anyNA(compartment)) {
    stop("sample_id and compartment must align with coords.", call. = FALSE)
  }

  bins <- sweep(coords, 2L, origin, FUN = "-")
  bins <- floor(sweep(bins, 2L, block_size, FUN = "/"))
  key_parts <- c(list(sample_id, compartment),
                 lapply(seq_len(dimension), function(j) bins[, j]))
  block_id <- do.call(paste, c(key_parts, sep = "|"))
  block_levels <- sort(unique(block_id), method = "radix")
  block_code <- match(block_id, block_levels)
  block_cells <- tabulate(block_code, nbins = length(block_levels))
  structure(
    list(
      block_code = as.integer(block_code),
      block_id = block_id,
      block_levels = block_levels,
      block_cells = block_cells,
      cell_id = cell_id,
      parameters = list(block_size = block_size, origin = origin,
                        dimensions = dimension,
                        preserves = c("sample_id", "compartment", "spatial_block"))
    ),
    class = "SpatialESSBlocks"
  )
}

#' Permute complete cell profiles within spatial exchangeability blocks
#'
#' The result maps each fixed target position to one source cell. Singleton
#' blocks remain fixed and no source cell crosses a block boundary.
#'
#' @param blocks A result from [build_spatial_blocks()].
#' @param seed One non-negative integer random seed.
#' @export
permute_within_spatial_blocks <- function(blocks, seed = 1L) {
  if (!inherits(blocks, "SpatialESSBlocks")) {
    stop("blocks must come from build_spatial_blocks().", call. = FALSE)
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed) || seed < 0L) stop("seed must be one integer.", call. = FALSE)
  set.seed(seed)
  permutation <- permute_within_strata_cpp(
    blocks$block_code, length(blocks$block_levels)
  )
  names(permutation) <- blocks$cell_id
  permutation
}

#' Prepare a referenced-gene sparse triMean kernel
#'
#' The expensive gene projection, normalization and CSC transpose are performed
#' once and reused across observed and permuted group summaries.
#'
#' @param expression A non-negative gene-by-cell matrix.
#' @param genes Optional unique genes to retain.
#' @param normalize Whether to divide by the global expression maximum.
#' @export
prepare_sparse_trimean <- function(expression, genes = NULL, normalize = TRUE) {
  if (!inherits(expression, "Matrix")) {
    expression <- Matrix::Matrix(expression, sparse = TRUE)
  }
  if (is.null(rownames(expression)) || is.null(colnames(expression))) {
    stop("expression must have gene and cell names.", call. = FALSE)
  }
  if (any(!is.finite(expression@x)) || any(expression@x < 0)) {
    stop("expression must contain finite non-negative values.", call. = FALSE)
  }
  if (is.null(genes)) genes <- rownames(expression)
  genes <- as.character(genes)
  if (!length(genes) || anyDuplicated(genes) ||
      any(!genes %in% rownames(expression))) {
    stop("genes must be unique identifiers present in expression.", call. = FALSE)
  }
  expression_max <- if (length(expression@x)) max(expression@x) else 0
  if (normalize && expression_max <= 0) {
    stop("expression has no positive values.", call. = FALSE)
  }
  selected <- expression[genes, , drop = FALSE]
  if (normalize) selected <- selected / expression_max
  cell_by_gene <- methods::as(Matrix::t(selected), "dgCMatrix")
  structure(
    list(cell_by_gene = cell_by_gene, genes = genes,
         cell_id = colnames(expression), expression_max = expression_max,
         normalize = normalize),
    class = "SpatialESSPreparedTriMean"
  )
}

#' Summarize a prepared sparse matrix by fixed target groups
#'
#' `source_for_target` optionally assigns permuted source profiles to target
#' positions without copying or reordering the expression matrix.
#'
#' @param prepared A result from [prepare_sparse_trimean()].
#' @param group One fixed target-group label per cell.
#' @param source_for_target Optional source-profile permutation.
#' @param group_levels Optional output group names.
#' @export
summarize_prepared_trimean <- function(prepared, group,
                                       source_for_target = NULL,
                                       group_levels = NULL) {
  if (!inherits(prepared, "SpatialESSPreparedTriMean")) {
    stop("prepared must come from prepare_sparse_trimean().", call. = FALSE)
  }
  info <- .spatialess_group_info(group)
  cell_count <- nrow(prepared$cell_by_gene)
  if (length(info$code) != cell_count) {
    stop("group must provide one label per prepared cell.", call. = FALSE)
  }
  if (is.null(group_levels)) group_levels <- info$levels[seq_len(info$count)]
  if (length(group_levels) != info$count || anyNA(group_levels)) {
    stop("group_levels must name every target group.", call. = FALSE)
  }
  if (is.null(source_for_target)) {
    average <- group_tri_mean_dgc_cpp(
      prepared$cell_by_gene, info$code, info$count
    )
  } else {
    source_for_target <- as.integer(source_for_target)
    if (length(source_for_target) != cell_count) {
      stop("source_for_target must provide one source per target cell.",
           call. = FALSE)
    }
    average <- group_tri_mean_dgc_permuted_cpp(
      prepared$cell_by_gene, info$code, source_for_target, info$count
    )
  }
  rownames(average) <- prepared$genes
  colnames(average) <- as.character(group_levels)
  structure(
    list(average = average, genes = prepared$genes,
         group_code = info$code, group_levels = as.character(group_levels),
         expression_max = prepared$expression_max,
         normalize = prepared$normalize),
    class = "SpatialESSGroupExpression"
  )
}

.spatialess_group_score_matrix <- function(average, group_support, components,
                                            Kh, n) {
  indices <- lapply(components[c(
    "ligand", "receptor", "co_a", "co_i", "agonist", "antagonist"
  )], .spatialess_index_matrix, genes = rownames(average))
  support <- group_support$group_pairs
  score_cellchat_group_matrix_cpp(
    average, support$sender_group, support$receiver_group,
    indices$ligand, indices$receptor, indices$co_a, indices$co_i,
    indices$agonist, indices$antagonist,
    components$has_agonist, components$has_antagonist, Kh, n
  )
}

#' Spatial block-permutation reference for CellChat-compatible scores
#'
#' Complete cell profiles are permuted within sample-compartment-spatial blocks.
#' The cell graph, target groups and sparse group-support domain remain fixed.
#' P-values use CellChat's strict-tail empirical convention by default.
#'
#' @param prepared A prepared sparse triMean kernel.
#' @param group Fixed target-group labels.
#' @param group_support Spatially supported group pairs.
#' @param components Resolved CellChat LR components.
#' @param blocks Spatial exchangeability blocks.
#' @param nperm Number of permutations.
#' @param seed First non-negative permutation seed.
#' @param Kh,n Positive Hill-function parameters.
#' @param finite_correction Whether to use the plus-one Monte Carlo correction.
#' @param retain_null_moments Whether to retain null means and variances.
#' @param tail Empirical upper-tail convention.
#' @param max_score_entries Maximum dense score-workspace entries.
#' @param verbose Whether to report permutation progress.
#' @export
permutation_cellchat_group_support <- function(
    prepared, group, group_support, components, blocks,
    nperm = 100L, seed = 1L, Kh = 0.5, n = 1,
    finite_correction = FALSE, retain_null_moments = FALSE,
    tail = c("strict_greater", "greater_equal"),
    max_score_entries = 1e8,
    verbose = interactive()) {
  if (!inherits(prepared, "SpatialESSPreparedTriMean") ||
      !inherits(group_support, "SpatialESSGroupSupport") ||
      !inherits(components, "SpatialESSLRComponents") ||
      !inherits(blocks, "SpatialESSBlocks")) {
    stop("Invalid prepared, group-support, LR-component or block object.",
         call. = FALSE)
  }
  if (!identical(prepared$cell_id, blocks$cell_id)) {
    stop("Prepared expression and spatial blocks have different cells/order.",
         call. = FALSE)
  }
  if (any(!components$genes %in% prepared$genes)) {
    stop("Prepared triMean kernel is missing referenced LR genes.", call. = FALSE)
  }
  nperm <- as.integer(nperm)
  seed <- as.integer(seed)
  tail <- match.arg(tail)
  if (length(nperm) != 1L || is.na(nperm) || nperm < 1L ||
      length(seed) != 1L || is.na(seed) ||
      as.double(seed) + as.double(nperm) - 1 > .Machine$integer.max) {
    stop("nperm and seed must define a valid positive integer seed sequence.",
         call. = FALSE)
  }
  info <- .spatialess_group_info(group)
  if (length(info$code) != length(prepared$cell_id) ||
      info$count != length(group_support$group_levels)) {
    stop("Target groups do not align with prepared cells/group support.",
         call. = FALSE)
  }
  group_levels <- group_support$group_levels

  pair_key <- blocks$block_code +
    (as.double(info$code) - 1) * length(blocks$block_levels)
  first_pair <- !duplicated(pair_key)
  groups_per_block <- tabulate(
    blocks$block_code[first_pair], nbins = length(blocks$block_levels)
  )
  exchangeable <- groups_per_block[blocks$block_code] > 1L
  exchangeable_fraction <- mean(exchangeable)
  if (!any(exchangeable)) {
    stop("No spatial block spans multiple target groups; this null cannot change group summaries. Increase block_size.",
         call. = FALSE)
  }

  max_score_entries <- as.double(max_score_entries)
  score_entries <- as.double(nrow(group_support$group_pairs)) * nrow(components$lr)
  if (!is.finite(max_score_entries) || max_score_entries < 1 ||
      score_entries > max_score_entries) {
    stop(sprintf("Score matrix requires %.0f entries; max_score_entries is %.0f.",
                 score_entries, max_score_entries), call. = FALSE)
  }

  observed_expression <- summarize_prepared_trimean(
    prepared, info$code, group_levels = group_levels
  )
  observed <- .spatialess_group_score_matrix(
    observed_expression$average, group_support, components, Kh, n
  )
  reject <- matrix(0L, nrow = nrow(observed), ncol = ncol(observed))
  if (retain_null_moments) {
    null_sum <- matrix(0, nrow = nrow(observed), ncol = ncol(observed))
    null_sum_squares <- matrix(0, nrow = nrow(observed), ncol = ncol(observed))
  }
  moved_fraction <- numeric(nperm)
  permutation_seconds <- numeric(nperm)

  for (b in seq_len(nperm)) {
    start <- proc.time()[["elapsed"]]
    source_for_target <- permute_within_spatial_blocks(
      blocks, seed = seed + b - 1L
    )
    moved_fraction[b] <- mean(info$code[source_for_target] != info$code)
    null_expression <- summarize_prepared_trimean(
      prepared, info$code, source_for_target, group_levels
    )
    null_score <- .spatialess_group_score_matrix(
      null_expression$average, group_support, components, Kh, n
    )
    reject <- reject + if (tail == "strict_greater") {
      null_score > observed
    } else {
      null_score >= observed
    }
    if (retain_null_moments) {
      null_sum <- null_sum + null_score
      null_sum_squares <- null_sum_squares + null_score * null_score
    }
    permutation_seconds[b] <- proc.time()[["elapsed"]] - start
    if (verbose && (b == 1L || b == nperm || b %% 10L == 0L)) {
      message(sprintf("Spatial permutation %d/%d", b, nperm))
    }
  }

  pvalue <- if (finite_correction) {
    (reject + 1) / (nperm + 1)
  } else {
    reject / nperm
  }
  pvalue[observed == 0] <- 1
  active <- which(observed > 0, arr.ind = TRUE)
  support <- group_support$group_pairs
  lr <- components$lr
  interaction <- if ("interaction_name" %in% colnames(lr)) {
    as.character(lr$interaction_name)
  } else {
    paste(lr$ligand, lr$receptor, sep = "_")
  }
  output <- data.frame(
    lr_index = active[, "col"],
    interaction_name = interaction[active[, "col"]],
    ligand = as.character(lr$ligand)[active[, "col"]],
    receptor = as.character(lr$receptor)[active[, "col"]],
    support_pair_index = active[, "row"],
    sender_group = support$sender_group[active[, "row"]],
    receiver_group = support$receiver_group[active[, "row"]],
    sender_group_name = group_levels[support$sender_group[active[, "row"]]],
    receiver_group_name = group_levels[support$receiver_group[active[, "row"]]],
    supported_edges = support$supported_edges[active[, "row"]],
    probability = observed[active],
    pvalue = pvalue[active],
    n_reject = reject[active],
    stringsAsFactors = FALSE
  )
  result <- list(
    group_pairs = output,
    parameters = list(
      nperm = nperm, seed = seed, Kh = Kh, n = n,
      finite_correction = isTRUE(finite_correction),
      score_entries = score_entries, max_score_entries = max_score_entries,
      null = "complete cell-profile permutation within sample x compartment x spatial block",
      tail = tail
    ),
    diagnostics = list(
      block_count = length(blocks$block_levels),
      exchangeable_cell_fraction = exchangeable_fraction,
      moved_group_fraction = moved_fraction,
      permutation_seconds = permutation_seconds,
      mean_permutation_seconds = mean(permutation_seconds)
    )
  )
  if (retain_null_moments) {
    null_mean <- null_sum / nperm
    null_variance <- pmax(0, null_sum_squares / nperm - null_mean * null_mean)
    result$null_moments <- list(mean = null_mean, variance = null_variance)
  }
  structure(result, class = "SpatialESSPermutation")
}
