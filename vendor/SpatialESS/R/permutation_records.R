.spatialess_group_score_records <- function(
    average, group_support, components, record_pairs, Kh, n) {
  indices <- lapply(components[c(
    "ligand", "receptor", "co_a", "co_i", "agonist", "antagonist"
  )], .spatialess_index_matrix, genes = rownames(average))
  support <- group_support$group_pairs
  support_index <- record_pairs$support_pair_index
  score_cellchat_group_records_cpp(
    average,
    support$sender_group[support_index],
    support$receiver_group[support_index],
    record_pairs$lr_index,
    indices$ligand, indices$receptor, indices$co_a, indices$co_i,
    indices$agonist, indices$antagonist,
    components$has_agonist, components$has_antagonist, Kh, n
  )
}

#' Spatial block permutation for selected LR/group-support records
#'
#' This record-sliced reference evaluates only explicitly requested
#' `(support_pair_index, lr_index)` combinations. It uses the same profile
#' permutation, exact sparse triMean, CellChat molecular formula and empirical
#' tail conventions as [permutation_cellchat_group_support()].
#'
#' @param prepared A prepared sparse triMean kernel.
#' @param group Fixed target-group labels.
#' @param group_support Spatially supported group pairs.
#' @param components Resolved CellChat LR components.
#' @param blocks Spatial exchangeability blocks.
#' @param record_pairs Data frame with one-based `support_pair_index` and
#'   `lr_index` columns. Duplicate records are not allowed.
#' @param nperm Number of permutations.
#' @param seed First non-negative permutation seed.
#' @param Kh,n Positive Hill-function parameters.
#' @param finite_correction Whether to use the plus-one correction.
#' @param tail Empirical upper-tail convention.
#' @param max_records Maximum selected records.
#' @param verbose Whether to report permutation progress.
#' @export
permutation_cellchat_records <- function(
    prepared, group, group_support, components, blocks, record_pairs,
    nperm = 100L, seed = 1L, Kh = 0.5, n = 1,
    finite_correction = FALSE,
    tail = c("strict_greater", "greater_equal"),
    max_records = 1e7, verbose = interactive()) {
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
    stop("Prepared triMean kernel is missing referenced LR genes.",
         call. = FALSE)
  }
  if (!is.data.frame(record_pairs) ||
      !all(c("support_pair_index", "lr_index") %in% colnames(record_pairs))) {
    stop("record_pairs must contain support_pair_index and lr_index columns.",
         call. = FALSE)
  }
  record_pairs <- data.frame(
    support_pair_index = as.integer(record_pairs$support_pair_index),
    lr_index = as.integer(record_pairs$lr_index)
  )
  if (!nrow(record_pairs) || anyNA(record_pairs) ||
      any(record_pairs$support_pair_index < 1L) ||
      any(record_pairs$support_pair_index > nrow(group_support$group_pairs)) ||
      any(record_pairs$lr_index < 1L) ||
      any(record_pairs$lr_index > nrow(components$lr)) ||
      anyDuplicated(record_pairs)) {
    stop("record_pairs contains invalid or duplicate record indices.",
         call. = FALSE)
  }
  max_records <- as.double(max_records)
  if (!is.finite(max_records) || max_records < 1 ||
      nrow(record_pairs) > max_records) {
    stop("Selected record count exceeds max_records.", call. = FALSE)
  }
  nperm <- as.integer(nperm)
  seed <- as.integer(seed)
  tail <- match.arg(tail)
  if (length(nperm) != 1L || is.na(nperm) || nperm < 1L ||
      length(seed) != 1L || is.na(seed) ||
      as.double(seed) + nperm - 1 > .Machine$integer.max) {
    stop("nperm and seed must define a valid positive integer seed sequence.",
         call. = FALSE)
  }
  info <- .spatialess_group_info(group)
  if (length(info$code) != length(prepared$cell_id) ||
      info$count != length(group_support$group_levels)) {
    stop("Target groups do not align with prepared cells/group support.",
         call. = FALSE)
  }
  pair_key <- blocks$block_code +
    (as.double(info$code) - 1) * length(blocks$block_levels)
  first_pair <- !duplicated(pair_key)
  groups_per_block <- tabulate(
    blocks$block_code[first_pair], nbins = length(blocks$block_levels)
  )
  exchangeable <- groups_per_block[blocks$block_code] > 1L
  if (!any(exchangeable)) {
    stop("No spatial block spans multiple target groups; increase block_size.",
         call. = FALSE)
  }
  group_levels <- group_support$group_levels
  observed_expression <- summarize_prepared_trimean(
    prepared, info$code, group_levels = group_levels
  )
  observed <- .spatialess_group_score_records(
    observed_expression$average, group_support, components,
    record_pairs, Kh, n
  )
  reject <- integer(length(observed))
  moved_fraction <- numeric(nperm)
  permutation_seconds <- numeric(nperm)
  for (permutation_index in seq_len(nperm)) {
    start <- proc.time()[["elapsed"]]
    source_for_target <- permute_within_spatial_blocks(
      blocks, seed = seed + permutation_index - 1L
    )
    moved_fraction[[permutation_index]] <- mean(
      info$code[source_for_target] != info$code
    )
    null_expression <- summarize_prepared_trimean(
      prepared, info$code, source_for_target, group_levels
    )
    null_score <- .spatialess_group_score_records(
      null_expression$average, group_support, components,
      record_pairs, Kh, n
    )
    reject <- reject + if (tail == "strict_greater") {
      null_score > observed
    } else {
      null_score >= observed
    }
    permutation_seconds[[permutation_index]] <-
      proc.time()[["elapsed"]] - start
    if (verbose && (permutation_index == 1L || permutation_index == nperm ||
                    permutation_index %% 10L == 0L)) {
      message(sprintf("Record-sliced spatial permutation %d/%d",
                      permutation_index, nperm))
    }
  }
  pvalue <- if (finite_correction) {
    (reject + 1) / (nperm + 1)
  } else {
    reject / nperm
  }
  pvalue[observed == 0] <- 1
  support <- group_support$group_pairs
  lr <- components$lr
  interaction <- if ("interaction_name" %in% colnames(lr)) {
    as.character(lr$interaction_name)
  } else {
    paste(lr$ligand, lr$receptor, sep = "_")
  }
  support_index <- record_pairs$support_pair_index
  lr_index <- record_pairs$lr_index
  output <- data.frame(
    lr_index = lr_index,
    interaction_name = interaction[lr_index],
    ligand = as.character(lr$ligand)[lr_index],
    receptor = as.character(lr$receptor)[lr_index],
    support_pair_index = support_index,
    sender_group = support$sender_group[support_index],
    receiver_group = support$receiver_group[support_index],
    sender_group_name = group_levels[support$sender_group[support_index]],
    receiver_group_name = group_levels[support$receiver_group[support_index]],
    supported_edges = support$supported_edges[support_index],
    probability = observed, pvalue = pvalue, n_reject = reject,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      group_pairs = output,
      parameters = list(
        nperm = nperm, seed = seed, Kh = Kh, n = n,
        finite_correction = isTRUE(finite_correction), tail = tail,
        mode = "record_sliced_spatial_permutation",
        selected_records = nrow(record_pairs)
      ),
      diagnostics = list(
        block_count = length(blocks$block_levels),
        exchangeable_cell_fraction = mean(exchangeable),
        moved_group_fraction = moved_fraction,
        permutation_seconds = permutation_seconds,
        mean_permutation_seconds = mean(permutation_seconds)
      )
    ),
    class = "SpatialESSRecordPermutation"
  )
}
