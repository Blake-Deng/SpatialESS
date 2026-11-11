#' Prepare count-scale Visium HD bins
#'
#' Aligns spatial positions by barcode and expresses bin centres in micrometres
#' on the capture-array grid. Image pixel coordinates are never used as distances.
#' Duplicate gene names are summed without changing gene identifiers.
#'
#' @param counts Non-negative integer gene-by-bin count matrix with dimnames.
#' @param positions Data frame containing `barcode`, zero-based `array_row`,
#'   `array_col` and `in_tissue` (zero or one). A missing tissue mask is only
#'   permitted when `in_tissue = FALSE` and is recorded in provenance.
#' @param bin_size One positive even integer bin width in micrometres.
#' @param in_tissue Keep only bins with `in_tissue == 1`.
#' @return A `SpatialESSVisiumHD` list with sparse `counts`, bin-aligned
#'   `coordinates` in micrometres, `metadata`, bin sizes and `provenance`.
#' @export
prepare_visium_hd <- function(counts, positions, bin_size, in_tissue = TRUE) {
  bin_size <- .visium_hd_bin_size(bin_size)
  .visium_hd_flag(in_tissue, "in_tissue")
  counts <- .visium_hd_counts(counts)
  position_info <- .visium_hd_positions(positions, bin_size, require_tissue = in_tissue)
  all_positions <- position_info$positions
  index <- match(colnames(counts), all_positions$barcode)
  if (anyNA(index)) stop("Every count-matrix barcode requires a spatial position.", call. = FALSE)
  positions <- all_positions[index, , drop = FALSE]
  n_input <- ncol(counts)
  if (in_tissue) {
    keep <- positions$in_tissue == 1L
    counts <- counts[, keep, drop = FALSE]
    positions <- positions[keep, , drop = FALSE]
  }
  if (!ncol(counts)) stop("No bins remain after tissue filtering.", call. = FALSE)
  coordinates <- cbind(x = (positions$array_col + 0.5) * bin_size,
                       y = (positions$array_row + 0.5) * bin_size)
  rownames(coordinates) <- colnames(counts)
  structure(list(counts = counts, coordinates = coordinates, metadata = positions,
                 all_positions = all_positions,
                 bin_size = bin_size, source_bin_size = bin_size,
                 provenance = list(coordinate_system = "array_grid_bin_centres_um",
                                   tissue_filtered = in_tissue, input_bins = n_input,
                                   counts_removed_by_tissue_filter = n_input - ncol(counts),
                                   tissue_mask_available = position_info$tissue_mask_available,
                                   tissue_rule = if (position_info$tissue_mask_available)
                                     "provided_source_mask" else "assumed_all_tissue_no_mask",
                                   position_universe = "all supplied positions; slide-design completeness not independently verified",
                                   retained_bins = ncol(counts), aggregation = NULL)),
            class = "SpatialESSVisiumHD")
}

.visium_hd_positions <- function(positions, bin_size, require_tissue = FALSE) {
  positions <- as.data.frame(positions, stringsAsFactors = FALSE)
  required <- c("barcode", "array_row", "array_col")
  if (!all(required %in% names(positions))) {
    stop("positions requires barcode, array_row and array_col columns.", call. = FALSE)
  }
  positions$barcode <- as.character(positions$barcode)
  if (anyNA(positions$barcode) || any(!nzchar(positions$barcode)) || anyDuplicated(positions$barcode)) {
    stop("positions barcodes must be non-empty and unique.", call. = FALSE)
  }
  for (key in c("array_row", "array_col")) {
    value <- positions[[key]]
    if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) || any(value < 0 | value != floor(value))) {
      stop("array_row and array_col must be finite non-negative integers.", call. = FALSE)
    }
  }
  tissue_mask_available <- "in_tissue" %in% names(positions)
  if (!tissue_mask_available && require_tissue) {
    stop("positions requires in_tissue when tissue filtering is requested.", call. = FALSE)
  }
  if (!tissue_mask_available) positions$in_tissue <- 1L
  if (anyNA(positions$in_tissue) || any(!positions$in_tissue %in% c(0, 1))) {
    stop("positions in_tissue must contain only zero or one.", call. = FALSE)
  }
  positions$in_tissue <- as.integer(positions$in_tissue)
  suffix <- .visium_hd_suffix(positions$barcode)
  grid_key <- paste(positions$array_row, positions$array_col, suffix, sep = ":")
  if (anyDuplicated(grid_key)) {
    stop("Spatial positions contain duplicate array-grid bins within a barcode suffix.", call. = FALSE)
  }
  known <- grepl("^s_[0-9]+um_[0-9]+_[0-9]+-[0-9]+$", positions$barcode)
  if (any(known)) {
    encoded_size <- as.numeric(sub("^s_([0-9]+)um_.*$", "\\1", positions$barcode[known]))
    if (any(encoded_size != bin_size)) {
      stop("bin_size differs from the resolution encoded in the barcodes.", call. = FALSE)
    }
    encoded_row <- as.numeric(sub("^s_[0-9]+um_([0-9]+)_[0-9]+-[0-9]+$", "\\1", positions$barcode[known]))
    encoded_col <- as.numeric(sub("^s_[0-9]+um_[0-9]+_([0-9]+)-[0-9]+$", "\\1", positions$barcode[known]))
    if (any(encoded_row != positions$array_row[known] | encoded_col != positions$array_col[known])) {
      stop("Array positions differ from the grid indices encoded in the barcodes.", call. = FALSE)
    }
  }
  positions$bin_size_um <- bin_size
  rownames(positions) <- positions$barcode
  list(positions = positions, tissue_mask_available = tissue_mask_available)
}

#' Sum Visium HD counts into larger square bins
#'
#' Uses the zero-based capture-array grid and sparse matrix multiplication.
#' Counts are summed before any normalization or deconvolution. A target bin is
#' in tissue if any supplied source position is in tissue, including positions
#' without an observed count-matrix column. The result
#' contains only source bins present in `hd`; filtered input cannot recover
#' counts from omitted bins at tissue boundaries.
#'
#' @param hd A `SpatialESSVisiumHD` object.
#' @param bin_size Target width in micrometres, an integer multiple of the source.
#' @return A `SpatialESSVisiumHD` with the summed counts and target-bin centres.
#' @export
aggregate_visium_hd <- function(hd, bin_size) {
  .visium_hd_validate(hd)
  bin_size <- .visium_hd_bin_size(bin_size)
  factor <- bin_size / hd$bin_size
  if (factor < 1 || factor != floor(factor)) {
    stop("Target bin_size must be an integer multiple of the source bin_size; finer bins cannot be reconstructed.", call. = FALSE)
  }
  if (factor == 1) return(hd)
  if (is.null(hd$all_positions)) {
    stop("hd must retain all_positions from the current Visium HD reader/preparer before aggregation.", call. = FALSE)
  }
  source_positions <- .visium_hd_positions(hd$all_positions, hd$bin_size)$positions
  row <- floor(source_positions$array_row / factor)
  col <- floor(source_positions$array_col / factor)
  suffix <- .visium_hd_suffix(source_positions$barcode)
  barcode <- sprintf("s_%03dum_%05d_%05d-%s", bin_size, row, col, suffix)
  first <- !duplicated(barcode)
  order <- order(suffix[first], row[first], col[first])
  all_target_ids <- barcode[first][order]
  position_index <- match(barcode, all_target_ids)
  position_assignment <- Matrix::sparseMatrix(
    i = seq_along(position_index), j = position_index, x = 1,
    dims = c(length(position_index), length(all_target_ids)))
  # Count columns and tissue positions have different universes: Space Ranger
  # raw matrices omit zero-read barcodes, while tissue_positions contains them.
  source_index <- match(colnames(hd$counts), source_positions$barcode)
  if (anyNA(source_index)) stop("Count bins are missing from all_positions.", call. = FALSE)
  count_target_barcodes <- barcode[source_index]
  target_ids <- all_target_ids[all_target_ids %in% count_target_barcodes]
  index <- match(count_target_barcodes, target_ids)
  assignment <- Matrix::sparseMatrix(i = seq_along(index), j = index, x = 1,
                                     dims = c(length(index), length(target_ids)),
                                     dimnames = list(colnames(hd$counts), target_ids))
  counts <- methods::as(hd$counts %*% assignment, "dgCMatrix")
  representatives <- match(all_target_ids, barcode)
  tissue <- as.numeric(Matrix::crossprod(position_assignment, source_positions$in_tissue)) > 0
  position_counts <- as.integer(Matrix::colSums(position_assignment))
  positions <- data.frame(barcode = all_target_ids, array_row = row[representatives],
                          array_col = col[representatives], in_tissue = as.integer(tissue),
                          source_bin_count = tabulate(match(count_target_barcodes, all_target_ids),
                                                      nbins = length(all_target_ids)),
                          source_position_count = position_counts,
                          expected_source_positions = factor^2,
                          source_positions_complete = position_counts == factor^2,
                          stringsAsFactors = FALSE)
  result <- prepare_visium_hd(counts, positions, bin_size, in_tissue = FALSE)
  result$source_bin_size <- hd$source_bin_size
  result$provenance <- hd$provenance
  potential_partial <- isTRUE(hd$provenance$incomplete_count_source) ||
    isTRUE(hd$provenance$potential_incomplete_source) ||
    identical(hd$provenance$matrix_kind, "filtered_feature_bc_matrix") ||
    isTRUE(hd$provenance$counts_removed_by_tissue_filter > 0)
  result$provenance$potential_incomplete_source <- potential_partial
  result$provenance$tissue_rule <- if (isTRUE(hd$provenance$tissue_mask_available))
    "any_supplied_source_position_in_tissue" else "assumed_all_tissue_no_mask"
  result$provenance$retained_bins <- ncol(counts)
  result$provenance$aggregation <- c(hd$provenance$aggregation, list(list(
    from_bin_size = hd$bin_size, to_bin_size = bin_size,
    source_bins = ncol(hd$counts), target_bins = ncol(counts),
    expected_sources_per_complete_bin = factor^2,
    position_universe = "all supplied source positions, including zero-read barcodes",
    count_scope = "sum of supplied source count columns; omitted columns cannot be recovered",
    potential_incomplete_source = potential_partial,
    fully_populated_target_position_grids = sum(positions$source_positions_complete),
    observed_target_position_grids = nrow(positions),
    native_output_equivalence = "not empirically verified",
    source_counts = "sum of all available source bins")))
  result
}

#' Read Space Ranger Visium HD outputs at a chosen bin size
#'
#' Accepts a run directory, `outs`, `binned_outputs`, or one `square_NNNum`
#' directory. Native requested bins are used when present; otherwise the largest
#' available finer divisor is summed. An explicitly selected square directory
#' requires `source_bin_size` when rebinning to another resolution. For rebinning,
#' raw counts are required by default and tissue filtering follows summation.
#' The filtered feature policy applies the source filtered feature IDs to raw
#' counts before gene-name mapping and aggregation, preserving probe exclusions.
#'
#' @param data_dir Space Ranger output directory.
#' @param bin_size Requested width in micrometres, for example 8, 16 or 32.
#' @param source_bin_size Optional explicit source resolution in micrometres.
#' @param in_tissue Keep bins marked in tissue after aggregation.
#' @param gene_column Select feature names or feature IDs. Duplicate values are
#'   summed. Feature types other than Gene Expression are excluded.
#' @param positions_file Optional CSV or Parquet override for the selected
#'   source resolution. CSV supports a header or the six Space Ranger columns.
#' @param feature_policy `"filtered"` (default) uses the Gene Expression feature
#'   IDs in the source filtered matrix metadata at every resolution; `"raw"`
#'   uses all raw Gene Expression features at every resolution.
#' @param allow_filtered_aggregation Explicitly permit partial-count coarsening
#'   if raw input is unavailable. This emits a warning and records incompleteness.
#' @return A `SpatialESSVisiumHD` object. Parquet requires optional `arrow`;
#'   HDF5-only input requires optional `hdf5r`. MEX plus CSV uses base R and Matrix.
#' @export
read_visium_hd <- function(data_dir, bin_size = 8, source_bin_size = NULL,
                           in_tissue = TRUE, gene_column = c("name", "id"), positions_file = NULL,
                           feature_policy = c("filtered", "raw"),
                           allow_filtered_aggregation = FALSE) {
  bin_size <- .visium_hd_bin_size(bin_size)
  .visium_hd_flag(in_tissue, "in_tissue")
  .visium_hd_flag(allow_filtered_aggregation, "allow_filtered_aggregation")
  gene_column <- match.arg(gene_column)
  feature_policy <- match.arg(feature_policy)
  if (length(data_dir) != 1L || is.na(data_dir) || !dir.exists(data_dir)) {
    stop("data_dir must be an existing directory.", call. = FALSE)
  }
  if (!is.null(source_bin_size)) source_bin_size <- .visium_hd_bin_size(source_bin_size)
  selected <- .visium_hd_directory(data_dir, bin_size, source_bin_size,
                                    feature_policy, allow_filtered_aggregation)
  source_size <- selected$bin_size
  if (source_size > bin_size || bin_size %% source_size != 0) {
    stop("Requested bin_size must be an integer multiple of source_bin_size.", call. = FALSE)
  }
  rebin <- source_size != bin_size
  raw_source <- .visium_hd_matrix_source(selected$path, "raw_feature_bc_matrix")
  filtered_source <- .visium_hd_matrix_source(selected$path, "filtered_feature_bc_matrix")
  filtered_metadata <- .visium_hd_feature_source(selected$path, "filtered_feature_bc_matrix")
  keep_ids <- NULL
  feature_id_source <- NULL
  if (feature_policy == "filtered") {
    if (is.null(filtered_metadata)) {
      stop("feature_policy='filtered' requires source filtered feature metadata to preserve target/probe exclusions; provide filtered features.tsv/H5 metadata or explicitly choose feature_policy='raw'.",
           call. = FALSE)
    }
    allowed <- .visium_hd_read_features(filtered_metadata)
    keep_ids <- allowed$id[allowed$type == "Gene Expression"]
    if (!length(keep_ids)) stop("Filtered feature metadata contains no Gene Expression features.", call. = FALSE)
    feature_id_source <- normalizePath(filtered_metadata$path)
  }
  source <- if (feature_policy == "raw" || rebin) raw_source else filtered_source
  if (is.null(source) && feature_policy == "filtered" && !rebin) source <- raw_source
  incomplete_count_source <- FALSE
  if (is.null(source) && rebin && !is.null(filtered_source) &&
      feature_policy == "filtered" && allow_filtered_aggregation) {
    source <- filtered_source
    incomplete_count_source <- TRUE
    warning("Filtered-only aggregation uses incomplete source counts; missing out-of-tissue counts cannot be recovered. This is not guaranteed equivalent to native Space Ranger bins.",
            call. = FALSE)
  }
  if (is.null(source)) {
    stop(if (rebin || feature_policy == "raw")
      "Raw source counts are required for this resolution/policy; filtered-only aggregation requires explicit allow_filtered_aggregation=TRUE with feature_policy='filtered'."
      else "No usable filtered or raw count matrix found in the selected bin directory.",
      call. = FALSE)
  }
  counts <- if (source$format == "mex")
    .visium_hd_read_mex(source$path, gene_column, keep_ids) else
      .visium_hd_read_h5(source$path, gene_column, keep_ids)
  kept_ids <- attr(counts, "visium_hd_feature_ids", exact = TRUE)
  attr(counts, "visium_hd_feature_ids") <- NULL
  if (is.null(feature_id_source)) feature_id_source <- normalizePath(source$path)
  if (is.null(positions_file)) {
    candidates <- file.path(selected$path, "spatial", c("tissue_positions.parquet", "tissue_positions.csv", "tissue_positions_list.csv"))
    positions_file <- candidates[file.exists(candidates)][1L]
  }
  if (length(positions_file) != 1L || is.na(positions_file) || !file.exists(positions_file)) {
    stop("No tissue positions found; provide positions_file for the selected source bin size.", call. = FALSE)
  }
  positions <- .visium_hd_read_positions(positions_file)
  tissue_mask_available <- "in_tissue" %in% names(positions)
  if (in_tissue && !tissue_mask_available) {
    stop("positions requires in_tissue when tissue filtering is requested.", call. = FALSE)
  }
  hd <- prepare_visium_hd(counts, positions, source_size, in_tissue = FALSE)
  hd$provenance$data_dir <- normalizePath(data_dir)
  hd$provenance$source_directory <- normalizePath(selected$path)
  hd$provenance$matrix_source <- normalizePath(source$path)
  hd$provenance$matrix_kind <- source$kind
  hd$provenance$feature_policy <- feature_policy
  hd$provenance$source_feature_id_source <- feature_id_source
  hd$provenance$feature_ids_kept <- kept_ids
  hd$provenance$feature_id_count <- length(kept_ids)
  hd$provenance$feature_filter_order <- "feature ID selection before gene-name mapping and duplicate-symbol summation"
  hd$provenance$incomplete_count_source <- incomplete_count_source
  hd$provenance$tissue_mask_available <- tissue_mask_available
  hd$provenance$positions_file <- normalizePath(positions_file)
  hd$provenance$gene_column <- gene_column
  if (rebin) hd <- aggregate_visium_hd(hd, bin_size)
  n_before_tissue_filter <- ncol(hd$counts)
  if (in_tissue) {
    keep <- hd$metadata$in_tissue == 1L
    if (!any(keep)) stop("No bins remain after tissue filtering.", call. = FALSE)
    hd$counts <- hd$counts[, keep, drop = FALSE]
    hd$coordinates <- hd$coordinates[keep, , drop = FALSE]
    hd$metadata <- hd$metadata[keep, , drop = FALSE]
  }
  hd$provenance$tissue_filtered <- in_tissue
  hd$provenance$counts_removed_by_tissue_filter <- n_before_tissue_filter - ncol(hd$counts)
  hd$provenance$native_output_equivalence <- "not empirically verified"
  hd$provenance$retained_bins <- ncol(hd$counts)
  hd
}

.visium_hd_bin_size <- function(value) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) || !is.finite(value) ||
      value < 2 || value %% 2 != 0 || value > .Machine$integer.max) {
    stop("bin_size must be one positive even integer width in micrometres.", call. = FALSE)
  }
  as.numeric(value)
}
.visium_hd_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) stop(paste0(name, " must be TRUE or FALSE."), call. = FALSE)
}
.visium_hd_suffix <- function(barcode) {
  known <- grepl("-[0-9]+$", barcode)
  suffix <- rep("1", length(barcode))
  suffix[known] <- sub("^.*-([0-9]+)$", "\\1", barcode[known])
  suffix
}
.visium_hd_counts <- function(counts) {
  if (!(is.matrix(counts) || inherits(counts, "Matrix")) || length(dim(counts)) != 2L || any(dim(counts) == 0L)) {
    stop("counts must be a non-empty gene-by-bin matrix.", call. = FALSE)
  }
  genes <- rownames(counts)
  barcodes <- colnames(counts)
  if (is.null(genes) || is.null(barcodes) || anyNA(genes) || anyNA(barcodes) || any(!nzchar(genes)) || any(!nzchar(barcodes)) || anyDuplicated(barcodes)) {
    stop("counts requires non-empty gene names and unique non-empty barcodes.", call. = FALSE)
  }
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  counts <- methods::as(methods::as(counts, "generalMatrix"), "CsparseMatrix")
  counts <- methods::as(counts, "dMatrix")
  if (any(!is.finite(counts@x)) || any(counts@x < 0) || any(counts@x != floor(counts@x))) {
    stop("counts must contain finite non-negative integer raw counts; do not normalize or round input.", call. = FALSE)
  }
  if (anyDuplicated(genes)) {
    unique_genes <- unique(genes)
    collapse <- Matrix::sparseMatrix(i = match(genes, unique_genes), j = seq_along(genes), x = 1,
                                     dims = c(length(unique_genes), length(genes)),
                                     dimnames = list(unique_genes, genes))
    counts <- methods::as(collapse %*% counts, "dgCMatrix")
  }
  Matrix::drop0(counts)
}
.visium_hd_validate <- function(hd) {
  if (!inherits(hd, "SpatialESSVisiumHD") || !inherits(hd$counts, "dgCMatrix") ||
      !identical(colnames(hd$counts), rownames(hd$coordinates)) ||
      !identical(colnames(hd$counts), rownames(hd$metadata)) || nrow(hd$metadata) != ncol(hd$counts) ||
      !all(c("barcode", "array_row", "array_col", "in_tissue") %in% names(hd$metadata))) {
    stop("hd must be an aligned SpatialESSVisiumHD object from prepare_visium_hd() or read_visium_hd().", call. = FALSE)
  }
  invisible(TRUE)
}
.visium_hd_directory <- function(data_dir, bin_size, source_bin_size,
                                   feature_policy = "filtered",
                                   allow_filtered_aggregation = FALSE) {
  data_dir <- normalizePath(data_dir)
  base <- basename(data_dir)
  direct <- grepl("^square_[0-9]+um$", base)
  if (direct) {
    native <- as.numeric(sub("^square_([0-9]+)um$", "\\1", base))
    if (is.null(source_bin_size) && bin_size != native) {
      stop("A direct square directory differs from bin_size; set source_bin_size explicitly to request aggregation.", call. = FALSE)
    }
    if (!is.null(source_bin_size) && source_bin_size != native) stop("source_bin_size differs from the selected square directory.", call. = FALSE)
    return(list(path = data_dir, bin_size = native))
  }
  roots <- unique(c(data_dir, file.path(data_dir, "binned_outputs"), file.path(data_dir, "outs", "binned_outputs")))
  choices <- unlist(lapply(roots, function(root) {
    if (!dir.exists(root)) return(character())
    dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
    dirs[grepl("^square_[0-9]+um$", basename(dirs))]
  }), use.names = FALSE)
  choices <- unique(choices)
  if (!length(choices)) {
    if (!is.null(source_bin_size) && any(file.exists(file.path(data_dir, c("filtered_feature_bc_matrix", "filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix", "raw_feature_bc_matrix.h5"))))) {
      return(list(path = data_dir, bin_size = source_bin_size))
    }
    stop("No square_NNNum bin directories found (use source_bin_size for a renamed bin directory).", call. = FALSE)
  }
  sizes <- as.numeric(sub("^square_([0-9]+)um$", "\\1", basename(choices)))
  target <- if (is.null(source_bin_size)) {
    valid <- sizes <= bin_size & bin_size %% sizes == 0
    if (!any(valid)) stop("No available source resolution divides the requested bin_size.", call. = FALSE)
    has_raw <- vapply(choices, function(path)
      !is.null(.visium_hd_matrix_source(path, "raw_feature_bc_matrix")), logical(1))
    has_filtered_features <- vapply(choices, function(path)
      !is.null(.visium_hd_feature_source(path, "filtered_feature_bc_matrix")), logical(1))
    has_filtered <- vapply(choices, function(path)
      !is.null(.visium_hd_matrix_source(path, "filtered_feature_bc_matrix")), logical(1))
    native_policy <- if (feature_policy == "raw") has_raw else
      has_filtered | (has_raw & has_filtered_features)
    native_usable <- valid & sizes == bin_size & native_policy
    if (any(native_usable)) {
      bin_size
    } else {
      usable <- valid & sizes < bin_size & has_raw &
        (feature_policy == "raw" | has_filtered_features)
      if (any(usable)) {
        max(sizes[usable])
      } else if (allow_filtered_aggregation && feature_policy == "filtered") {
        partial <- valid & has_filtered & has_filtered_features
        if (!any(partial)) stop("No usable raw source or filtered source with feature metadata is available for aggregation.", call. = FALSE)
        max(sizes[partial])
      } else {
        stop("No usable raw source is available for aggregation with this feature policy. Default filtered policy also requires filtered feature metadata; provide these files, explicitly select feature_policy='raw', or opt into incomplete filtered aggregation.",
             call. = FALSE)
      }
    }
  } else source_bin_size
  match <- which(sizes == target)
  if (!length(match)) stop("Requested source_bin_size directory is unavailable.", call. = FALSE)
  if (length(match) > 1L) stop("Multiple matching bin directories found; select binned_outputs or a square directory explicitly.", call. = FALSE)
  list(path = choices[match], bin_size = target)
}
.visium_hd_matrix_source <- function(directory, kind) {
  mex <- file.path(directory, kind)
  if (dir.exists(mex) && any(file.exists(file.path(mex, c("matrix.mtx.gz", "matrix.mtx"))))) {
    return(list(path = mex, format = "mex", kind = kind))
  }
  h5 <- file.path(directory, paste0(kind, ".h5"))
  if (file.exists(h5)) return(list(path = h5, format = "h5", kind = kind))
  NULL
}
.visium_hd_feature_source <- function(directory, kind) {
  mex <- file.path(directory, kind)
  if (dir.exists(mex) && any(file.exists(file.path(mex,
    c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))))) {
    return(list(path = mex, format = "mex", kind = kind))
  }
  h5 <- file.path(directory, paste0(kind, ".h5"))
  if (file.exists(h5)) return(list(path = h5, format = "h5", kind = kind))
  NULL
}
.visium_hd_existing <- function(directory, names) {
  paths <- file.path(directory, names)
  existing <- paths[file.exists(paths)]
  if (!length(existing)) stop(paste0("Missing input file: ", paste(names, collapse = " or "), "."), call. = FALSE)
  existing[1L]
}
.visium_hd_read_table <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  utils::read.delim(connection, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE, na.strings = character())
}
.visium_hd_normalize_features <- function(features) {
  features <- as.data.frame(features, stringsAsFactors = FALSE)
  if (!nrow(features) || !ncol(features)) stop("Feature metadata must not be empty.", call. = FALSE)
  result <- data.frame(
    id = as.character(features[[1L]]),
    name = as.character(features[[if (ncol(features) >= 2L) 2L else 1L]]),
    type = if (ncol(features) >= 3L) as.character(features[[3L]]) else "Gene Expression",
    stringsAsFactors = FALSE
  )
  if (anyNA(result$id) || any(!nzchar(result$id)) || anyDuplicated(result$id) ||
      anyNA(result$type) || any(!nzchar(result$type))) {
    stop("Feature IDs must be unique/nonempty and feature types must be present.", call. = FALSE)
  }
  result
}
.visium_hd_read_features <- function(source) {
  if (source$format == "mex") {
    file <- .visium_hd_existing(source$path,
      c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))
    return(.visium_hd_normalize_features(.visium_hd_read_table(file)))
  }
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("HDF5 feature metadata requires optional package 'hdf5r'.", call. = FALSE)
  }
  h5 <- hdf5r::H5File$new(source$path, mode = "r")
  on.exit(h5$close_all())
  .visium_hd_h5_features(h5)
}
.visium_hd_h5_features <- function(h5) {
  if (!h5$exists("matrix/features/id") || !h5$exists("matrix/features/name")) {
    stop("Expected 10x HDF5 matrix/features/id and name metadata.", call. = FALSE)
  }
  features <- data.frame(id = as.character(h5[["matrix/features/id"]]$read()),
                         name = as.character(h5[["matrix/features/name"]]$read()),
                         stringsAsFactors = FALSE)
  if (h5$exists("matrix/features/feature_type"))
    features$type <- as.character(h5[["matrix/features/feature_type"]]$read())
  .visium_hd_normalize_features(features)
}
.visium_hd_feature_matrix <- function(counts, features, barcodes, gene_column,
                                       keep_feature_ids = NULL) {
  features <- .visium_hd_normalize_features(features)
  if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) stop("Feature/barcode dimensions differ from the count matrix.", call. = FALSE)
  keep <- features$type == "Gene Expression"
  if (!is.null(keep_feature_ids)) {
    if (length(setdiff(keep_feature_ids, features$id))) {
      stop("Source raw matrix is missing IDs from the filtered feature universe.", call. = FALSE)
    }
    keep <- keep & features$id %in% keep_feature_ids
  }
  if (!any(keep)) stop("No Gene Expression features remain under the selected feature policy.", call. = FALSE)
  counts <- counts[keep, , drop = FALSE]
  features <- features[keep, , drop = FALSE]
  rownames(counts) <- features[[gene_column]]
  colnames(counts) <- as.character(barcodes)
  counts <- .visium_hd_counts(counts)
  attr(counts, "visium_hd_feature_ids") <- features$id
  counts
}
.visium_hd_read_mex <- function(directory, gene_column, keep_feature_ids = NULL) {
  matrix_file <- .visium_hd_existing(directory, c("matrix.mtx.gz", "matrix.mtx"))
  feature_file <- .visium_hd_existing(directory, c("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))
  barcode_file <- .visium_hd_existing(directory, c("barcodes.tsv.gz", "barcodes.tsv"))
  connection <- if (grepl("\\.gz$", matrix_file)) gzfile(matrix_file, "rt") else file(matrix_file, "rt")
  on.exit(close(connection))
  counts <- Matrix::readMM(connection)
  features <- .visium_hd_read_table(feature_file)
  barcodes <- .visium_hd_read_table(barcode_file)[[1L]]
  .visium_hd_feature_matrix(counts, features, barcodes, gene_column, keep_feature_ids)
}
.visium_hd_read_h5 <- function(path, gene_column, keep_feature_ids = NULL) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("HDF5 input requires optional package 'hdf5r'; alternatively provide the MEX matrix directory.", call. = FALSE)
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  if (!h5$exists("matrix")) stop("Expected a 10x HDF5 'matrix' group.", call. = FALSE)
  mat <- h5[["matrix"]]
  counts <- Matrix::sparseMatrix(i = as.integer(mat[["indices"]]$read()), p = as.integer(mat[["indptr"]]$read()),
                                 x = as.numeric(mat[["data"]]$read()), dims = as.integer(mat[["shape"]]$read()), index1 = FALSE)
  features <- .visium_hd_h5_features(h5)
  .visium_hd_feature_matrix(counts, features, mat[["barcodes"]]$read(), gene_column,
                             keep_feature_ids)
}
.visium_hd_read_positions <- function(path) {
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("Parquet positions require optional package 'arrow'; alternatively pass positions_file with a CSV export.", call. = FALSE)
    return(as.data.frame(arrow::read_parquet(path)))
  }
  first <- readLines(path, n = 1L, warn = FALSE)
  has_header <- grepl("(^|,)\"?barcode\"?(,|$)", first)
  positions <- utils::read.csv(path, header = has_header, stringsAsFactors = FALSE, check.names = FALSE)
  if (!has_header) {
    if (ncol(positions) != 6L) stop("Headerless tissue positions require the six Space Ranger columns.", call. = FALSE)
    names(positions) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  }
  positions
}
