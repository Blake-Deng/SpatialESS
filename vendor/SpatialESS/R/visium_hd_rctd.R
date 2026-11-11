#' Deconvolve Visium HD bins using the original spacexr RCTD implementation
#'
#' Fits each resolution independently, learning platform effects and sigma once
#' across the complete resolution before fitting pixels in bounded batches.
#' Untransformed integer UMI counts are required. This interface uses the
#' original dmcable/spacexr API, not the incompatible Bioconductor API.
#'
#' @param hd A SpatialESSVisiumHD object returned by read_visium_hd().
#' @param reference_counts Gene-by-cell raw reference count matrix.
#' @param cell_types Named reference-cell annotations, aligned by barcode.
#' @param doublet_mode RCTD mode: full (default), doublet, or multi.
#' @param max_cores Number of RCTD workers.
#' @param min_umi,min_genes Minimum total UMI and detected genes per input bin.
#' @param create_args Additional named create.RCTD arguments, including marker,
#'   reference-cell and sigma thresholds. UMI_min defaults to min_umi.
#' @param run_args Named pixel-fitting configuration overrides: MAX_MULTI_TYPES,
#'   CONFIDENCE_THRESHOLD, DOUBLET_THRESHOLD, or MIN_CHANGE_REG.
#' @param seed Random seed for reference processing and sigma estimation.
#' @param batch_size Maximum pixels per fitPixels call; NULL fits all at once.
#' @param keep_fit Retain the fitted upstream RCTD object, including raw results.
#' @return A SpatialESSVisiumHDRCTD list with sparse bin-by-type weights,
#'   barcode-aligned annotations, native classification information, diagnostics,
#'   provenance and optional fit. Normalized full-mode weights describe relative
#'   RNA/UMI contributions, not cell-number fractions. In doublet mode,
#'   \code{weights} are derived annotation weights: singlets are one-hot and
#'   doublets use normalized pair weights. Filtered, unfit and rejected bins have
#'   zero rows. Native \code{results_df} (including confidence and convergence)
#'   and unmodified \code{weights_doublet} remain in \code{classification}, even
#'   when \code{keep_fit = FALSE}; these contain only fitted bins. Multi-mode
#'   classification retains native per-bin results except unrestricted
#'   all-type weights. Ambiguous bins and
#'   tied full-mode maxima have no dominant_type. Bins remain spatial bins,
#'   rather than inferred individual cells.
#' @export
run_visium_hd_rctd <- function(
    hd, reference_counts, cell_types,
    doublet_mode = c("full", "doublet", "multi"), max_cores = 1L,
    min_umi = 100, min_genes = 10L, create_args = list(), run_args = list(),
    seed = 1L, batch_size = 5000L, keep_fit = FALSE) {
  doublet_mode <- match.arg(doublet_mode)
  if (!inherits(hd, "SpatialESSVisiumHD")) {
    stop("hd must be a SpatialESSVisiumHD object.", call. = FALSE)
  }
  counts <- .hd_rctd_counts(hd$counts, "hd$counts")
  reference_counts <- .hd_rctd_counts(reference_counts, "reference_counts")
  barcodes <- colnames(counts)
  coordinates <- as.matrix(hd$coordinates)
  if (!is.numeric(coordinates) || !identical(dim(coordinates), c(ncol(counts), 2L)) ||
      any(!is.finite(coordinates)) ||
      (!is.null(rownames(coordinates)) &&
       !identical(rownames(coordinates), barcodes))) {
    stop("hd coordinates must be finite bin-by-2 values aligned with counts.", call. = FALSE)
  }
  rownames(coordinates) <- barcodes
  colnames(coordinates) <- c("x", "y")
  if (!is.null(hd$metadata) &&
      (nrow(hd$metadata) != ncol(counts) ||
       ("barcode" %in% names(hd$metadata) &&
        !identical(as.character(hd$metadata$barcode), barcodes)))) {
    stop("hd metadata must follow the count-matrix barcode order.", call. = FALSE)
  }
  if (is.null(names(cell_types)) || anyNA(names(cell_types)) ||
      anyDuplicated(names(cell_types)) ||
      !setequal(names(cell_types), colnames(reference_counts))) {
    stop("cell_types must have exactly one uniquely named label per reference barcode.", call. = FALSE)
  }
  cell_types <- as.character(cell_types[colnames(reference_counts)])
  if (anyNA(cell_types) || any(!nzchar(trimws(cell_types))) ||
      any(grepl("/", cell_types, fixed = TRUE))) {
    stop("Reference cell types must be nonempty, nonmissing labels without '/'.", call. = FALSE)
  }
  cell_types <- factor(cell_types)
  names(cell_types) <- colnames(reference_counts)
  types <- levels(cell_types)
  if (length(types) < 2L) stop("RCTD requires at least two reference cell types.", call. = FALSE)
  max_cores <- .hd_rctd_number(max_cores, "max_cores", 1, integer = TRUE)
  min_umi <- .hd_rctd_number(min_umi, "min_umi", 0)
  min_genes <- .hd_rctd_number(min_genes, "min_genes", 0, integer = TRUE)
  seed <- .hd_rctd_number(seed, "seed", 0, integer = TRUE)
  if (!is.null(batch_size)) {
    batch_size <- .hd_rctd_number(batch_size, "batch_size", 1, integer = TRUE)
  }
  if (!is.logical(keep_fit) || length(keep_fit) != 1L || is.na(keep_fit)) {
    stop("keep_fit must be TRUE or FALSE.", call. = FALSE)
  }
  .hd_rctd_args(create_args, "create_args")
  .hd_rctd_args(run_args, "run_args")
  reserved <- intersect(names(create_args), c("spatialRNA", "reference", "max_cores"))
  if (length(reserved)) stop("create_args cannot override spatialRNA, reference or max_cores.", call. = FALSE)
  run_allowed <- c("MAX_MULTI_TYPES", "CONFIDENCE_THRESHOLD", "DOUBLET_THRESHOLD", "MIN_CHANGE_REG")
  if (length(setdiff(names(run_args), run_allowed))) {
    stop("run_args only supports MAX_MULTI_TYPES, CONFIDENCE_THRESHOLD, DOUBLET_THRESHOLD and MIN_CHANGE_REG.", call. = FALSE)
  }
  for (key in names(run_args)) {
    .hd_rctd_number(run_args[[key]], paste0("run_args$", key),
                    if (key == "MAX_MULTI_TYPES") 1 else 0,
                    integer = key == "MAX_MULTI_TYPES")
  }

  # Preserve complete-library denominators before matching reference genes.
  nUMI <- Matrix::colSums(counts)
  nUMI_ref <- Matrix::colSums(reference_counts)
  if (any(nUMI_ref <= 0)) stop("Reference cells must have positive total counts.", call. = FALSE)
  n_genes <- Matrix::colSums(counts > 0)
  eligible <- nUMI >= min_umi & n_genes >= min_genes & nUMI > 0
  if (!any(eligible)) stop("No bins pass min_umi and min_genes.", call. = FALSE)
  common <- intersect(rownames(counts), rownames(reference_counts))
  if (length(common) < 10L) {
    stop("RCTD requires at least 10 shared genes; check reference gene identifiers.", call. = FALSE)
  }
  backend <- .hd_rctd_backend()
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  selected <- barcodes[eligible]
  reference <- backend$Reference(reference_counts[common, , drop = FALSE],
                                  cell_types, nUMI_ref)
  # Record the actual reference used to learn profiles before create.RCTD()
  # replaces its stored reference with a small plotting-only subsample.
  reference_diagnostics <- list(
    input_cells = ncol(reference_counts),
    retained_cells = ncol(reference@counts),
    removed_cells = ncol(reference_counts) - ncol(reference@counts),
    input_cells_per_type = table(cell_types),
    retained_cells_per_type = table(reference@cell_types),
    retained_umi_range = range(reference@nUMI),
    stage = "Reference constructor, before create.RCTD reference subsampling"
  )
  spatial <- backend$SpatialRNA(
    coords = as.data.frame(coordinates[selected, , drop = FALSE]),
    counts = counts[common, selected, drop = FALSE], nUMI = nUMI[selected]
  )
  if (!"UMI_min" %in% names(create_args)) create_args$UMI_min <- min_umi
  model <- do.call(backend$create.RCTD,
                   c(list(spatialRNA = spatial, reference = reference,
                          max_cores = max_cores), create_args))
  model@config$RCTDmode <- doublet_mode
  for (key in names(run_args)) model@config[[key]] <- run_args[[key]]
  model <- backend$fitBulk(model)
  # Upstream can restrict the requested types or use supplied type profiles.
  # Its fitted profiles, rather than the original factor levels, define output.
  types <- as.character(model@cell_type_info$renorm[[2L]])
  if (length(types) < 2L || anyNA(types) || any(!nzchar(types)) ||
      anyDuplicated(types)) {
    stop("spacexr fitted profiles require at least two uniquely named cell types.",
         call. = FALSE)
  }
  model <- backend$choose_sigma_c(model)
  fitted_barcodes <- colnames(model@spatialRNA@counts)
  if (!length(fitted_barcodes) || anyDuplicated(fitted_barcodes) ||
      !all(fitted_barcodes %in% selected)) {
    stop("spacexr returned an empty or invalid set of fitted barcodes.", call. = FALSE)
  }
  if (is.null(batch_size)) batch_size <- length(fitted_barcodes)
  batches <- split(fitted_barcodes,
                   ceiling(seq_along(fitted_barcodes) / batch_size))
  out_weights <- vector("list", length(batches))
  out_annotations <- vector("list", length(batches))
  out_classification <- if (doublet_mode != "full") vector("list", length(batches)) else NULL
  raw_results <- if (keep_fit) vector("list", length(batches)) else NULL
  for (i in seq_along(batches)) {
    ids <- batches[[i]]
    chunk <- model
    chunk@spatialRNA <- .hd_rctd_subset_spatial(chunk@spatialRNA, ids)
    chunk@originalSpatialRNA <- .hd_rctd_subset_spatial(chunk@originalSpatialRNA, ids)
    chunk <- backend$fitPixels(chunk, doublet_mode = doublet_mode)
    parsed <- .hd_rctd_parse(chunk@results, ids, types, doublet_mode)
    out_weights[[i]] <- parsed$weights
    out_annotations[[i]] <- parsed$annotations
    if (doublet_mode != "full") {
      out_classification[[i]] <- .hd_rctd_native_classification(
        chunk@results, ids, doublet_mode
      )
    }
    if (keep_fit) raw_results[[i]] <- chunk@results
  }
  fitted_weights <- do.call(rbind, out_weights)
  fitted_annotations <- do.call(rbind, out_annotations)
  # Sparse selector expands results to every original bin without dense padding.
  selector <- Matrix::sparseMatrix(
    i = match(rownames(fitted_weights), barcodes),
    j = seq_len(nrow(fitted_weights)), x = 1,
    dims = c(length(barcodes), nrow(fitted_weights))
  )
  weights <- methods::as(selector %*% fitted_weights, "dgCMatrix")
  dimnames(weights) <- list(barcodes, types)
  annotations <- data.frame(
    barcode = barcodes, status = ifelse(eligible, "not_fit", "filtered"),
    filter_reason = ifelse(!eligible,
                           ifelse(nUMI <= 0, "zero_counts",
                                  ifelse(nUMI < min_umi, "low_umi", "low_genes")),
                           "upstream_filter"),
    nUMI = as.numeric(nUMI), n_genes = as.integer(n_genes),
    first_type = NA_character_, second_type = NA_character_,
    spot_class = NA_character_, dominant_type = NA_character_,
    max_weight = 0, stringsAsFactors = FALSE, row.names = barcodes
  )
  pos <- match(fitted_annotations$barcode, barcodes)
  for (key in setdiff(names(fitted_annotations), "barcode")) {
    annotations[pos, key] <- fitted_annotations[[key]]
  }
  annotations$filter_reason[pos] <- NA_character_
  annotations$filter_reason[annotations$status == "not_fit" &
    is.na(annotations$filter_reason)] <- "zero_fitted_weights"
  if (keep_fit) {
    model@results <- .hd_rctd_combine_results(raw_results, doublet_mode)
    model@internal_vars$cell_types_assigned <- TRUE
  }
  structure(list(
    weights = weights, annotations = annotations,
    classification = .hd_rctd_combine_classification(out_classification, doublet_mode),
    parameters = list(doublet_mode = doublet_mode, max_cores = max_cores,
                      min_umi = min_umi, min_genes = min_genes, seed = seed,
                      batch_size = batch_size, create_args = create_args,
                      run_args = run_args, backend = "dmcable/spacexr",
                      weight_units = "relative RNA/UMI contribution; not cell-number fraction",
                      upstream_config = model@config,
                      weight_policy = switch(doublet_mode,
                        full = "normalized_unrestricted_estimated_RNA_UMI_contributions",
                        doublet = "derived_annotation_weights_singlet_one_hot_reject_zero_native_pair_in_classification",
                        multi = "normalized_multi_mode_sub_weights_uncertainty_audited")),
    diagnostics = list(input_bins = ncol(counts), eligible_bins = sum(eligible),
                       fitted_bins = length(fitted_barcodes), shared_genes = length(common),
                       spatial_genes = nrow(counts), reference_genes = nrow(reference_counts),
                       reference = reference_diagnostics,
                       batches = length(batches), status_counts = table(annotations$status),
                       spacexr_version = backend$version),
    provenance = list(bin_size = hd$bin_size, source_bin_size = hd$source_bin_size,
                      input = hd$provenance),
    fit = if (keep_fit) model else NULL
  ), class = "SpatialESSVisiumHDRCTD")
}

.hd_rctd_native_classification <- function(results, ids, mode) {
  if (mode == "doublet") {
    return(list(
      results_df = results$results_df[ids, , drop = FALSE],
      weights_doublet = results$weights_doublet[ids, , drop = FALSE]
    ))
  }
  # Multi-mode results may also contain an unrestricted vector for every type.
  # Retain the selected mixture, scores and convergence without that large copy.
  if (!is.null(names(results))) results <- results[ids]
  pixel_results <- lapply(results, function(x) x[setdiff(names(x), "all_weights")])
  names(pixel_results) <- ids
  list(pixel_results = pixel_results)
}

.hd_rctd_combine_classification <- function(results, mode) {
  if (mode == "full") return(NULL)
  if (mode == "doublet") {
    return(list(
      results_df = do.call(rbind, lapply(results, `[[`, "results_df")),
      weights_doublet = do.call(rbind, lapply(results, `[[`, "weights_doublet"))
    ))
  }
  list(pixel_results = unlist(lapply(results, `[[`, "pixel_results"), recursive = FALSE))
}

.hd_rctd_number <- function(x, name, minimum, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum ||
      (integer && (x != floor(x) || x > .Machine$integer.max))) {
    stop(name, " must be a finite ", if (integer) "integer " else "number ",
         ">= ", minimum, ".", call. = FALSE)
  }
  if (integer) as.integer(x) else x
}

.hd_rctd_args <- function(x, name) {
  if (!is.list(x) || (length(x) &&
      (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
       anyDuplicated(names(x))))) {
    stop(name, " must be a list with unique, nonempty argument names.", call. = FALSE)
  }
}

.hd_rctd_counts <- function(x, name) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop(name, " must be a numeric gene-by-barcode matrix.", call. = FALSE)
  }
  x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
  if (any(dim(x) == 0L) || is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(name, " needs unique, nonempty gene and barcode names.", call. = FALSE)
  }
  if (any(!is.finite(x@x)) || any(x@x < 0) ||
      any(abs(x@x - round(x@x)) > 1e-8)) {
    stop(name, " must contain finite nonnegative raw integer counts.", call. = FALSE)
  }
  Matrix::drop0(x)
}

.hd_rctd_backend <- function() {
  required <- c("Reference", "SpatialRNA", "create.RCTD", "fitBulk",
                "choose_sigma_c", "fitPixels")
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    stop("Install the original dmcable/spacexr package to run Visium HD RCTD.", call. = FALSE)
  }
  if (!all(required %in% getNamespaceExports("spacexr"))) {
    stop("Visium HD RCTD requires the original dmcable/spacexr API (Reference, SpatialRNA, create.RCTD, fitBulk, choose_sigma_c, fitPixels). The Bioconductor spacexr API is not supported by this interface.", call. = FALSE)
  }
  functions <- lapply(required, function(name) getExportedValue("spacexr", name))
  names(functions) <- required
  c(functions, list(version = as.character(utils::packageVersion("spacexr"))))
}

.hd_rctd_subset_spatial <- function(spatial, ids) {
  spatial@counts <- spatial@counts[, ids, drop = FALSE]
  spatial@coords <- spatial@coords[ids, , drop = FALSE]
  spatial@nUMI <- spatial@nUMI[ids]
  spatial
}

.hd_rctd_normalize <- function(weights, ids, types) {
  weights <- methods::as(Matrix::Matrix(as.matrix(weights), sparse = TRUE), "dgCMatrix")
  if (is.null(rownames(weights)) || is.null(colnames(weights)) ||
      anyDuplicated(rownames(weights)) || anyDuplicated(colnames(weights)) ||
      !setequal(rownames(weights), ids) || !setequal(colnames(weights), types) ||
      any(!is.finite(weights@x)) || any(weights@x < 0)) {
    stop("spacexr returned invalid weights or mismatched barcode/type names.", call. = FALSE)
  }
  weights <- weights[ids, types, drop = FALSE]
  totals <- Matrix::rowSums(weights)
  normalized <- Matrix::Diagonal(x = ifelse(totals > 0, 1 / totals, 0)) %*% weights
  dimnames(normalized) <- list(ids, types)
  methods::as(Matrix::drop0(normalized), "dgCMatrix")
}

.hd_rctd_parse <- function(results, ids, types, mode) {
  annotation <- data.frame(
    barcode = ids, status = "ok", first_type = NA_character_,
    second_type = NA_character_, spot_class = NA_character_,
    dominant_type = NA_character_, max_weight = 0,
    stringsAsFactors = FALSE, row.names = ids
  )
  if (mode == "full") {
    weights <- .hd_rctd_normalize(results$weights, ids, types)
  } else if (mode == "doublet") {
    res <- as.data.frame(results$results_df)
    if (!all(c("spot_class", "first_type", "second_type") %in% names(res)) ||
        anyDuplicated(rownames(res)) || !setequal(rownames(res), ids)) {
      stop("Invalid doublet results_df returned by spacexr.", call. = FALSE)
    }
    res <- res[ids, , drop = FALSE]
    annotation$spot_class <- as.character(res$spot_class)
    if (anyNA(annotation$spot_class) || !all(annotation$spot_class %in%
        c("reject", "singlet", "doublet_certain", "doublet_uncertain"))) {
      stop("Invalid doublet spot_class returned by spacexr.", call. = FALSE)
    }
    annotation$first_type <- as.character(res$first_type)
    annotation$second_type <- as.character(res$second_type)
    rejected <- annotation$spot_class == "reject"
    certain <- annotation$spot_class == "doublet_certain"
    annotation$status[rejected] <- "reject"
    annotation$status[annotation$spot_class == "doublet_uncertain"] <- "ambiguous"
    dw <- as.matrix(results$weights_doublet)
    if (is.null(rownames(dw)) || anyDuplicated(rownames(dw)) ||
        !setequal(rownames(dw), ids) ||
        !all(c("first_type", "second_type") %in% colnames(dw))) {
      stop("Invalid weights_doublet returned by spacexr.", call. = FALSE)
    }
    dw <- dw[ids, c("first_type", "second_type"), drop = FALSE]
    if (any(!is.finite(dw[!rejected, , drop = FALSE])) || any(dw[!rejected, , drop = FALSE] < 0) ||
        any(!annotation$first_type[!rejected] %in% types) ||
        any(!annotation$second_type[!rejected & annotation$spot_class != "singlet"] %in% types)) {
      stop("Invalid confident doublet labels or weights returned by spacexr.", call. = FALSE)
    }
    singlet <- annotation$spot_class == "singlet"
    if (any(singlet)) dw[singlet, ] <- cbind(rep(1, sum(singlet)), rep(0, sum(singlet)))
    dw[rejected, ] <- 0
    ij <- which(dw > 0, arr.ind = TRUE)
    labels <- cbind(annotation$first_type, annotation$second_type)[ij]
    weights <- Matrix::sparseMatrix(
      i = ij[, 1L], j = match(labels, types), x = dw[ij],
      dims = c(length(ids), length(types)), dimnames = list(ids, types)
    )
    weights <- .hd_rctd_normalize(weights, ids, types)
    annotation$first_type[rejected] <- NA_character_
    annotation$second_type[rejected | annotation$spot_class == "singlet"] <- NA_character_
  } else {
    if (length(results) != length(ids)) stop("Invalid multi-mode results returned by spacexr.", call. = FALSE)
    if (!is.null(names(results))) {
      if (anyDuplicated(names(results)) || !setequal(names(results), ids)) {
        stop("Multi-mode result barcodes do not match fitted bins.", call. = FALSE)
      }
      results <- results[ids]
    }
    entries <- lapply(seq_along(ids), function(i) {
      result <- results[[i]]
      labels <- as.character(result$cell_type_list)
      w <- as.numeric(result$sub_weights)
      confidence <- as.logical(result$conf_list)
      if (!length(labels) || length(w) != length(labels) ||
          length(confidence) != length(labels) || anyNA(confidence) ||
          anyDuplicated(labels) || any(!labels %in% types) ||
          any(!is.finite(w)) || any(w < 0)) {
        stop("Invalid multi-mode labels, weights or confidence returned by spacexr.", call. = FALSE)
      }
      if (any(!confidence)) annotation$status[i] <<- "ambiguous"
      list(i = rep(i, length(labels)), j = match(labels, types), x = w)
    })
    weights <- Matrix::sparseMatrix(
      i = unlist(lapply(entries, `[[`, "i")),
      j = unlist(lapply(entries, `[[`, "j")),
      x = unlist(lapply(entries, `[[`, "x")),
      dims = c(length(ids), length(types)), dimnames = list(ids, types)
    )
    weights <- .hd_rctd_normalize(weights, ids, types)
  }
  for (i in seq_along(ids)) {
    w <- as.numeric(weights[i, ])
    maximum <- max(w)
    annotation$max_weight[i] <- maximum
    if (maximum <= 0 && annotation$status[i] == "ok") annotation$status[i] <- "not_fit"
    tied <- sum(abs(w - maximum) <= 1e-10) > 1L
    if (maximum > 0 && tied && annotation$status[i] == "ok") annotation$status[i] <- "ambiguous"
    if (maximum > 0 && !tied && annotation$status[i] == "ok") {
      annotation$dominant_type[i] <- types[which.max(w)]
    }
  }
  list(weights = weights, annotations = annotation)
}

.hd_rctd_combine_results <- function(results, mode) {
  if (mode == "multi") return(unlist(results, recursive = FALSE))
  fields <- names(results[[1L]])
  combined <- lapply(fields, function(key) {
    values <- lapply(results, `[[`, key)
    if (key %in% c("weights", "weights_doublet", "results_df")) do.call(rbind, values)
    else unlist(values, recursive = FALSE)
  })
  names(combined) <- fields
  combined
}
