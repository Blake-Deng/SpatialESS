#' Prepare dominant-type Visium HD bins for SpatialESS
#'
#' Retains accepted RCTD bins whose largest type proportion meets
#' \code{min_weight}. These are bin-level dominant-type labels, not reconstructed
#' single cells or cell-type-specific expression. The complete RCTD proportions
#' remain available in the deconvolution result.
#'
#' @param hd A \code{SpatialESSVisiumHD} object.
#' @param deconvolution A result from \code{run_visium_hd_rctd()}.
#' @param min_weight Minimum dominant-type proportion, between zero and one.
#' @param normalization \code{"log1p_cp10k"} normalizes each bin by its
#'   all-gene UMI total before taking log1p; \code{"none"} retains raw counts.
#' @return A list with aligned sparse \code{expression}, micrometre
#'   \code{coordinates}, \code{group}, retained \code{metadata}, and an
#'   all-bin inclusion \code{audit}.
#' @export
prepare_visium_hd_spatialess <- function(
    hd, deconvolution, min_weight = 0.5,
    normalization = c("log1p_cp10k", "none")) {
  normalization <- match.arg(normalization)
  if (!inherits(hd, "SpatialESSVisiumHD")) {
    stop("hd must be a SpatialESSVisiumHD object.", call. = FALSE)
  }
  if (length(min_weight) != 1L || !is.numeric(min_weight) ||
      !is.finite(min_weight) || min_weight < 0 || min_weight > 1) {
    stop("min_weight must be one finite number between zero and one.",
         call. = FALSE)
  }
  ids <- colnames(hd$counts)
  if (is.null(ids) || anyDuplicated(ids) ||
      !identical(rownames(hd$coordinates), ids) ||
      !identical(rownames(hd$metadata), ids)) {
    stop("hd counts, coordinates and metadata must have identical bin order.",
         call. = FALSE)
  }
  if (length(unique(.visium_hd_suffix(ids))) > 1L) {
    stop("Run bin-level SpatialESS separately for each barcode suffix/sample to avoid cross-sample spatial edges.",
         call. = FALSE)
  }
  if (!is.null(deconvolution$provenance$bin_size) &&
      !identical(as.numeric(deconvolution$provenance$bin_size),
                 as.numeric(hd$bin_size))) {
    stop("hd and deconvolution have different bin sizes.", call. = FALSE)
  }
  annotations <- deconvolution$annotations
  weights <- deconvolution$weights
  required <- c("barcode", "status", "dominant_type", "max_weight")
  if (!is.data.frame(annotations) ||
      !all(required %in% names(annotations))) {
    stop("deconvolution annotations require barcode, status, dominant_type and max_weight.",
         call. = FALSE)
  }
  if (anyNA(annotations$barcode) || anyDuplicated(annotations$barcode) ||
      !setequal(as.character(annotations$barcode), ids)) {
    stop("deconvolution annotations must cover every HD bin exactly once.",
         call. = FALSE)
  }
  if (is.null(dim(weights)) || is.null(rownames(weights)) ||
      anyDuplicated(rownames(weights)) ||
      !setequal(rownames(weights), ids) ||
      is.null(colnames(weights)) || anyDuplicated(colnames(weights))) {
    stop("deconvolution weights require unique bin rows and cell-type columns aligned to hd.",
         call. = FALSE)
  }
  annotations <- annotations[match(ids, annotations$barcode), , drop = FALSE]
  if ("nUMI" %in% names(annotations) &&
      (!is.numeric(annotations$nUMI) || anyNA(annotations$nUMI) ||
       any(abs(annotations$nUMI - Matrix::colSums(hd$counts)) > 1e-6))) {
    stop("hd all-gene UMI totals differ from the deconvolution input.", call. = FALSE)
  }
  weights <- methods::as(Matrix::Matrix(weights[ids, , drop = FALSE], sparse = TRUE),
                         "dgCMatrix")
  if (anyNA(weights@x) || any(!is.finite(weights@x)) || any(weights@x < 0)) {
    stop("deconvolution weights must be finite and non-negative.", call. = FALSE)
  }
  totals <- Matrix::rowSums(weights)
  nonzero <- totals > 0
  if (any(abs(totals[nonzero] - 1) > 1e-6)) {
    stop("nonzero deconvolution weight rows must sum to one.", call. = FALSE)
  }
  accepted <- !is.na(annotations$status) & annotations$status == "ok"
  labels <- as.character(annotations$dominant_type)
  labelled <- !is.na(labels) & nzchar(labels) & labels %in% colnames(weights)
  if (!is.numeric(annotations$max_weight)) {
    stop("annotation max_weight must be numeric.", call. = FALSE)
  }
  valid_weight <- is.finite(annotations$max_weight) & annotations$max_weight >= 0 &
    annotations$max_weight <= 1
  # Check the claimed dominant label against the actual proportions. An old or
  # edited annotation file must not silently annotate a different mixture.
  check <- which(accepted & labelled & valid_weight & nonzero)
  if (length(check)) {
    observed <- as.numeric(weights[cbind(check, match(labels[check], colnames(weights)))])
    # Column-wise sparse maxima avoid one R submatrix allocation per HD bin.
    all_max <- numeric(nrow(weights))
    for (j in seq_len(ncol(weights))) {
      if (weights@p[j] == weights@p[j + 1L]) next
      at <- seq.int(weights@p[j] + 1L, weights@p[j + 1L])
      rows <- weights@i[at] + 1L
      all_max[rows] <- pmax(all_max[rows], weights@x[at])
    }
    row_max <- all_max[check]
    if (any(abs(observed - annotations$max_weight[check]) > 1e-6) ||
        any(abs(observed - row_max) > 1e-6)) {
      stop("dominant_type or max_weight disagrees with deconvolution weights.",
           call. = FALSE)
    }
  }
  keep <- accepted & labelled & valid_weight & nonzero &
    annotations$max_weight >= min_weight
  reason <- rep("retained", length(ids))
  reason[!accepted] <- "rctd_not_accepted"
  reason[accepted & (!labelled | !valid_weight | !nonzero)] <- "no_confident_type"
  reason[accepted & labelled & valid_weight & nonzero &
           annotations$max_weight < min_weight] <- "below_min_weight"
  audit <- annotations
  rownames(audit) <- ids
  audit$included <- keep
  audit$exclusion_reason <- reason
  if (!any(keep)) {
    stop("No accepted dominant-type bins meet min_weight; inspect RCTD annotations and proportions.",
         call. = FALSE)
  }
  expression <- methods::as(hd$counts[, keep, drop = FALSE], "dgCMatrix")
  library_size <- Matrix::colSums(hd$counts)
  if (any(!is.finite(library_size[keep]) | library_size[keep] <= 0)) {
    stop("retained HD bins must have positive finite all-gene UMI totals.",
         call. = FALSE)
  }
  if (normalization == "log1p_cp10k") {
    expression <- expression %*% Matrix::Diagonal(x = 1e4 / library_size[keep])
    expression <- methods::as(expression, "dgCMatrix")
    dimnames(expression) <- list(rownames(hd$counts), ids[keep])
    expression@x <- log1p(expression@x)
  }
  group <- factor(stats::setNames(labels[keep], ids[keep]))
  list(
    expression = expression,
    coordinates = hd$coordinates[keep, , drop = FALSE],
    group = group,
    metadata = hd$metadata[keep, , drop = FALSE],
    audit = audit,
    parameters = list(
      bin_size = hd$bin_size, source_bin_size = hd$source_bin_size,
      coordinate_units = "micrometres", min_weight = min_weight,
      normalization = normalization, normalization_genes = "all input genes",
      annotation = "RCTD dominant type per bin",
      expression_model = "observed bin expression; no cell-type unmixing"
    )
  )
}

#' Run SpatialESS on accepted dominant-type Visium HD bins
#'
#' This is a bin-level analysis using the observed expression of each bin.
#' RCTD proportions are retained for inspection and are not used to split a bin
#' into cells or infer unobserved cell-type-specific expression.
#'
#' @inheritParams prepare_visium_hd_spatialess
#' @param lr,complex,cofactor CellChat-compatible database component tables.
#' @param ... Other arguments to \code{spatialess()}. Coordinates are already
#'   in micrometres, so \code{ratio} must equal one.
#' @return A \code{SpatialESSResult}, with \code{visium_hd} containing the bin
#'   audit, metadata, proportions and preparation parameters.
#' @export
spatialess_from_visium_hd <- function(
    hd, deconvolution, lr, complex, cofactor, min_weight = 0.5,
    normalization = c("log1p_cp10k", "none"), ...) {
  prepared <- prepare_visium_hd_spatialess(
    hd, deconvolution, min_weight, normalization
  )
  dots <- list(...)
  if ("ratio" %in% names(dots)) {
    if (!identical(as.numeric(dots$ratio), 1)) {
      stop("Visium HD coordinates are in micrometres; ratio must equal 1.",
           call. = FALSE)
    }
    dots$ratio <- NULL
  }
  if (any(names(dots) %in% c("expression", "coordinates", "group"))) {
    stop("expression, coordinates and group are prepared from hd and RCTD.",
         call. = FALSE)
  }
  result <- do.call(spatialess, c(list(
    expression = prepared$expression, coordinates = prepared$coordinates,
    group = prepared$group, lr = lr, complex = complex, cofactor = cofactor,
    ratio = 1
  ), dots))
  result$visium_hd <- list(
    parameters = prepared$parameters, audit = prepared$audit,
    metadata = prepared$metadata,
    weights = deconvolution$weights,
    deconvolution_parameters = deconvolution$parameters,
    input_provenance = hd$provenance
  )
  result
}
