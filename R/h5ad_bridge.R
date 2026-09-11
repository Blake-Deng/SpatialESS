#' Export an H5AD expression matrix for SpatialESS
#'
#' The Python bridge computes per-cell library sizes over all genes before
#' selecting the requested signaling genes. This avoids changing CP10K
#' denominators when a gene list is supplied.
#'
#' @param h5ad Path to a SPARKLE-corrected H5AD file.
#' @param group_col Column containing shared cell-group labels.
#' @param coordinate_cols Two or three spatial-coordinate columns.
#' @param gene_list Genes to export. When NULL, all genes are exported.
#' @param output_dir Destination bundle directory.
#' @param normalization Either code{"log1p_cp10k"} or code{"none"}.
#' @param metadata_file Optional TSV containing cell labels and coordinates.
#' @param cell_col Cell identifier column in code{metadata_file}.
#' @param layer H5AD layer to export; NULL uses code{X}.
#' @param python Python executable with anndata, numpy, pandas and scipy.
#' @param bridge_script Optional explicit path to the Python bridge.
#' @param overwrite Replace an existing bundle.
#' @return Normalized bundle path, invisibly.
#' @export
export_spatialess_h5ad_bundle <- function(
    h5ad, group_col, coordinate_cols,
    gene_list = NULL, output_dir,
    normalization = c("log1p_cp10k", "none"),
    metadata_file = NULL, cell_col = "cell", layer = NULL,
    python = Sys.which("python3"), bridge_script = NULL,
    overwrite = FALSE) {
  normalization <- match.arg(normalization)
  if (!file.exists(h5ad)) stop("H5AD file does not exist.", call. = FALSE)
  if (length(coordinate_cols) < 2L || length(coordinate_cols) > 3L) {
    stop("coordinate_cols must contain two or three columns.", call. = FALSE)
  }
  if (!nzchar(python)) stop("Cannot find a Python executable.", call. = FALSE)
  if (is.null(bridge_script)) {
    bridge_script <- system.file(
      "python", "h5ad_to_spatialess.py", package = "SpatialESSSPARKLE"
    )
  }
  if (!nzchar(bridge_script) || !file.exists(bridge_script)) {
    stop("Cannot locate h5ad_to_spatialess.py.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  args <- c(
    bridge_script,
    normalizePath(h5ad, mustWork = TRUE),
    normalizePath(output_dir, mustWork = TRUE),
    "--group-col", group_col,
    "--coordinate-cols", paste(coordinate_cols, collapse = ","),
    "--normalization", normalization
  )
  temporary_gene_file <- NULL
  if (!is.null(gene_list)) {
    gene_list <- unique(as.character(gene_list))
    gene_list <- gene_list[!is.na(gene_list) & nzchar(gene_list)]
    if (!length(gene_list)) stop("gene_list is empty.", call. = FALSE)
    temporary_gene_file <- tempfile(fileext = ".txt")
    writeLines(gene_list, temporary_gene_file)
    args <- c(args, "--gene-list", temporary_gene_file)
  }
  on.exit({
    if (!is.null(temporary_gene_file)) unlink(temporary_gene_file)
  }, add = TRUE)
  if (!is.null(metadata_file)) {
    if (!file.exists(metadata_file)) {
      stop("metadata_file does not exist.", call. = FALSE)
    }
    args <- c(
      args, "--metadata-tsv", normalizePath(metadata_file, mustWork = TRUE),
      "--cell-col", cell_col
    )
  }
  if (!is.null(layer) && nzchar(layer)) args <- c(args, "--layer", layer)
  if (isTRUE(overwrite)) args <- c(args, "--overwrite")

  output <- system2(
    python, args = shQuote(args), stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      paste(c("H5AD bridge failed:", output), collapse = "\n"),
      call. = FALSE
    )
  }
  path <- normalizePath(output_dir, mustWork = TRUE)
  attr(path, "bridge_output") <- output
  invisible(path)
}

#' Read a SpatialESS bridge bundle
#'
#' @param path Bundle directory produced by
#'   code{export_spatialess_h5ad_bundle()}.
#' @return A list containing expression, coordinates, group and metadata.
#' @export
read_spatialess_bundle <- function(path) {
  required <- c("expression.mtx.gz", "genes.tsv", "cells.tsv", "metadata.tsv")
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) {
    stop("Bundle is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  genes <- readLines(file.path(path, "genes.tsv"))
  cells <- readLines(file.path(path, "cells.tsv"))
  expression <- methods::as(
    Matrix::readMM(gzfile(file.path(path, "expression.mtx.gz"))),
    "CsparseMatrix"
  )
  if (!identical(dim(expression), c(length(genes), length(cells)))) {
    stop("Bundle matrix dimensions do not match genes/cells.", call. = FALSE)
  }
  rownames(expression) <- genes
  colnames(expression) <- cells

  metadata <- utils::read.delim(
    file.path(path, "metadata.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!identical(as.character(metadata$cell), cells)) {
    stop("Bundle metadata and expression cell order differ.", call. = FALSE)
  }
  coordinate_columns <- grep("^coordinate_[123]$", names(metadata), value = TRUE)
  if (length(coordinate_columns) < 2L) {
    stop("Bundle does not contain two spatial coordinates.", call. = FALSE)
  }
  coordinates <- as.matrix(metadata[, coordinate_columns, drop = FALSE])
  storage.mode(coordinates) <- "double"
  rownames(coordinates) <- cells
  group <- droplevels(factor(metadata$group))
  names(group) <- cells
  list(
    expression = expression,
    coordinates = coordinates,
    group = group,
    metadata = metadata,
    manifest = file.path(path, "manifest.json")
  )
}

#' Run SpatialESS directly from a SPARKLE H5AD file
#'
#' @param h5ad,group_col,coordinate_cols,normalization,metadata_file,cell_col,
#'   layer,python,bridge_script Arguments passed to
#'   code{export_spatialess_h5ad_bundle()}.
#' @param lr,complex,cofactor CellChat-compatible database tables.
#' @param gene_list Genes to export. If NULL, they are obtained using
#'   code{CellChat::extractGene()}.
#' @param output_dir Optional persistent bridge directory.
#' @param keep_bundle Keep a temporary bridge directory.
#' @param filter_unresolved Legacy option: TRUE filters incomplete LR with a
#'   warning and a full audit; FALSE uses strict checking.
#' @param lr_missing Explicitly select `"filter"` or `"error"`, as in SpatialESS.
#'   NULL preserves `filter_unresolved`. Conflicting explicit options are rejected.
#' @param ... Arguments passed to
#'   code{SpatialESS::spatialess()}.
#' @return A code{SpatialESSResult} with bridge diagnostics.
#' @export
spatialess_from_h5ad <- function(
    h5ad, group_col, coordinate_cols, lr, complex, cofactor,
    gene_list = NULL,
    normalization = c("log1p_cp10k", "none"),
    metadata_file = NULL, cell_col = "cell", layer = NULL,
    python = Sys.which("python3"), bridge_script = NULL,
    output_dir = NULL, keep_bundle = FALSE,
    filter_unresolved = TRUE, lr_missing = NULL, ...) {
  normalization <- match.arg(normalization)
  if (!is.logical(filter_unresolved) || length(filter_unresolved) != 1L ||
      is.na(filter_unresolved)) {
    stop("filter_unresolved must be TRUE or FALSE.", call. = FALSE)
  }
  legacy_mode <- if (filter_unresolved) "filter" else "error"
  if (is.null(lr_missing)) {
    lr_missing <- legacy_mode
  } else {
    lr_missing <- match.arg(lr_missing, c("error", "filter"))
    if (!missing(filter_unresolved) && lr_missing != legacy_mode) {
      stop("filter_unresolved and lr_missing specify conflicting policies.", call. = FALSE)
    }
  }
  if (is.null(gene_list)) {
    if (!requireNamespace("CellChat", quietly = TRUE)) {
      stop("CellChat is required when gene_list is NULL.", call. = FALSE)
    }
    gene_list <- CellChat::extractGene(list(
      interaction = lr, complex = complex, cofactor = cofactor
    ))
  }
  temporary_bundle <- is.null(output_dir)
  if (temporary_bundle) output_dir <- tempfile("spatialess_h5ad_")
  if (temporary_bundle && !isTRUE(keep_bundle)) {
    on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  export_spatialess_h5ad_bundle(
    h5ad = h5ad,
    group_col = group_col,
    coordinate_cols = coordinate_cols,
    gene_list = gene_list,
    output_dir = output_dir,
    normalization = normalization,
    metadata_file = metadata_file,
    cell_col = cell_col,
    layer = layer,
    python = python,
    bridge_script = bridge_script,
    overwrite = TRUE
  )
  bundle <- read_spatialess_bundle(output_dir)

  lr <- as.data.frame(lr, stringsAsFactors = FALSE)
  result <- SpatialESS::spatialess(
    expression = bundle$expression,
    coordinates = bundle$coordinates,
    group = bundle$group,
    lr = lr,
    complex = complex,
    cofactor = cofactor,
    lr_missing = lr_missing,
    ...
  )
  coverage <- result$lr_filter
  excluded_lr <- if (is.null(coverage)) 0L else nrow(coverage$excluded)
  if (!is.null(coverage)) {
    utils::write.table(coverage$audit, file.path(output_dir, "lr_audit.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(coverage$missing_genes, file.path(output_dir, "lr_missing_genes.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  result$bridge <- list(
    source_h5ad = normalizePath(h5ad, mustWork = TRUE),
    normalization = normalization,
    exported_genes = nrow(bundle$expression),
    cells = ncol(bundle$expression),
    groups = nlevels(bundle$group),
    computable_lr = nrow(lr) - excluded_lr,
    excluded_lr = excluded_lr,
    lr_missing = lr_missing,
    bundle = if (!temporary_bundle || isTRUE(keep_bundle)) {
      normalizePath(output_dir, mustWork = TRUE)
    } else {
      NULL
    }
  )
  result
}
