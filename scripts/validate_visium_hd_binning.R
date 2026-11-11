#!/usr/bin/env Rscript
# Compare against genuinely existing native target outputs, never another
# automatic fallback. This validates binning, not RCTD biological accuracy.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1L && args %in% c("--help", "-h")) {
  cat("Usage: Rscript scripts/validate_visium_hd_binning.R OUTS SOURCE_UM TARGETS OUTPUT_TSV [filtered|raw]\n",
      "Example: ... /data/sample/outs 8 16,32 binning_check.tsv filtered\n",
      "Each target square_NNNum directory must already exist (e.g. Space Ranger --custom-bin-size=32).\n",
      "The source needs raw counts, complete tissue positions, and filtered feature metadata for filtered policy.\n",
      sep = "")
  quit(status = 0L)
}
if (!length(args) %in% c(4L, 5L)) stop("Expected 4 or 5 arguments; use --help.", call. = FALSE)
if (!requireNamespace("SpatialESS", quietly = TRUE) ||
    !"feature_policy" %in% names(formals(SpatialESS::read_visium_hd))) {
  stop("Install SpatialESS 0.1.4.9002 or newer with the explicit feature policy API.", call. = FALSE)
}
root <- normalizePath(args[1L], mustWork = TRUE)
source_size <- suppressWarnings(as.numeric(args[2L]))
targets <- suppressWarnings(as.numeric(strsplit(args[3L], ",", fixed = TRUE)[[1L]]))
policy <- if (length(args) == 5L) match.arg(args[5L], c("filtered", "raw")) else "filtered"
sizes <- c(source_size, targets)
if (anyNA(sizes) || any(!is.finite(sizes)) || any(sizes < 2 | sizes %% 2 != 0) ||
    !length(targets) || any(targets <= source_size | targets %% source_size != 0) ||
    anyDuplicated(targets)) stop("Sizes require distinct even target multiples larger than SOURCE_UM.", call. = FALSE)
output <- path.expand(args[4L])
if (file.exists(output)) stop("Output already exists; choose a new TSV path.", call. = FALSE)
roots <- unique(c(root, file.path(root, "binned_outputs"),
                  file.path(root, "outs", "binned_outputs")))
square <- function(size) {
  name <- sprintf("square_%03dum", size)
  paths <- unique(c(if (basename(root) == name) root,
                    file.path(roots, name)))
  paths <- paths[dir.exists(paths)]
  if (length(paths) != 1L) stop("Require exactly one existing native ", name,
                               " directory; no automatic substitute is permitted.", call. = FALSE)
  paths[[1L]]
}
source_dir <- square(source_size)
rows <- list()
for (target in targets) {
  row <- data.frame(source_um = source_size, target_um = target,
    feature_policy = policy, gene_ids_equal = FALSE, count_barcodes_equal = FALSE,
    counts_exact = FALSE, maximum_absolute_count_difference = NA_real_,
    tissue_barcodes_equal = FALSE, tissue_mask_equal = FALSE,
    native_bins = NA_integer_, rebinned_bins = NA_integer_,
    status = "failed", error = "", stringsAsFactors = FALSE)
  problem <- tryCatch({
    native_dir <- square(target)
    native <- SpatialESS::read_visium_hd(native_dir, bin_size = target,
      gene_column = "id", feature_policy = policy, in_tissue = policy == "filtered")
    rebuilt <- SpatialESS::read_visium_hd(source_dir, bin_size = target,
      source_bin_size = source_size, gene_column = "id", feature_policy = policy,
      in_tissue = policy == "filtered", allow_filtered_aggregation = FALSE)
    row$native_bins <- ncol(native$counts)
    row$rebinned_bins <- ncol(rebuilt$counts)
    row$gene_ids_equal <- setequal(rownames(native$counts), rownames(rebuilt$counts))
    row$count_barcodes_equal <- setequal(colnames(native$counts), colnames(rebuilt$counts))
    if (row$gene_ids_equal && row$count_barcodes_equal) {
      delta <- methods::as(native$counts - rebuilt$counts[
        rownames(native$counts), colnames(native$counts), drop = FALSE], "dgCMatrix")
      row$maximum_absolute_count_difference <- if (length(delta@x)) max(abs(delta@x)) else 0
      row$counts_exact <- row$maximum_absolute_count_difference == 0
    }
    a <- native$all_positions
    b <- rebuilt$all_positions
    if (!isTRUE(native$provenance$tissue_mask_available) ||
        !isTRUE(rebuilt$provenance$tissue_mask_available)) {
      stop("A supplied tissue mask is required; an assumed all-tissue mask cannot establish native equivalence.")
    }
    if (is.null(a) || is.null(b) || !"in_tissue" %in% names(a) ||
        !"in_tissue" %in% names(b)) stop("Both objects need complete positions for tissue-mask validation.")
    row$tissue_barcodes_equal <- setequal(a$barcode, b$barcode)
    if (row$tissue_barcodes_equal) {
      row$tissue_mask_equal <- identical(as.integer(a$in_tissue),
        as.integer(b$in_tissue[match(a$barcode, b$barcode)]))
    }
    if (!all(row$gene_ids_equal, row$count_barcodes_equal, row$counts_exact,
             row$tissue_barcodes_equal, row$tissue_mask_equal)) {
      stop("Native count/feature/barcode/tissue-mask equality was not established.")
    }
    NULL
  }, error = function(e) e)
  if (is.null(problem)) row$status <- "pass" else row$error <- conditionMessage(problem)
  rows[[as.character(target)]] <- row
}
result <- do.call(rbind, rows)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.table(result, output, sep = "\t", quote = TRUE, row.names = FALSE)
print(result)
if (any(result$status != "pass")) stop("Native equivalence is incomplete; inspect ", output, ".", call. = FALSE)
message("Exact count and tissue-mask checks passed for the supplied native outputs; this does not validate cell-type accuracy.")
