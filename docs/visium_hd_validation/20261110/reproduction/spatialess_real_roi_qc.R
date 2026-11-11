#!/usr/bin/env Rscript

# Read-only scientific QC of completed full-mode RCTD ROI results. This script
# writes summary tables and plots; it never changes counts, labels, fits or QC.
args <- commandArgs(trailingOnly = TRUE)
usage <- paste0(
  "Usage: Rscript spatialess_real_roi_qc.R RUN_DIR [--overwrite]\n",
  "       RUN_DIR=/absolute/run/path Rscript spatialess_real_roi_qc.R\n\n",
  "Requires completed RUN_DIR/real_roi_cli_results/bin8, bin16 and bin32\n",
  "with visium_hd_input.rds and rctd_result.rds; full RCTD mode only.\n",
  "Also requires reports/roi_retained_feature_ids.tsv (id,name,type) and\n",
  "inputs/ovarian_reference.rds (list(counts,cell_types)).\n",
  "Writes reports/real_roi_qc_*.tsv, reference/RNA contribution tables, and\n",
  "three headless PNG figures under reports/plots.\n",
  "These are internal QC and RNA-contribution summaries, not cell counts,\n",
  "calibrated confidence or independent biological accuracy validation.\n"
)
if (any(args %in% c("-h", "--help"))) {
  cat(usage)
  quit(status = 0L)
}
overwrite <- "--overwrite" %in% args
positional <- args[args != "--overwrite"]
if (length(positional) > 1L || any(startsWith(positional, "--"))) {
  stop(usage, call. = FALSE)
}
run_dir <- if (length(positional)) positional[[1L]] else Sys.getenv("RUN_DIR", unset = "")
if (!nzchar(run_dir)) stop(usage, call. = FALSE)
run_dir <- normalizePath(path.expand(run_dir), mustWork = TRUE)
if (!dir.exists(run_dir)) stop("RUN_DIR must be a directory.", call. = FALSE)
report_dir <- file.path(run_dir, "reports")
plot_dir <- file.path(report_dir, "plots")
feature_file <- file.path(report_dir, "roi_retained_feature_ids.tsv")
reference_file <- file.path(run_dir, "inputs", "ovarian_reference.rds")
bin_sizes <- c(8L, 16L, 32L)
bin_dirs <- file.path(run_dir, "real_roi_cli_results", paste0("bin", bin_sizes))
required <- c(feature_file, reference_file,
              file.path(bin_dirs, "visium_hd_input.rds"),
              file.path(bin_dirs, "rctd_result.rds"))
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing required completed-run files:\n", paste(missing, collapse = "\n"), call. = FALSE)
}
outputs <- c(
  file.path(report_dir, c("real_roi_qc_summary.tsv", "real_roi_qc_status_counts.tsv",
                         "real_roi_reference_type_counts.tsv",
                         "real_roi_umi_weighted_rna_contributions.tsv",
                         "real_roi_gene_identifier_audit.tsv", "real_roi_qc_notes.txt")),
  file.path(plot_dir, c("01_depth_genes_and_retention.png",
                       "02_umi_weighted_rna_contributions.png",
                       "03_dominant_rna_type_bins.png"))
)
if (!overwrite && any(file.exists(outputs))) {
  stop("QC outputs already exist. Use --overwrite to replace only these named outputs.", call. = FALSE)
}
if (!requireNamespace("Matrix", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("This QC script requires the existing Matrix and ggplot2 packages.", call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = TRUE,
                     row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
assert_ids <- function(x, name) {
  if (is.null(x) || anyNA(x) || any(!nzchar(as.character(x))) || anyDuplicated(x)) {
    stop(name, " requires unique, nonmissing identifiers.", call. = FALSE)
  }
}
quantile_values <- function(x, prefix) {
  probs <- c(0, .05, .25, .5, .75, .95, 1)
  values <- if (length(x)) stats::quantile(x, probs, names = FALSE, na.rm = FALSE) else rep(NA_real_, length(probs))
  stats::setNames(as.list(as.numeric(values)),
                  paste0(prefix, c("_min", "_p05", "_p25", "_median", "_p75", "_p95", "_max")))
}
sparse_row_max <- function(x) {
  answer <- numeric(nrow(x))
  for (j in seq_len(ncol(x))) {
    if (x@p[j] == x@p[j + 1L]) next
    at <- seq.int(x@p[j] + 1L, x@p[j + 1L])
    rows <- x@i[at] + 1L
    answer[rows] <- pmax(answer[rows], x@x[at])
  }
  answer
}
open_png <- function(path, width, height, res = 160) {
  if (isTRUE(capabilities("cairo"))) {
    grDevices::png(path, width = width, height = height, res = res, type = "cairo")
  } else if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = width, height = height, units = "px", res = res)
  } else {
    stop("Headless PNG output needs an R build with Cairo support or optional ragg.", call. = FALSE)
  }
}
save_plot <- function(plot, path, width, height) {
  open_png(path, width, height)
  tryCatch(print(plot), finally = grDevices::dev.off())
}

feature_info <- utils::read.delim(feature_file, stringsAsFactors = FALSE,
                                 check.names = FALSE, quote = "\"", comment.char = "")
if (!all(c("id", "name", "type") %in% names(feature_info))) {
  stop("Feature audit TSV must contain id, name and type columns.", call. = FALSE)
}
feature_info <- feature_info[feature_info$type == "Gene Expression", , drop = FALSE]
if (!nrow(feature_info) || anyNA(feature_info$id) || anyNA(feature_info$name)) {
  stop("Feature audit TSV has no usable Gene Expression feature mapping.", call. = FALSE)
}
reference <- readRDS(reference_file)
if (!is.list(reference) || !all(c("counts", "cell_types") %in% names(reference)) ||
    is.null(dim(reference$counts)) || is.null(colnames(reference$counts))) {
  stop("Reference RDS must contain gene-by-cell counts and named cell_types.", call. = FALSE)
}
assert_ids(colnames(reference$counts), "Reference count columns")
assert_ids(names(reference$cell_types), "Reference cell_types names")
if (!setequal(names(reference$cell_types), colnames(reference$counts))) {
  stop("Reference annotations do not match its count columns.", call. = FALSE)
}
reference_types <- reference$cell_types[colnames(reference$counts)]
if (anyNA(reference_types) || any(!nzchar(trimws(as.character(reference_types))))) {
  stop("Reference cell types contain missing or empty labels.", call. = FALSE)
}
input_reference_counts <- table(reference_types)
reference_input_cells <- ncol(reference$counts)
# Counts can be large; only names and the observed annotation census are needed.
rm(reference)
invisible(gc())

summary_rows <- metrics <- status_rows <- reference_rows <- contribution_rows <- gene_rows <- list()
for (i in seq_along(bin_sizes)) {
  size <- bin_sizes[[i]]
  message("[QC bin", size, "] Reading completed full-mode results")
  hd <- readRDS(file.path(bin_dirs[[i]], "visium_hd_input.rds"))
  fit <- readRDS(file.path(bin_dirs[[i]], "rctd_result.rds"))
  if (!inherits(hd, "SpatialESSVisiumHD") || is.null(dim(hd$counts)) ||
      is.null(hd$coordinates) || !is.data.frame(hd$metadata)) {
    stop("bin", size, " has an invalid Visium HD input object.", call. = FALSE)
  }
  if (!identical(as.numeric(hd$bin_size), as.numeric(size))) {
    stop("bin", size, " input records a different bin size.", call. = FALSE)
  }
  if (!identical(fit$parameters$doublet_mode, "full") ||
      grepl("derived|one_hot", paste(fit$parameters$weight_policy, collapse = ";"))) {
    stop("This RNA-contribution QC requires full-mode mixture weights. Doublet annotation weights cannot be treated as raw RNA fractions.", call. = FALSE)
  }
  ids <- colnames(hd$counts)
  genes <- rownames(hd$counts)
  assert_ids(ids, paste0("bin", size, " count barcodes"))
  assert_ids(genes, paste0("bin", size, " count genes"))
  if (!inherits(hd$counts, "Matrix")) hd$counts <- Matrix::Matrix(hd$counts, sparse = TRUE)
  if (any(!is.finite(hd$counts@x)) || any(hd$counts@x < 0) ||
      any(hd$counts@x != floor(hd$counts@x))) {
    stop("bin", size, " input is not finite nonnegative raw integer counts.", call. = FALSE)
  }
  coordinates <- as.matrix(hd$coordinates)
  assert_ids(rownames(coordinates), paste0("bin", size, " coordinate barcodes"))
  if (!setequal(rownames(coordinates), ids) || ncol(coordinates) != 2L ||
      !is.numeric(coordinates) || any(!is.finite(coordinates))) {
    stop("bin", size, " coordinates must be aligned, finite micrometre pairs.", call. = FALSE)
  }
  coordinates <- coordinates[ids, , drop = FALSE]
  metadata_ids <- if ("barcode" %in% names(hd$metadata)) as.character(hd$metadata$barcode) else rownames(hd$metadata)
  assert_ids(metadata_ids, paste0("bin", size, " metadata barcodes"))
  if (!setequal(metadata_ids, ids) || !"in_tissue" %in% names(hd$metadata)) {
    stop("bin", size, " lacks a complete keyed tissue-mask audit.", call. = FALSE)
  }
  metadata <- hd$metadata[match(ids, metadata_ids), , drop = FALSE]
  if (!all(c("array_row", "array_col") %in% names(metadata)) ||
      !is.numeric(metadata$array_row) || !is.numeric(metadata$array_col) ||
      any(!is.finite(metadata$array_row)) || any(!is.finite(metadata$array_col))) {
    stop("bin", size, " lacks finite array-grid positions for checking micrometre units.", call. = FALSE)
  }
  expected_xy <- cbind((metadata$array_col + .5) * size,
                        (metadata$array_row + .5) * size)
  if (any(abs(unname(coordinates) - expected_xy) > 1e-6)) {
    stop("bin", size, " coordinates do not match array-grid bin-centre micrometres.", call. = FALSE)
  }
  if (anyNA(metadata$in_tissue) || any(!metadata$in_tissue %in% c(0, 1))) {
    stop("bin", size, " tissue mask must contain zero or one.", call. = FALSE)
  }
  annotations <- fit$annotations
  if (!is.data.frame(annotations) ||
      !all(c("barcode", "status", "dominant_type", "max_weight") %in% names(annotations))) {
    stop("bin", size, " lacks the expected per-bin RCTD annotation audit.", call. = FALSE)
  }
  assert_ids(as.character(annotations$barcode), paste0("bin", size, " annotation barcodes"))
  if (!setequal(as.character(annotations$barcode), ids)) {
    stop("bin", size, " annotations must cover every input bin.", call. = FALSE)
  }
  annotations <- annotations[match(ids, annotations$barcode), , drop = FALSE]
  weights <- fit$weights
  if (is.null(dim(weights))) stop("bin", size, " lacks a weight matrix.", call. = FALSE)
  assert_ids(rownames(weights), paste0("bin", size, " weight barcodes"))
  assert_ids(colnames(weights), paste0("bin", size, " weight types"))
  if (!setequal(rownames(weights), ids)) stop("bin", size, " weight rows do not match input bins.", call. = FALSE)
  weights <- methods::as(Matrix::Matrix(weights[ids, , drop = FALSE], sparse = TRUE), "dgCMatrix")
  if (any(!is.finite(weights@x)) || any(weights@x < 0)) {
    stop("bin", size, " weights must be finite and nonnegative.", call. = FALSE)
  }
  weight_sum <- as.numeric(Matrix::rowSums(weights))
  fitted <- weight_sum > 0
  if (any(abs(weight_sum[fitted] - 1) > 1e-6)) {
    stop("bin", size, " nonzero weight rows are not normalized.", call. = FALSE)
  }
  umi <- as.numeric(Matrix::colSums(hd$counts))
  detected <- as.numeric(Matrix::colSums(hd$counts > 0))
  if ("nUMI" %in% names(annotations) &&
      (!is.numeric(annotations$nUMI) || any(!is.finite(annotations$nUMI)) ||
       any(abs(annotations$nUMI - umi) > 1e-6))) {
    stop("bin", size, " count totals differ from the RCTD input audit.", call. = FALSE)
  }
  if (any(fitted & umi <= 0)) stop("bin", size, " has positive fitted weights on a zero-UMI bin.", call. = FALSE)
  actual_max <- sparse_row_max(weights)
  if (!is.numeric(annotations$max_weight) || any(!is.finite(annotations$max_weight)) ||
      any(abs(annotations$max_weight - actual_max) > 1e-6)) {
    stop("bin", size, " max_weight disagrees with the recorded mixture.", call. = FALSE)
  }
  status <- as.character(annotations$status)
  if (anyNA(status) || any(!nzchar(status))) stop("Missing bin numerical status.", call. = FALSE)
  dominant <- as.character(annotations$dominant_type)
  labelled <- !is.na(dominant) & nzchar(dominant)
  if (any(!dominant[labelled] %in% colnames(weights))) {
    stop("bin", size, " dominant labels are absent from the fitted profiles.", call. = FALSE)
  }
  if (any(labelled & (!fitted | status != "ok"))) {
    stop("bin", size, " has dominant labels on bins without accepted numerical fits.", call. = FALSE)
  }
  if (any(labelled)) {
    at <- which(labelled)
    labelled_weight <- as.numeric(weights[cbind(at, match(dominant[at], colnames(weights)))])
    if (any(abs(labelled_weight - actual_max[at]) > 1e-6)) {
      stop("bin", size, " dominant labels disagree with the largest RNA contribution.", call. = FALSE)
    }
  }
  known_id <- all(genes %in% feature_info$id)
  known_name <- all(genes %in% feature_info$name)
  recorded_scheme <- hd$provenance$gene_column
  if (!is.null(recorded_scheme) && recorded_scheme %in% c("id", "name")) {
    scheme <- recorded_scheme
    mapped <- if (scheme == "id") genes %in% feature_info$id else genes %in% feature_info$name
    if (!all(mapped)) stop("bin", size, " genes conflict with its recorded ID/name scheme.", call. = FALSE)
  } else if (known_id && !known_name) scheme <- "id"
  else if (known_name && !known_id) scheme <- "name"
  else if (known_id && known_name) scheme <- "id_and_name_ambiguous"
  else stop("bin", size, " genes cannot be fully mapped to the retained feature ID/name audit.", call. = FALSE)
  gene_rows[[i]] <- data.frame(
    bin_size_um = size, identifier_scheme = scheme, expression_gene_rows = length(genes),
    retained_feature_ids = length(unique(feature_info$id)),
    retained_feature_symbols = length(unique(feature_info$name)),
    mapped_rows_by_id = sum(genes %in% feature_info$id),
    mapped_rows_by_name = sum(genes %in% feature_info$name),
    duplicate_symbols_in_feature_map = sum(duplicated(feature_info$name)),
    example_identifiers = paste(utils::head(genes, 5L), collapse = ";"),
    stringsAsFactors = FALSE
  )
  ref_audit <- fit$diagnostics$reference
  if (!is.null(ref_audit$input_cells) && ref_audit$input_cells != reference_input_cells) {
    stop("bin", size, " recorded reference input size differs from ovarian_reference.rds.", call. = FALSE)
  }
  if (!is.null(ref_audit$input_cells_per_type)) {
    observed <- ref_audit$input_cells_per_type
    positive <- observed[observed > 0]
    expected <- input_reference_counts[input_reference_counts > 0]
    if (!setequal(names(positive), names(expected)) ||
        any(as.numeric(positive[names(expected)]) != as.numeric(expected))) {
      stop("bin", size, " recorded reference type census differs from ovarian_reference.rds.", call. = FALSE)
    }
  }
  if (is.null(ref_audit$retained_cells) || is.null(ref_audit$retained_cells_per_type) ||
      is.null(names(ref_audit$retained_cells_per_type))) {
    stop("bin", size, " lacks actual Reference-constructor retention diagnostics; these cannot be inferred from the plotting reference or original input alone.", call. = FALSE)
  }
  retained_counts <- as.integer(ref_audit$retained_cells_per_type)
  names(retained_counts) <- names(ref_audit$retained_cells_per_type)
  if (anyNA(retained_counts) || any(retained_counts < 0) ||
      sum(retained_counts) != ref_audit$retained_cells) {
    stop("bin", size, " reference retention audit is inconsistent.", call. = FALSE)
  }
  type_universe <- union(union(names(input_reference_counts), names(retained_counts)), colnames(weights))
  reference_rows[[i]] <- data.frame(
    bin_size_um = size, reference_type = type_universe,
    input_reference_cells = as.integer(input_reference_counts[match(type_universe, names(input_reference_counts))]),
    constructor_retained_cells = as.integer(retained_counts[match(type_universe, names(retained_counts))]),
    included_in_fitted_profiles = type_universe %in% colnames(weights),
    reference_audit_stage = as.character(ref_audit$stage), stringsAsFactors = FALSE
  )
  # Only nonzero fitted rows enter this denominator. Do not convert RNA shares
  # into cell numbers or average fractions over a different observation mask.
  fitted_umi <- sum(umi[fitted])
  contribution <- rep(NA_real_, ncol(weights))
  if (fitted_umi > 0) {
    contribution <- as.numeric(Matrix::crossprod(weights[fitted, , drop = FALSE], umi[fitted])) / fitted_umi
    if (any(!is.finite(contribution)) || abs(sum(contribution) - 1) > 1e-6) {
      stop("bin", size, " UMI-weighted contribution summary failed mass balance.", call. = FALSE)
    }
  }
  expanded <- rep(NA_real_, length(type_universe))
  expanded[match(colnames(weights), type_universe)] <- contribution
  contribution_rows[[i]] <- data.frame(
    bin_size_um = size, reference_type = type_universe,
    umi_weighted_rna_contribution = expanded,
    included_in_fitted_profiles = type_universe %in% colnames(weights),
    fitted_umi_denominator = fitted_umi,
    fitted_bins = sum(fitted), input_umi = sum(umi),
    fitted_umi_fraction = if (sum(umi) > 0) fitted_umi / sum(umi) else NA_real_,
    units = "estimated RNA/UMI contribution; not cell-number fraction",
    stringsAsFactors = FALSE
  )
  filter_reason <- if ("filter_reason" %in% names(annotations)) as.character(annotations$filter_reason) else rep(NA_character_, length(ids))
  audit_state <- ifelse(umi == 0, "zero_UMI",
                         ifelse(fitted, paste0("mixture:", status), paste0("no_mixture:", status)))
  tile_key <- paste0("status:", status)
  tile_key[labelled] <- paste0("type:", dominant[labelled])
  tile_key[!fitted & status == "filtered"] <- "status:filtered"
  tile_key[!fitted & status == "not_fit"] <- "status:not_fit"
  tile_key[umi == 0] <- "status:zero_UMI"
  tile_key[metadata$in_tissue == 0] <- "status:outside_tissue"
  metrics[[i]] <- data.frame(
    barcode = ids, bin_size_um = size, x_um = coordinates[, 1L], y_um = coordinates[, 2L],
    umi = umi, detected_genes = detected, status = status, filter_reason = filter_reason,
    in_tissue = as.integer(metadata$in_tissue), has_mixture = fitted,
    dominant_type = dominant, max_weight = actual_max,
    audit_state = audit_state, tile_key = tile_key,
    stringsAsFactors = FALSE
  )
  cross_status <- data.frame(status = status,
                             filter_reason = ifelse(is.na(filter_reason), "none", filter_reason),
                             in_tissue = as.integer(metadata$in_tissue),
                             zero_umi = umi == 0, has_mixture = fitted,
                             n = 1L, stringsAsFactors = FALSE)
  counts <- stats::aggregate(n ~ status + filter_reason + in_tissue + zero_umi + has_mixture,
                             data = cross_status, FUN = sum)
  counts$bin_size_um <- size
  status_rows[[i]] <- counts[, c("bin_size_um", setdiff(names(counts), "bin_size_um"))]
  base_row <- list(
    bin_size_um = size, source_bin_size_um = hd$source_bin_size,
    gene_identifier_scheme = scheme, gene_rows = nrow(hd$counts), input_bins = length(ids),
    in_tissue_bins = sum(metadata$in_tissue == 1), outside_tissue_bins = sum(metadata$in_tissue == 0),
    positive_umi_bins = sum(umi > 0), zero_umi_bins = sum(umi == 0), total_umi = sum(umi),
    fitted_nonzero_weight_bins = sum(fitted), zero_weight_bins = sum(!fitted),
    dominant_rna_label_bins = sum(labelled), fitted_umi_denominator = fitted_umi,
    fitted_umi_fraction = if (sum(umi) > 0) fitted_umi / sum(umi) else NA_real_,
    reference_input_cells = reference_input_cells,
    reference_input_types_with_cells = sum(input_reference_counts > 0),
    reference_constructor_retained_cells = ref_audit$retained_cells,
    reference_constructor_types_with_cells = sum(retained_counts > 0),
    fitted_reference_types = ncol(weights),
    rctd_mode = fit$parameters$doublet_mode,
    spacexr_version = fit$diagnostics$spacexr_version,
    numerical_status_note = "ok denotes a valid numerical mixture, not calibrated confidence or biological accuracy"
  )
  for (state in c("ok", "ambiguous", "reject", "filtered", "not_fit")) {
    base_row[[paste0("status_", state, "_bins")]] <- sum(status == state)
  }
  base_row$status_other_bins <- sum(!status %in% c("ok", "ambiguous", "reject", "filtered", "not_fit"))
  summary_rows[[i]] <- as.data.frame(c(
    base_row, quantile_values(umi, "umi"), quantile_values(detected, "detected_genes"),
    quantile_values(actual_max[fitted], "max_rna_weight_fitted")
  ), stringsAsFactors = FALSE)
  rm(hd, fit, weights)
  invisible(gc())
}

summary_table <- do.call(rbind, summary_rows)
metric_table <- do.call(rbind, metrics)
status_table <- do.call(rbind, status_rows)
reference_table <- do.call(rbind, reference_rows)
contribution_table <- do.call(rbind, contribution_rows)
gene_table <- do.call(rbind, gene_rows)
metric_table$bin_label <- factor(paste0(metric_table$bin_size_um, " um"),
                                 levels = paste0(bin_sizes, " um"))
contribution_table$bin_label <- factor(paste0(contribution_table$bin_size_um, " um"),
                                       levels = paste0(bin_sizes, " um"))
plot_types <- unique(contribution_table$reference_type)
mean_contribution <- vapply(plot_types, function(type) {
  x <- contribution_table$umi_weighted_rna_contribution[contribution_table$reference_type == type]
  if (all(is.na(x))) -Inf else mean(x, na.rm = TRUE)
}, numeric(1))
plot_types <- plot_types[order(-mean_contribution, plot_types)]
contribution_table$reference_type <- factor(contribution_table$reference_type,
                                            levels = rev(plot_types))
contribution_table$value_label <- ifelse(
  is.na(contribution_table$umi_weighted_rna_contribution), "not fitted",
  sprintf("%.1f%%", 100 * contribution_table$umi_weighted_rna_contribution)
)

# Do not create partially reassuring figures until every required result and
# alignment/normalization/reference audit has passed format validation.
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(plot_dir)) stop("Cannot create QC plot directory.", call. = FALSE)
write_tsv(summary_table, outputs[[1L]])
write_tsv(status_table, outputs[[2L]])
write_tsv(reference_table, outputs[[3L]])
write_tsv(contribution_table, outputs[[4L]])
write_tsv(gene_table, outputs[[5L]])
writeLines(c(
  "Completed-result internal QC only. No independent truth labels are available to this script.",
  "RCTD was required to be full mode. Weights are estimated relative RNA/UMI contributions, not cell-number fractions.",
  "The heatmap uses sum(bin UMI * type weight) / sum(bin UMI), restricted to nonzero fitted rows at each resolution.",
  "Fit masks and covered UMI fractions may differ across resolutions; comparisons are conditional on those reported masks.",
  "Numerical status ok and max_weight are not calibrated classification probabilities or biological accuracy measurements.",
  "Dominant-type tiles depict capture bins, not segmented cells. No cell-cell contact or true cell count is inferred.",
  "Coordinate axes are array-grid micrometres; row coordinates increase downward. All map panels use common limits inferred from actual bin extents.",
  "Outside-tissue, zero-count, filtered, rejected and unfit states remain in the summaries and maps where present.",
  "Reference retention counts come from the actual Reference constructor before create.RCTD plotting-reference subsampling.",
  "No biological pass/fail threshold or marker-based accuracy test is imposed. Histology/marker/orthogonal validation remains separate."
), outputs[[6L]], useBytes = TRUE)

# Figure 1: depth, detected features and mutually exclusive recorded fit states.
open_png(outputs[[7L]], 3000, 2350)
tryCatch({
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    4L, 3L, heights = grid::unit(c(.15, 1, 1, 1.2), "null")
  )))
  grid::grid.text("ROI input depth and numerical fit retention (bins, not cells)",
                  vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1:3),
                  gp = grid::gpar(fontsize = 17, fontface = "bold"))
  for (i in seq_along(bin_sizes)) {
    size <- bin_sizes[[i]]
    df <- metric_table[metric_table$bin_size_um == size, , drop = FALSE]
    depth <- ggplot(df, aes(x = log10(umi + 1))) +
      geom_histogram(bins = 35, fill = "#377EB8", color = "white", linewidth = .1) +
      labs(title = paste0(size, " um: ", format(nrow(df), big.mark = ","), " input bins"),
           x = "log10(UMI + 1)", y = "Number of bins") + theme_minimal(base_size = 12)
    detected <- ggplot(df, aes(x = log10(detected_genes + 1))) +
      geom_histogram(bins = 35, fill = "#4DAF4A", color = "white", linewidth = .1) +
      labs(x = "log10(detected genes + 1)", y = "Number of bins") + theme_minimal(base_size = 12)
    retention <- as.data.frame(table(df$audit_state), stringsAsFactors = FALSE)
    names(retention) <- c("audit_state", "bins")
    retention$audit_state <- factor(retention$audit_state)
    retention_plot <- ggplot(retention, aes(x = "Input bins", y = bins, fill = audit_state)) +
      geom_col(width = .6) + geom_text(aes(label = bins), position = position_stack(vjust = .5), size = 3) +
      labs(title = paste0(sum(df$has_mixture), " nonzero fitted mixtures"),
           x = NULL, y = "Number of bins", fill = "Recorded state") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
      guides(fill = guide_legend(ncol = 1, byrow = TRUE))
    print(depth, vp = grid::viewport(layout.pos.row = 2L, layout.pos.col = i))
    print(detected, vp = grid::viewport(layout.pos.row = 3L, layout.pos.col = i))
    print(retention_plot, vp = grid::viewport(layout.pos.row = 4L, layout.pos.col = i))
  }
}, finally = grDevices::dev.off())

# Figure 2: fitted-UMI-weighted RNA contributions, preserving absent-model types.
heatmap <- ggplot(contribution_table, aes(x = bin_label, y = reference_type,
                                         fill = umi_weighted_rna_contribution)) +
  geom_tile(color = "white", linewidth = .5) +
  geom_text(aes(label = value_label,
                color = ifelse(umi_weighted_rna_contribution > .55, "light", "dark")),
             size = 3.3, show.legend = FALSE) +
  scale_color_manual(values = c(light = "white", dark = "#222222"), na.value = "#222222") +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 1),
                       na.value = "#BDBDBD", labels = function(x) paste0(round(100 * x), "%")) +
  scale_y_discrete(labels = function(x) vapply(x, function(v) paste(strwrap(v, width = 34), collapse = "\n"), character(1))) +
  labs(title = "UMI-weighted estimated RNA contributions by reference type",
       subtitle = "Full-mode mixtures; each column uses only nonzero fitted rows",
       x = "Bin width", y = "Reference type", fill = "RNA contribution",
       caption = "RNA contributions are not cell-number fractions. Covered UMI fractions and fit masks are reported separately.") +
  theme_minimal(base_size = 12) + theme(panel.grid = element_blank(), plot.caption = element_text(hjust = 0))
save_plot(heatmap, outputs[[8L]], 1800, max(1550, 85 * length(plot_types) + 350))

# Figure 3: observed bin footprints on common physical axes; no reconstructed cells.
x_limits <- range(metric_table$x_um - metric_table$bin_size_um / 2,
                   metric_table$x_um + metric_table$bin_size_um / 2)
y_limits <- range(metric_table$y_um - metric_table$bin_size_um / 2,
                   metric_table$y_um + metric_table$bin_size_um / 2)
type_keys <- paste0("type:", plot_types)
type_palette <- stats::setNames(grDevices::hcl.colors(length(type_keys), palette = "Dark 3"), type_keys)
state_keys <- setdiff(unique(metric_table$tile_key), type_keys)
state_palette <- stats::setNames(rep("#D9D9D9", length(state_keys)), state_keys)
known_states <- c("status:zero_UMI" = "#FAFAFA", "status:filtered" = "#BDBDBD",
                   "status:not_fit" = "#737373", "status:reject" = "#222222",
                   "status:ambiguous" = "#E6AB02", "status:outside_tissue" = "#FFFFFF",
                   "status:ok" = "#969696")
for (key in intersect(names(state_palette), names(known_states))) state_palette[[key]] <- known_states[[key]]
key_labels <- function(x) {
  answer <- sub("^type:", "", x)
  replacements <- c("status:zero_UMI" = "Zero UMI", "status:filtered" = "Filtered / no mixture",
                    "status:not_fit" = "Not fitted", "status:reject" = "Rejected",
                    "status:ambiguous" = "Ambiguous / no dominant label",
                    "status:outside_tissue" = "Outside tissue mask",
                    "status:ok" = "Fitted / no dominant label")
  at <- match(x, names(replacements))
  answer[!is.na(at)] <- replacements[at[!is.na(at)]]
  answer
}
map <- ggplot(metric_table, aes(x = x_um, y = y_um, fill = tile_key)) +
  geom_tile(aes(width = bin_size_um, height = bin_size_um)) +
  facet_wrap(~bin_label, nrow = 1, drop = FALSE) +
  scale_fill_manual(values = c(type_palette, state_palette), labels = key_labels,
                    drop = FALSE, na.value = "#D9D9D9") +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_reverse(limits = rev(y_limits), expand = c(0, 0)) + coord_fixed() +
  labs(title = "Dominant estimated RNA contributor per spatial bin",
       subtitle = "Capture bins, not segmented cells; these maps establish neither calibrated confidence nor biological accuracy",
       x = "Array-grid x (micrometres)", y = "Array-grid y (micrometres; rows increase downward)",
       fill = "RNA contributor / audit state",
       caption = "All panels share physical coordinate limits. White gaps may be absent input bins; explicit outside-tissue and zero-UMI bins are separately labelled.") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.text = element_text(size = 9),
         plot.caption = element_text(hjust = 0)) + guides(fill = guide_legend(ncol = 5, byrow = TRUE))
save_plot(map, outputs[[9L]], 3300, 1550)
message("QC tables and PNG figures saved under ", report_dir)
message("No biological accuracy, confidence threshold or cell-count claim was inferred.")
