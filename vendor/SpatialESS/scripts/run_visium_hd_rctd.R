#!/usr/bin/env Rscript

# Configurations and references are RDS data; this script never sources them.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1L && args %in% c("-h", "--help")) {
  cat(paste0(
    "Usage: Rscript scripts/run_visium_hd_rctd.R config.rds\n\n",
    "Run independent Visium HD RCTD fits for each configured bin size.\n",
    "See examples/visium_hd_rctd_config.R and docs/VISIUM_HD_RCTD.md.\n",
    "All relative paths are resolved from the configuration file directory.\n"
  ))
  quit(status = 0L)
}
if (length(args) != 1L) {
  stop("Expected one configuration RDS path. Use --help for usage.", call. = FALSE)
}

config_file <- normalizePath(args[[1L]], mustWork = TRUE)
config_dir <- dirname(config_file)
config <- readRDS(config_file)
if (!is.list(config) || is.null(names(config)) || anyNA(names(config)) ||
    any(!nzchar(names(config))) || anyDuplicated(names(config))) {
  stop("The configuration must be a uniquely named list.", call. = FALSE)
}
allowed <- c(
  "data_dir", "reference_file", "output_dir", "bin_sizes", "doublet_mode",
  "max_cores", "min_umi", "min_genes", "create_args", "run_args", "seed",
  "batch_size", "keep_fit", "read_args", "read_args_by_bin", "fit_args_by_bin", "overwrite",
  "save_inputs", "lr_file", "min_weight", "normalization", "spatialess_args"
)
unknown <- setdiff(names(config), allowed)
if (length(unknown)) {
  stop("Unknown configuration fields: ", paste(unknown, collapse = ", "), call. = FALSE)
}
for (field in c("data_dir", "reference_file", "output_dir")) {
  if (!is.character(config[[field]]) || length(config[[field]]) != 1L ||
      is.na(config[[field]]) || !nzchar(config[[field]])) {
    stop("Configuration field '", field, "' must be one path.", call. = FALSE)
  }
}
if (!is.null(config$lr_file) &&
    (!is.character(config$lr_file) || length(config$lr_file) != 1L ||
     is.na(config$lr_file) || !nzchar(config$lr_file))) {
  stop("Configuration field 'lr_file' must be one path.", call. = FALSE)
}
resolve_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(config_dir, path)
  path
}
for (field in c("data_dir", "reference_file", "output_dir", "lr_file")) {
  if (!is.null(config[[field]])) config[[field]] <- resolve_path(config[[field]])
}
for (field in c("overwrite", "save_inputs", "keep_fit")) {
  if (!is.null(config[[field]]) &&
      (!is.logical(config[[field]]) || length(config[[field]]) != 1L ||
       is.na(config[[field]]))) {
    stop("Configuration field '", field, "' must be TRUE or FALSE.", call. = FALSE)
  }
}
overwrite <- isTRUE(config$overwrite)
save_inputs <- if (is.null(config$save_inputs)) TRUE else config$save_inputs
bin_sizes <- if (is.null(config$bin_sizes)) c(8, 16, 32) else config$bin_sizes
if (!is.numeric(bin_sizes) || !length(bin_sizes) || anyNA(bin_sizes) ||
    any(!is.finite(bin_sizes)) || any(bin_sizes <= 0 | bin_sizes != floor(bin_sizes)) ||
    anyDuplicated(bin_sizes)) {
  stop("bin_sizes must be distinct positive integer micrometre sizes.", call. = FALSE)
}
check_named_list <- function(x, name) {
  if (is.null(x)) return(list())
  if (!is.list(x) || (length(x) &&
      (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
       anyDuplicated(names(x))))) {
    stop(name, " must be a uniquely named list.", call. = FALSE)
  }
  x
}
read_args <- check_named_list(config$read_args, "read_args")
per_bin <- check_named_list(config$read_args_by_bin, "read_args_by_bin")
bad_keys <- setdiff(names(per_bin), paste0("bin", bin_sizes))
if (length(bad_keys)) {
  stop("read_args_by_bin has unrequested bin keys: ", paste(bad_keys, collapse = ", "),
       call. = FALSE)
}
for (key in names(per_bin)) per_bin[[key]] <- check_named_list(per_bin[[key]], key)
read_options <- function(bin_size) {
  overrides <- per_bin[[paste0("bin", bin_size)]]
  if (is.null(overrides)) overrides <- list()
  result <- utils::modifyList(read_args, overrides)
  if (any(c("data_dir", "bin_size") %in% names(result))) {
    stop("Use top-level data_dir and bin_sizes; do not repeat them in read_args.",
         call. = FALSE)
  }
  if (!is.null(result$positions_file)) {
    result$positions_file <- resolve_path(result$positions_file)
  }
  result
}
fit_fields <- c(
  "doublet_mode", "max_cores", "min_umi", "min_genes", "create_args",
  "run_args", "seed", "batch_size", "keep_fit"
)
fit_args <- config[intersect(fit_fields, names(config))]
fit_args <- fit_args[!vapply(fit_args, is.null, logical(1L))]
if (is.null(fit_args$doublet_mode)) fit_args$doublet_mode <- "full"
fit_by_bin <- check_named_list(config$fit_args_by_bin, "fit_args_by_bin")
if (length(setdiff(names(fit_by_bin), paste0("bin", bin_sizes)))) {
  stop("fit_args_by_bin contains an unrequested bin key.", call. = FALSE)
}
for (key in names(fit_by_bin)) {
  fit_by_bin[[key]] <- check_named_list(fit_by_bin[[key]], key)
  if (length(setdiff(names(fit_by_bin[[key]]), fit_fields))) {
    stop("Unsupported fitting option in fit_args_by_bin$", key, ".", call. = FALSE)
  }
}
fit_options <- function(bin_size) {
  overrides <- fit_by_bin[[paste0("bin", bin_size)]]
  if (is.null(overrides)) overrides <- list()
  utils::modifyList(fit_args, overrides, keep.null = TRUE)
}

if (!requireNamespace("SpatialESS", quietly = TRUE)) {
  stop("Install this SpatialESS test build first: R CMD INSTALL .", call. = FALSE)
}
needed <- c("read_visium_hd", "run_visium_hd_rctd", "spatialess_from_visium_hd")
if (!all(needed %in% getNamespaceExports("SpatialESS"))) {
  stop("Installed SpatialESS lacks the Visium HD API. Install this test build first.",
       call. = FALSE)
}
if (!"feature_policy" %in% names(formals(SpatialESS::read_visium_hd))) {
  stop("Installed SpatialESS predates the scientific binning fixes; install version 0.1.4.9002 or newer.",
       call. = FALSE)
}
if (!dir.exists(config$data_dir)) stop("data_dir does not exist.", call. = FALSE)
reference <- readRDS(config$reference_file)
if (!is.list(reference) || !all(c("counts", "cell_types") %in% names(reference))) {
  stop("reference_file must contain list(counts = ..., cell_types = ...).", call. = FALSE)
}
lr_db <- NULL
if (!is.null(config$lr_file)) {
  lr_db <- readRDS(config$lr_file)
  if (!is.list(lr_db) || !all(c("lr", "complex", "cofactor") %in% names(lr_db))) {
    stop("lr_file must contain list(lr = ..., complex = ..., cofactor = ...).", call. = FALSE)
  }
}
spatialess_args <- check_named_list(config$spatialess_args, "spatialess_args")
reserved <- c("hd", "deconvolution", "lr", "complex", "cofactor",
              "min_weight", "normalization")
if (any(names(spatialess_args) %in% reserved)) {
  stop("spatialess_args repeats a reserved bridge argument.", call. = FALSE)
}

output_dir <- config$output_dir
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                               no.. = TRUE)) && !overwrite) {
  stop("output_dir is not empty. Choose a new directory or set overwrite = TRUE.",
       call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Cannot create output_dir.", call. = FALSE)
saveRDS(config, file.path(output_dir, "resolved_config.rds"))
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
summary_rows <- list()
write_summary <- function() {
  if (length(summary_rows)) {
    utils::write.table(do.call(rbind, summary_rows),
                       file.path(output_dir, "run_summary.tsv"),
                       sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
  }
}

for (bin_size in bin_sizes) {
  key <- paste0("bin", bin_size)
  bin_dir <- file.path(output_dir, key)
  dir.create(bin_dir, recursive = TRUE, showWarnings = FALSE)
  if (overwrite) {
    # Remove only this command's known files; leave unrelated user files alone.
    generated <- c("visium_hd_input.rds", "rctd_result.rds", "proportions.rds",
                   "annotations.tsv", "spatialess_result.rds", "spatialess_records.tsv",
                   "resolved_fit_args.rds", "ERROR.txt")
    unlink(file.path(bin_dir, generated))
  }
  message("[", key, "] Reading counts and tissue positions")
  start <- proc.time()[["elapsed"]]
  row <- data.frame(
    bin_size_um = bin_size, source_bin_size_um = NA_real_, n_genes = NA_integer_,
    n_bins = NA_integer_, n_fitted_bins = NA_integer_, n_annotated_bins = NA_integer_,
    total_umi = NA_real_, feature_policy = NA_character_, matrix_kind = NA_character_,
    doublet_mode = NA_character_, min_umi = NA_real_, min_genes = NA_real_,
    elapsed_seconds = NA_real_, status = "running",
    error = "", stringsAsFactors = FALSE
  )
  error <- tryCatch({
    hd <- do.call(SpatialESS::read_visium_hd,
                  c(list(data_dir = config$data_dir, bin_size = bin_size),
                    read_options(bin_size)))
    row$source_bin_size_um <- hd$source_bin_size
    row$n_genes <- nrow(hd$counts)
    row$n_bins <- ncol(hd$counts)
    row$total_umi <- sum(hd$counts)
    row$feature_policy <- hd$provenance$feature_policy
    row$matrix_kind <- hd$provenance$matrix_kind
    if (save_inputs) saveRDS(hd, file.path(bin_dir, "visium_hd_input.rds"))
    bin_fit_args <- fit_options(bin_size)
    saveRDS(bin_fit_args, file.path(bin_dir, "resolved_fit_args.rds"))
    message("[", key, "] Fitting RCTD independently (", bin_fit_args$doublet_mode, " mode)")
    result <- do.call(
      SpatialESS::run_visium_hd_rctd,
      c(list(hd = hd, reference_counts = reference$counts,
             cell_types = reference$cell_types), bin_fit_args)
    )
    saveRDS(result, file.path(bin_dir, "rctd_result.rds"))
    saveRDS(result$weights, file.path(bin_dir, "proportions.rds"))
    utils::write.table(result$annotations, file.path(bin_dir, "annotations.tsv"),
                       sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
    row$n_fitted_bins <- sum(Matrix::rowSums(result$weights) > 0)
    row$n_annotated_bins <- sum(!is.na(result$annotations$dominant_type))
    row$doublet_mode <- result$parameters$doublet_mode
    row$min_umi <- result$parameters$min_umi
    row$min_genes <- result$parameters$min_genes
    if (!is.null(lr_db)) {
      message("[", key, "] Running SpatialESS with dominant-bin labels")
      bridge_args <- list(
        hd = hd, deconvolution = result, lr = lr_db$lr,
        complex = lr_db$complex, cofactor = lr_db$cofactor,
        min_weight = if (is.null(config$min_weight)) 0.5 else config$min_weight,
        normalization = if (is.null(config$normalization)) "log1p_cp10k" else config$normalization
      )
      communication <- do.call(SpatialESS::spatialess_from_visium_hd,
                               c(bridge_args, spatialess_args))
      saveRDS(communication, file.path(bin_dir, "spatialess_result.rds"))
      utils::write.table(communication$records, file.path(bin_dir, "spatialess_records.tsv"),
                         sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
    }
    NULL
  }, error = function(e) e)
  row$elapsed_seconds <- proc.time()[["elapsed"]] - start
  row$status <- if (is.null(error)) "ok" else "failed"
  if (!is.null(error)) row$error <- conditionMessage(error)
  summary_rows[[key]] <- row
  write_summary()
  if (!is.null(error)) {
    writeLines(conditionMessage(error), file.path(bin_dir, "ERROR.txt"))
    stop("[", key, "] Failed: ", conditionMessage(error),
         ". See run_summary.tsv; later resolutions were not run.", call. = FALSE)
  }
  # Release the preceding fit before reading another HD resolution.
  rm(hd, result)
  if (exists("communication", inherits = FALSE)) rm(communication)
  invisible(gc())
  message("[", key, "] Saved to ", bin_dir)
}
message("Completed all requested resolutions. Summary: ",
        file.path(output_dir, "run_summary.tsv"))
