# Shared SpatialESS per-dataset pipeline library.
# Sourced by workspace/<dataset>/run_spatialess_lmm.R drivers.
# Requires conda env "spatialess" (SpatialESS, lme4, SeuratObject, rhdf5).

# ---------------------------------------------------------------------------
# CLI parsing: --key value pairs -> named list
# ---------------------------------------------------------------------------
pipeline_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key, call. = FALSE)
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[name]] <- TRUE
    } else {
      out[[name]] <- args[[i + 1L]]
      i <- i + 1L
    }
    i <- i + 1L
  }
  out
}

arg_value <- function(args, name, default) {
  if (!name %in% names(args)) return(default)
  args[[name]]
}

arg_integer <- function(args, name, default) {
  if (!name %in% names(args)) return(default)
  value <- suppressWarnings(as.integer(args[[name]]))
  if (is.na(value) || value < 1L) {
    stop("--", name, " must be a positive integer.", call. = FALSE)
  }
  value
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log_event <- function(...) {
  cat(format(Sys.time(), "%FT%T%z"), ..., "\n")
  flush.console()
}

# ---------------------------------------------------------------------------
# CellChatDB: download once, verify, load
# ---------------------------------------------------------------------------
CELLCHATDB_URL <- paste0(
  "https://raw.githubusercontent.com/sqjin/CellChat/master/",
  "data/CellChatDB.human.rda"
)
CELLCHATDB_BYTES <- 756442L

ensure_cellchatdb <- function(cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  rda <- file.path(cache_dir, "CellChatDB.human.rda")
  if (!file.exists(rda)) {
    log_event("Downloading CellChatDB.human.rda ->", rda)
    download.file(CELLCHATDB_URL, destfile = rda, mode = "wb")
  }
  actual <- file.size(rda)
  if (actual != CELLCHATDB_BYTES) {
    stop(
      sprintf(
        "CellChatDB.human.rda size mismatch: got %d bytes, expected %d. Delete %s and retry.",
        actual, CELLCHATDB_BYTES, rda
      ),
      call. = FALSE
    )
  }
  env <- new.env(parent = emptyenv())
  load(rda, envir = env)
  db <- env$CellChatDB.human
  required <- c("interaction", "complex", "cofactor")
  if (!all(required %in% names(db))) {
    stop("CellChatDB.human lacks required tables: ",
         paste(setdiff(required, names(db)), collapse = ", "), call. = FALSE)
  }
  log_event("CellChatDB.human:", nrow(db$interaction), "interactions")
  db
}

# ---------------------------------------------------------------------------
# Seurat-object input (GEO RDS files)
# ---------------------------------------------------------------------------
# meta_cols: named list with entries section/cell_type/patient/condition and
# character vector covariates. Each entry is a character vector of candidate
# column names; the first candidate present in the object metadata wins and
# the mapping is logged. Entries may be NULL when not applicable.
detect_column <- function(candidates, columns, label, required = TRUE) {
  if (is.null(candidates)) return(NULL)
  hit <- candidates[candidates %in% columns]
  if (length(hit)) {
    log_event("Metadata mapping:", label, "->", hit[1L])
    return(hit[1L])
  }
  if (required) {
    stop("No candidate column for '", label, "' found. Candidates: ",
         paste(candidates, collapse = ", "), ". Available columns: ",
         paste(columns, collapse = ", "), call. = FALSE)
  }
  NULL
}

load_seurat_dataset <- function(rds_path, meta_cols) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required: conda install -n spatialess ",
         "-c conda-forge r-seuratobject", call. = FALSE)
  }
  log_event("Loading", rds_path)
  obj <- readRDS(rds_path)
  meta <- obj[[]]
  log_event("Cells:", nrow(meta), " metadata columns:",
            paste(colnames(meta), collapse = ", "))

  columns <- colnames(meta)
  meta_cols$section <- detect_column(meta_cols$section, columns, "section")
  meta_cols$cell_type <- detect_column(meta_cols$cell_type, columns, "cell_type")
  meta_cols$patient <- detect_column(meta_cols$patient, columns, "patient",
                                     required = FALSE)
  meta_cols$condition <- detect_column(meta_cols$condition, columns,
                                       "condition", required = FALSE)
  if (length(meta_cols$covariates)) {
    meta_cols$covariates <- Filter(
      function(cc) cc %in% columns,
      vapply(meta_cols$covariates, detect_column, character(1L),
             columns = columns, label = "covariate", required = FALSE,
             USE.NAMES = FALSE)
    )
  }
  for (label in c("section", "cell_type")) {
    if (is.null(meta_cols[[label]])) {
      stop("Required metadata '", label, "' could not be detected.",
           call. = FALSE)
    }
  }

  assay <- SeuratObject::DefaultAssay(obj)
  counts <- tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) SeuratObject::GetAssayData(obj, assay = assay,
                                                   slot = "counts")
  )
  counts <- methods::as(counts, "dgCMatrix")
  log_event("Counts:", nrow(counts), "genes x", ncol(counts), "cells, assay:",
            assay)

  coordinates <- .extract_coordinates(obj, meta)
  if (is.null(rownames(coordinates))) {
    rownames(coordinates) <- colnames(counts)
  }
  if (!identical(rownames(coordinates), colnames(counts))) {
    stop("Coordinate cell order differs from expression cell order.",
         call. = FALSE)
  }

  # Canonical column names downstream: section, cell_type, patient_id,
  # condition; covariates keep their original names.
  selected <- data.frame(
    section = meta[[meta_cols$section]],
    cell_type = meta[[meta_cols$cell_type]],
    stringsAsFactors = FALSE
  )
  if (!is.null(meta_cols$patient)) {
    selected$patient_id <- as.character(meta[[meta_cols$patient]])
  }
  if (!is.null(meta_cols$condition)) {
    selected$condition <- as.character(meta[[meta_cols$condition]])
  }
  for (cc in meta_cols$covariates) {
    selected[[cc]] <- as.character(meta[[cc]])
  }
  rownames(selected) <- rownames(meta)
  list(
    counts = counts,
    coordinates = coordinates,
    metadata = selected
  )
}

.extract_coordinates <- function(obj, meta) {
  candidates <- list(
    c("x_centroid", "y_centroid"),
    c("x_coord", "y_coord"),
    c("center_x", "center_y"),
    c("x", "y")
  )
  for (cc in candidates) {
    if (all(cc %in% colnames(meta))) {
      log_event("Coordinates from metadata columns:", paste(cc, collapse = ", "))
      return(as.matrix(meta[, cc]))
    }
  }
  log_event("No coordinate columns in metadata; trying FOV centroids")
  for (image in SeuratObject::Images(obj)) {
    coords <- tryCatch(
      SeuratObject::GetTissueCoordinates(obj[[image]]),
      error = function(e) NULL
    )
    if (is.null(coords)) {
      fov <- obj@images[[image]]
      coords <- tryCatch(fov@centroids@coords, error = function(e) NULL)
    }
    if (!is.null(coords) && nrow(coords)) {
      log_event("Coordinates from image:", image)
      return(as.matrix(coords[, seq_len(2L), drop = FALSE]))
    }
  }
  stop("Could not locate cell coordinates (metadata or FOV).", call. = FALSE)
}

# ---------------------------------------------------------------------------
# 10x Xenium h5 + cells.csv.gz input
# ---------------------------------------------------------------------------
load_xenium_dataset <- function(h5_path, cells_csv) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required: conda install -n spatialess ",
         "-c conda-forge -c bioconda r-rhdf5", call. = FALSE)
  }
  log_event("Reading h5:", h5_path)
  data <- rhdf5::h5read(h5_path, "matrix/data")
  indices <- rhdf5::h5read(h5_path, "matrix/indices")
  indptr <- rhdf5::h5read(h5_path, "matrix/indptr")
  shape <- rhdf5::h5read(h5_path, "matrix/shape")
  barcodes <- rhdf5::h5read(h5_path, "matrix/barcodes")
  feature_names <- rhdf5::h5read(h5_path, "matrix/features/name")
  if (is.list(feature_names)) feature_names <- feature_names[[1L]]
  counts <- Matrix::sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = as.integer(shape),
    dimnames = list(as.character(feature_names), as.character(barcodes))
  )
  is_control <- grepl(
    "^(NegControlProbe|NegControlCodeword|UnassignedCodeword|DeprecatedCodeword)",
    rownames(counts)
  )
  log_event("Dropped", sum(is_control), "negative/deprecated control features")
  counts <- counts[!is_control, , drop = FALSE]

  log_event("Reading cell centroids:", cells_csv)
  cells <- utils::read.csv(gzfile(cells_csv), stringsAsFactors = FALSE)
  match <- match(colnames(counts), cells$cell_id)
  if (anyNA(match)) {
    stop("cells.csv cell_id does not cover all h5 barcodes.", call. = FALSE)
  }
  coordinates <- as.matrix(cells[match, c("x_centroid", "y_centroid")])
  rownames(coordinates) <- colnames(counts)
  list(counts = methods::as(counts, "dgCMatrix"),
       coordinates = coordinates,
       metadata = NULL)
}

# ---------------------------------------------------------------------------
# Per-section spatialess inference with checkpointing
# ---------------------------------------------------------------------------
run_sections <- function(dataset, db, section_levels, out_dir,
                         max_cells = Inf, nboot = 100L, seed = 1L,
                         section_label = "section",
                         cell_type_label = "cell_type") {
  sections_dir <- file.path(out_dir, "sections")
  dir.create(sections_dir, recursive = TRUE, showWarnings = FALSE)
  stopifnot(!is.null(dataset$metadata))
  section_values <- as.character(dataset$metadata[[section_label]])
  cell_types <- as.character(dataset$metadata[[cell_type_label]])
  cell_types[is.na(cell_types) | !nzchar(cell_types)] <- "Unknown"

  results <- list()
  for (i in seq_along(section_levels)) {
    section <- section_levels[[i]]
    rds_file <- file.path(sections_dir, paste0(make.names(section), ".rds"))
    if (file.exists(rds_file)) {
      log_event("[", i, "/", length(section_levels), "] resume:", section)
      results[[section]] <- readRDS(rds_file)
      next
    }
    cells <- which(section_values == section)
    if (is.finite(max_cells) && length(cells) > max_cells) {
      set.seed(seed + i)
      cells <- sort(sample(cells, max_cells))
    }
    log_event("[", i, "/", length(section_levels), "] section", section,
              "cells:", length(cells), "nboot:", nboot)
    result <- tryCatch(
      SpatialESS::spatialess(
        dataset$counts[, cells, drop = FALSE],
        dataset$coordinates[cells, , drop = FALSE],
        factor(cell_types[cells]),
        lr = db$interaction,
        complex = db$complex,
        cofactor = db$cofactor,
        nboot = nboot,
        seed.use = seed + i,
        lr_missing = "filter"
      ),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      log_event("Section", section, "FAILED:", conditionMessage(result))
      next
    }
    log_event("Section", section, "records:", nrow(result$records),
              "LR used:", nrow(result$lr))
    saveRDS(result, rds_file)
    results[[section]] <- result
  }
  results
}

# ---------------------------------------------------------------------------
# Section-level metadata assembly
# ---------------------------------------------------------------------------
# metadata: per-cell canonical metadata from load_seurat_dataset().
# sections: character vector of section levels to include (order preserved).
# covariates: additional per-section constant columns to retain.
# Returns one row per section: sample_id, patient_id (when available),
# condition (when available), plus covariates. Patient and condition must be
# constant within a section.
build_section_metadata <- function(metadata, sections, covariates = character()) {
  stopifnot("section" %in% colnames(metadata))
  per_section <- function(column) {
    if (is.null(column) || !column %in% colnames(metadata)) return(NULL)
    vapply(sections, function(section) {
      values <- unique(metadata[[column]][metadata$section == section])
      if (length(values) != 1L) {
        stop("Column '", column, "' is not constant within section '",
             section, "'.", call. = FALSE)
      }
      as.character(values)
    }, character(1L))
  }
  out <- data.frame(
    sample_id = sections,
    patient_id = per_section("patient_id"),
    condition = per_section("condition"),
    stringsAsFactors = FALSE
  )
  for (cc in covariates) {
    values <- per_section(cc)
    if (!is.null(values)) out[[cc]] <- values
  }
  out
}

# ---------------------------------------------------------------------------
# Cohort modelling
# ---------------------------------------------------------------------------
write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = TRUE,
                     na = "NA")
  log_event("Wrote", path)
}

# section_metadata: one row per section with at least sample_id, patient_id,
# condition and any covariates. model is "lmm", "glm", or both.
run_cohort_models <- function(results, section_metadata, out_dir,
                              design = ~ condition,
                              coefficient,
                              patient_col = "patient_id",
                              model = c("lmm", "glm"),
                              covariates = character(),
                              prefix = "cohort",
                              min_units = 6L, min_nonzero = 2L) {
  model <- match.arg(model, several.ok = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  prepared <- SpatialESS::prepare_multisample_communication(
    results,
    section_metadata,
    unit = "sample",
    patient_col = patient_col,
    covariates = covariates
  )
  log_event("Prepared score matrix:", nrow(prepared$score), "sections x",
            ncol(prepared$score), "features")
  saveRDS(prepared, file.path(out_dir, paste0(prefix, "_prepared.rds")))

  tables <- list()
  if ("lmm" %in% model) {
    tables$lmm <- SpatialESS::fit_multisample_communication_lmm(
      prepared, design, coefficient, patient_col = patient_col,
      min_units = min_units, min_nonzero = min_nonzero
    )
    write_tsv(tables$lmm, file.path(out_dir, paste0(prefix, "_lmm.tsv")))
  }
  if ("glm" %in% model) {
    tables$glm <- SpatialESS::fit_multisample_communication_glm(
      prepared, design, coefficient,
      min_units = min_units, min_nonzero = min_nonzero
    )
    write_tsv(tables$glm, file.path(out_dir, paste0(prefix, "_glm.tsv")))
  }
  for (kind in names(tables)) {
    tab <- tables[[kind]]
    status <- as.data.frame(table(tab$status), stringsAsFactors = FALSE)
    names(status) <- c("status", "features")
    write_tsv(status, file.path(out_dir, paste0(prefix, "_", kind,
                                                "_fit_status.tsv")))
    ok <- tab$status == "ok"
    log_event(kind, "features:", nrow(tab), "ok:", sum(ok),
              "significant (q<0.05):",
              sum(ok & !is.na(tab$q_value) & tab$q_value < 0.05))
  }
  tables
}

write_run_summary <- function(out_dir, rows) {
  write_tsv(as.data.frame(rows, check.names = FALSE),
            file.path(out_dir, "run_summary.tsv"))
}
