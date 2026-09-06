#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
  library(SpatialESS)
})

env_number <- function(name, default, integer = FALSE) {
  value <- Sys.getenv(name, unset = as.character(default))
  parsed <- if (integer) suppressWarnings(as.integer(value)) else suppressWarnings(as.double(value))
  if (length(parsed) != 1L || is.na(parsed)) {
    stop(sprintf("%s must be numeric; received '%s'.", name, value), call. = FALSE)
  }
  parsed
}

data_root <- Sys.getenv("SECACT_ROOT", "/data/dzf/SecAct_2026")
archive_path <- Sys.getenv("SECACT_ARCHIVE", file.path(data_root, "3_Application.zip"))
data_path <- Sys.getenv(
  "SECACT_DATA",
  file.path(data_root, "LIHC_CosMx_data.rda")
)
nperm <- env_number("SPATIALESS_NPERM", 100L, integer = TRUE)
seed <- env_number("SPATIALESS_SEED", 20260731L, integer = TRUE)
radius_um <- env_number("SPATIALESS_RADIUS_UM", 20)
block_um <- env_number("SPATIALESS_BLOCK_UM", 100)
scale_factor <- env_number("SPATIALESS_SCALE_FACTOR", 1000)
min_genes <- env_number("SPATIALESS_MIN_GENES", 50L, integer = TRUE)
out_dir <- Sys.getenv(
  "SPATIALESS_OUT",
  file.path(data_root, "results", sprintf("secact_cosmx_radius%g_perm%d", radius_um, nperm))
)

if (nperm < 1L || seed < 0L || radius_um <= 0 || block_um <= radius_um ||
    scale_factor <= 0 || min_genes < 0L) {
  stop("Invalid benchmark settings: require nperm >= 1, seed >= 0, radius > 0, block > radius, scale factor > 0 and min_genes >= 0.", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  if (!file.exists(archive_path)) {
    stop(
      sprintf(
        "CosMx data not found. Place 3_Application.zip at %s or LIHC_CosMx_data.rda at %s.",
        archive_path, data_path
      ),
      call. = FALSE
    )
  }
  archive_files <- utils::unzip(archive_path, list = TRUE)$Name
  candidates <- archive_files[grepl("(^|/)LIHC_CosMx_data\\.rda$", archive_files)]
  if (length(candidates) != 1L) {
    stop(sprintf("Expected one LIHC_CosMx_data.rda in the archive; found %d.", length(candidates)), call. = FALSE)
  }
  message("Extracting only the HCC CosMx data from 3_Application.zip...")
  utils::unzip(archive_path, files = candidates, exdir = data_root)
  extracted <- file.path(data_root, candidates)
  if (!file.exists(extracted)) stop("CosMx extraction failed.", call. = FALSE)
  if (!identical(normalizePath(extracted), normalizePath(data_path, mustWork = FALSE))) {
    dir.create(dirname(data_path), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(extracted, data_path, overwrite = FALSE)) {
      stop(sprintf("Could not place extracted data at %s.", data_path), call. = FALSE)
    }
  }
}

stage <- list()
workflow_start <- proc.time()[["elapsed"]]
message("Loading SecAct HCC CosMx data...")
load_start <- proc.time()[["elapsed"]]
data_env <- new.env(parent = emptyenv())
loaded <- load(data_path, envir = data_env)
required <- c("counts", "spotCoordinates", "metaData")
if (!all(required %in% loaded)) {
  stop(sprintf("CosMx RDA is missing: %s.", paste(setdiff(required, loaded), collapse = ", ")), call. = FALSE)
}
counts <- data_env$counts
coords <- as.matrix(data_env$spotCoordinates)
metadata <- as.data.frame(data_env$metaData, stringsAsFactors = FALSE)
rm(data_env)
gc()

if (!inherits(counts, "Matrix")) counts <- Matrix(counts, sparse = TRUE)
counts <- as(counts, "dgCMatrix")
storage.mode(coords) <- "double"
if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
    anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
  stop("Counts require unique gene row names and cell column names.", call. = FALSE)
}
if (ncol(coords) < 2L || any(!is.finite(coords[, 1:2, drop = FALSE]))) {
  stop("spotCoordinates must contain at least two finite coordinate columns.", call. = FALSE)
}
coords <- coords[, 1:2, drop = FALSE]
if (is.null(rownames(coords)) || is.null(rownames(metadata))) {
  stop("Coordinates and metadata require cell row names.", call. = FALSE)
}
cell_id <- colnames(counts)
if (!all(cell_id %in% rownames(coords)) || !all(cell_id %in% rownames(metadata))) {
  stop("Counts, coordinates and metadata do not contain the same cells.", call. = FALSE)
}
coords <- coords[cell_id, , drop = FALSE]
metadata <- metadata[cell_id, , drop = FALSE]
if (!"cellType" %in% colnames(metadata)) {
  stop("SecAct CosMx metadata does not contain cellType.", call. = FALSE)
}

detected_genes <- Matrix::colSums(counts != 0)
library_size <- Matrix::colSums(counts)
keep <- detected_genes >= min_genes & library_size > 0 &
  !is.na(metadata$cellType) & nzchar(as.character(metadata$cellType))
if (!any(keep)) stop("No cells passed the CosMx quality-control filters.", call. = FALSE)
counts <- counts[, keep, drop = FALSE]
coords <- coords[keep, , drop = FALSE]
metadata <- metadata[keep, , drop = FALSE]
detected_genes <- detected_genes[keep]
library_size <- library_size[keep]

normalization <- counts %*% Diagonal(x = scale_factor / library_size)
dimnames(normalization) <- dimnames(counts)
normalization@x <- log1p(normalization@x)
expression <- drop0(as(normalization, "dgCMatrix"))
rm(normalization, counts)
gc()
stage$load_and_normalize_seconds <- proc.time()[["elapsed"]] - load_start

group <- factor(as.character(metadata$cellType), levels = sort(unique(as.character(metadata$cellType))))
sample_id <- rep("CancerousLiver", nrow(coords))
niche <- if ("niche" %in% colnames(metadata)) {
  value <- as.character(metadata$niche)
  value[is.na(value) | !nzchar(value)] <- "unassigned"
  value
} else {
  rep("all", nrow(coords))
}

message(sprintf("Building the %.1f um CSR radius graph for %s cells...", radius_um, format(ncol(expression), big.mark = ",")))
graph_start <- proc.time()[["elapsed"]]
graph <- build_radius_graph_csr(
  coords, radius = radius_um, sample_id = sample_id,
  weight = "binary", store_distance = FALSE, max_edges = 5e8
)
graph$cells$group_code <- as.integer(group)
graph$cells$group_levels <- levels(group)
checks <- validate_spatial_csr_graph(graph)
if (!isTRUE(attr(checks, "valid"))) stop("CSR graph validation failed.", call. = FALSE)
group_support <- build_group_support_csr(graph, group, max_group_pairs = 1e7)
stage$graph_and_support_seconds <- proc.time()[["elapsed"]] - graph_start

data(CellChatDB.human, package = "CellChat")
db <- CellChatDB.human
complex_columns <- grep("^subunit", colnames(db$complex), value = TRUE)
entity_available <- function(entity) {
  entity <- as.character(entity)
  if (is.na(entity) || !nzchar(entity)) return(FALSE)
  if (!entity %in% rownames(db$complex)) return(entity %in% rownames(expression))
  subunits <- unlist(db$complex[entity, complex_columns, drop = FALSE], use.names = FALSE)
  subunits <- as.character(subunits)
  subunits <- subunits[!is.na(subunits) & nzchar(subunits)]
  length(subunits) > 0L && all(subunits %in% rownames(expression))
}
lr_all <- db$interaction
valid_lr <- vapply(lr_all$ligand, entity_available, logical(1)) &
  vapply(lr_all$receptor, entity_available, logical(1))
lr <- lr_all[valid_lr, , drop = FALSE]
if (!nrow(lr)) stop("No CellChatDB ligand-receptor records are computable on this CosMx panel.", call. = FALSE)
components <- prepare_cellchat_lr_components(
  lr, rownames(expression), db$complex, db$cofactor
)

message(sprintf("Scoring %d computable CellChatDB interactions across %d cell types...", nrow(lr), nlevels(group)))
score_start <- proc.time()[["elapsed"]]
prepared <- prepare_sparse_trimean(expression, components$genes, normalize = TRUE)
group_expression <- summarize_prepared_trimean(
  prepared, group, group_levels = group_support$group_levels
)
observed <- score_cellchat_group_support(
  group_expression, group_support, components,
  Kh = 0.5, n = 1, max_output_records = 1e8
)
stage$observed_score_seconds <- proc.time()[["elapsed"]] - score_start

message(sprintf("Running %d profile permutations within sample x niche x %.1f um blocks...", nperm, block_um))
permutation_start <- proc.time()[["elapsed"]]
blocks <- build_spatial_blocks(
  coords, block_size = block_um, sample_id = sample_id, compartment = niche
)
permutation <- permutation_cellchat_group_support(
  prepared, group, group_support, components, blocks,
  nperm = nperm, seed = seed, Kh = 0.5, n = 1,
  finite_correction = TRUE, retain_null_moments = FALSE,
  tail = "greater_equal",
  max_score_entries = 1e8, verbose = TRUE
)
stage$permutation_seconds <- proc.time()[["elapsed"]] - permutation_start
stage$end_to_end_seconds <- proc.time()[["elapsed"]] - workflow_start

summary <- data.frame(
  source_article_doi = "10.1038/s41592-026-03172-0",
  source_data_url = "https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/",
  sample = "CancerousLiver",
  input_genes = nrow(expression),
  retained_cells = ncol(expression),
  median_detected_genes = stats::median(detected_genes),
  cell_types = nlevels(group),
  radius_um = radius_um,
  graph_edges = length(graph$neighbors),
  supported_group_pairs = nrow(group_support$group_pairs),
  computable_lr = nrow(lr),
  complex_lr = sum(lengths(components$ligand) > 1L | lengths(components$receptor) > 1L),
  referenced_genes = length(components$genes),
  block_um = block_um,
  spatial_blocks = length(blocks$block_levels),
  exchangeable_cell_fraction = permutation$diagnostics$exchangeable_cell_fraction,
  nperm = nperm,
  significant_records_0_05 = sum(permutation$group_pairs$pvalue <= 0.05),
  observed_active_records = nrow(permutation$group_pairs),
  load_and_normalize_seconds = stage$load_and_normalize_seconds,
  graph_and_support_seconds = stage$graph_and_support_seconds,
  observed_score_seconds = stage$observed_score_seconds,
  permutation_seconds = stage$permutation_seconds,
  end_to_end_seconds = stage$end_to_end_seconds,
  stringsAsFactors = FALSE
)

saveRDS(graph, file.path(out_dir, "spatial_graph_csr.rds"), compress = FALSE)
saveRDS(group_support, file.path(out_dir, "group_support.rds"), compress = FALSE)
saveRDS(components, file.path(out_dir, "cellchat_lr_components.rds"), compress = FALSE)
saveRDS(observed, file.path(out_dir, "observed_cellchat_compat.rds"), compress = FALSE)
saveRDS(permutation, file.path(out_dir, "spatial_permutation.rds"), compress = FALSE)
write.table(summary, file.path(out_dir, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(table(cell_type = group)), file.path(out_dir, "cell_type_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
significance_connection <- gzfile(
  file.path(out_dir, "significance_all.tsv.gz"), open = "wt"
)
write.table(
  permutation$group_pairs, significance_connection,
  sep = "\t", quote = FALSE, row.names = FALSE
)
close(significance_connection)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(
  c(
    sprintf("data_file=%s", normalizePath(data_path)),
    sprintf("data_md5=%s", unname(tools::md5sum(data_path))),
    sprintf("normalization=log1p(count / cell_library_size * %g)", scale_factor),
    sprintf("cell_qc=detected_genes >= %d and nonzero library size", min_genes),
    "group=SecAct metaData$cellType",
    sprintf("graph=undirected binary CSR radius graph, radius %g um", radius_um),
    sprintf("null=complete profiles permuted within CancerousLiver x niche x %g um block", block_um),
    "pvalue=(greater_or_equal_null_count + 1) / (nperm + 1)"
  ),
  file.path(out_dir, "provenance.txt")
)
print(summary)
