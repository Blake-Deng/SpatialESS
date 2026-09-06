#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})

input_path <- Sys.getenv(
  "SECACT_BRUKER_RDS",
  "/data/dzf/SecAct_2026/LiverDataReleaseSeurat_newUMAP.RDS"
)
output_path <- Sys.getenv(
  "SECACT_DATA",
  "/data/dzf/SecAct_2026/LIHC_CosMx_data.rda"
)
sample_column <- Sys.getenv("SECACT_SAMPLE_COLUMN", "Run_Tissue_name")
requested_sample <- Sys.getenv("SECACT_SAMPLE", "")

if (!file.exists(input_path)) {
  stop(sprintf("Official Bruker Seurat object not found: %s", input_path), call. = FALSE)
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

message("Reading the official Bruker/NanoString liver Seurat object...")
object <- readRDS(input_path)
if (!inherits(object, "Seurat")) {
  stop(sprintf("Expected a Seurat object; received: %s", paste(class(object), collapse = ", ")), call. = FALSE)
}

metadata <- object@meta.data
required_metadata <- c(sample_column, "niche", "cellType", "x_slide_mm", "y_slide_mm")
if (!all(required_metadata %in% colnames(metadata))) {
  stop(
    sprintf(
      "Seurat metadata is missing: %s",
      paste(setdiff(required_metadata, colnames(metadata)), collapse = ", ")
    ),
    call. = FALSE
  )
}

sample_values <- as.character(metadata[[sample_column]])
sample_counts <- sort(table(sample_values, useNA = "ifany"), decreasing = TRUE)
message("Samples in the official object:")
print(sample_counts)

available <- unique(sample_values[!is.na(sample_values) & nzchar(sample_values)])
if (nzchar(requested_sample)) {
  target <- available[tolower(available) == tolower(requested_sample)]
} else {
  target <- available[grepl("cancer|carcinoma|tumou?r", available, ignore.case = TRUE)]
}
if (length(target) != 1L) {
  stop(
    sprintf(
      "Could not identify exactly one cancer sample in %s. Candidates: %s",
      sample_column, paste(target, collapse = ", ")
    ),
    call. = FALSE
  )
}
keep <- !is.na(sample_values) & sample_values == target
if (!any(keep)) stop("The selected cancer sample contains no cells.", call. = FALSE)

rna <- object[["RNA"]]
if (is.null(rna)) stop("The Seurat object does not contain an RNA assay.", call. = FALSE)
if ("counts" %in% slotNames(rna)) {
  counts_all <- methods::slot(rna, "counts")
} else {
  counts_all <- LayerData(rna, layer = "counts")
}
if (!inherits(counts_all, "Matrix")) counts_all <- Matrix(counts_all, sparse = TRUE)
counts <- as(counts_all[, keep, drop = FALSE], "dgCMatrix")
metadata <- metadata[keep, , drop = FALSE]

if (!identical(colnames(counts), rownames(metadata))) {
  if (!all(colnames(counts) %in% rownames(metadata))) {
    stop("RNA counts and Seurat metadata do not contain the same cancer cells.", call. = FALSE)
  }
  metadata <- metadata[colnames(counts), , drop = FALSE]
}

niche_vec <- as.character(metadata$niche)
cell_type_vec <- as.character(metadata$cellType)
coordinate_mat <- as.matrix(metadata[, c("x_slide_mm", "y_slide_mm"), drop = FALSE]) * 1000
storage.mode(coordinate_mat) <- "double"
colnames(coordinate_mat) <- c("coordinate_x_um", "coordinate_y_um")

# Cell-type collapsing copied from the published SecAct preprocessing script.
cell_type_vec[grepl("Antibody.secreting.B.cells", cell_type_vec)] <- "B"
cell_type_vec[grepl("Mature.B.cells", cell_type_vec)] <- "B"
cell_type_vec[grepl("CD3+.alpha.beta.T.cells", cell_type_vec, fixed = TRUE)] <- "T.alpha.beta"
cell_type_vec[grepl("gamma.delta.T.cells.1", cell_type_vec)] <- "T.gamma.delta"
cell_type_vec[grepl("NK.like.cells", cell_type_vec)] <- "NK"
cell_type_vec[grepl("Hep", cell_type_vec)] <- "Hepatocyte"
cell_type_vec[grepl("Cholangiocytes", cell_type_vec)] <- "Cholangiocyte"
cell_type_vec[grepl("Erthyroid.cells", cell_type_vec)] <- "Erythrocyte"
cell_type_vec[grepl("Inflammatory.macrophages", cell_type_vec)] <- "Macrophage"
cell_type_vec[grepl("Non.inflammatory.macrophages", cell_type_vec)] <- "Macrophage"
cell_type_vec[grepl("Central.venous.LSECs", cell_type_vec)] <- "Endothelial"
cell_type_vec[grepl("Periportal.LSECs", cell_type_vec)] <- "Endothelial"
cell_type_vec[grepl("Portal.endothelial.cells", cell_type_vec)] <- "Endothelial"
cell_type_vec[grepl("Stellate.cells", cell_type_vec)] <- "Fibroblast"
cell_type_vec[niche_vec == "interface"] <- "Tumor_boundary"
cell_type_vec[grepl("tumor_1", cell_type_vec)] <- "Tumor_core"
cell_type_vec[grepl("tumor_2", cell_type_vec)] <- "Tumor_core"

spotCoordinates <- coordinate_mat
rownames(spotCoordinates) <- colnames(counts)
metaData <- data.frame(
  cellType = cell_type_vec,
  niche = niche_vec,
  row.names = colnames(counts),
  stringsAsFactors = FALSE
)

if (any(!is.finite(spotCoordinates))) stop("Cancer coordinates contain non-finite values.", call. = FALSE)
if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
  stop("Prepared counts contain duplicated gene or cell identifiers.", call. = FALSE)
}
if (!identical(colnames(counts), rownames(spotCoordinates)) ||
    !identical(colnames(counts), rownames(metaData))) {
  stop("Prepared counts, coordinates and metadata are not identically ordered.", call. = FALSE)
}

input_md5 <- unname(tools::md5sum(input_path))
rm(object, rna, counts_all, metadata, coordinate_mat)
gc()

message(sprintf(
  "Saving %s genes x %s cancer cells to %s...",
  format(nrow(counts), big.mark = ","),
  format(ncol(counts), big.mark = ","),
  output_path
))
save(counts, spotCoordinates, metaData, file = output_path, version = 3, compress = "gzip")

provenance_path <- paste0(output_path, ".provenance.txt")
writeLines(
  c(
    "source_provider=Bruker Spatial Biology / NanoString",
    "source_page=https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/human-liver-rna-ffpe-dataset/",
    "source_file_url=https://smi-public.objects.liquidweb.services/LiverDataReleaseSeurat_newUMAP.RDS",
    sprintf("source_file=%s", normalizePath(input_path)),
    sprintf("source_md5=%s", input_md5),
    sprintf("sample_column=%s", sample_column),
    sprintf("sample_value=%s", target),
    sprintf("genes=%d", nrow(counts)),
    sprintf("cells=%d", ncol(counts)),
    sprintf("nonzero_counts=%d", length(counts@x)),
    sprintf("cell_types_after_secact_mapping=%d", length(unique(metaData$cellType))),
    "coordinates=x_slide_mm,y_slide_mm converted from mm to um",
    "cell_type_mapping=published SecAct Secretome_s6_application_1_ST_2_cosmx_0_preprocess.R"
  ),
  provenance_path
)
write.table(
  data.frame(sample = names(sample_counts), cells = as.integer(sample_counts)),
  paste0(output_path, ".sample_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  as.data.frame(table(cellType = metaData$cellType, useNA = "ifany")),
  paste0(output_path, ".cell_type_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("Preparation completed.")

