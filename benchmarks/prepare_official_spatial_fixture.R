#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
out_dir <- file.path(project_dir, "benchmarks", "results", "official_spatial_fixture")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
})

set.seed(20260731)
coordinates_pixel <- as.matrix(expand.grid(x = 0:29, y = 0:19))
rownames(coordinates_pixel) <- paste0("cell", seq_len(nrow(coordinates_pixel)))
coordinates_um <- coordinates_pixel * 5
group <- factor(
  c("A", "B", "C")[pmin(3L, coordinates_pixel[, "x"] %/% 10L + 1L)],
  levels = c("A", "B", "C")
)
meta <- data.frame(labels = group, row.names = rownames(coordinates_pixel))

data(CellChatDB.human, package = "CellChat")
db <- CellChatDB.human
interaction <- db$interaction
empty_optional <- function(column) {
  !column %in% colnames(interaction) |
    is.na(interaction[[column]]) | !nzchar(as.character(interaction[[column]]))
}
simple <- !interaction$ligand %in% rownames(db$complex) &
  !interaction$receptor %in% rownames(db$complex) &
  empty_optional("agonist") & empty_optional("antagonist") &
  empty_optional("co_A_receptor") & empty_optional("co_I_receptor")
candidate <- interaction[simple, , drop = FALSE]
candidate <- candidate[!duplicated(paste(candidate$ligand, candidate$receptor)), , drop = FALSE]
annotation_order <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor")
selected <- unlist(lapply(annotation_order, function(annotation) {
  which(candidate$annotation == annotation)[seq_len(min(4L, sum(candidate$annotation == annotation)))]
}), use.names = FALSE)
lr <- candidate[selected, , drop = FALSE]
stopifnot(nrow(lr) == 12L)

genes <- unique(c(as.character(lr$ligand), as.character(lr$receptor)))
counts <- matrix(
  rpois(length(genes) * nrow(coordinates_pixel), lambda = 3),
  nrow = length(genes),
  dimnames = list(genes, rownames(coordinates_pixel))
)
for (index in seq_len(nrow(lr))) {
  sender <- levels(group)[(index - 1L) %% 3L + 1L]
  receiver <- levels(group)[index %% 3L + 1L]
  ligand <- as.character(lr$ligand[index])
  receptor <- as.character(lr$receptor[index])
  counts[ligand, group == sender] <- counts[ligand, group == sender] +
    rpois(sum(group == sender), lambda = 8)
  counts[receptor, group == receiver] <- counts[receptor, group == receiver] +
    rpois(sum(group == receiver), lambda = 8)
}
counts[matrix(runif(length(counts)) < 0.1, nrow = nrow(counts))] <- 0
counts <- Matrix(counts, sparse = TRUE)
library_size <- Matrix::colSums(counts)
expression <- counts %*% Diagonal(x = 1000 / library_size)
dimnames(expression) <- dimnames(counts)
expression@x <- log1p(expression@x)
expression <- drop0(as(expression, "dgCMatrix"))

fixture <- list(
  expression = expression,
  meta = meta,
  group = group,
  coordinates_pixel = coordinates_pixel,
  coordinates_um = coordinates_um,
  scale_factors = list(spot.diameter = 5, spot = 1),
  lr = lr,
  seed = 20260731L,
  interaction_length_um = 20,
  scale_distance = 0.05,
  k_min = 5L,
  nboot = 20L
)
saveRDS(fixture, file.path(out_dir, "fixture.rds"), compress = FALSE)
write.table(
  lr[, c("interaction_name", "pathway_name", "ligand", "receptor", "annotation")],
  file.path(out_dir, "lr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(cells = ncol(expression), genes = nrow(expression),
             groups = nlevels(group), lr = nrow(lr)),
  file.path(out_dir, "fixture_summary.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE
)

