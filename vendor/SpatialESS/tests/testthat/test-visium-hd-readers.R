# Real on-disk formats, rather than mocked reader return values. The input CSC
# layout intentionally contains zero entries, duplicate gene symbols and an
# antibody feature with much larger counts than any expression feature.
hd_reader_features <- function() {
  data.frame(
    id = c("ENSG_HLA_1", "ENSG_HLA_2", "ENSG_L", "ENSG_R", "ENSG_MT", "AB_CD3"),
    name = c("HLA-DRA", "HLA-DRA", "L", "R", "MT-CO1", "CD3 antibody"),
    feature_type = c(rep("Gene Expression", 5L), "Antibody Capture"),
    stringsAsFactors = FALSE
  )
}

hd_reader_counts <- function(barcodes) {
  value <- matrix(c(
    0, 1, 7, 2,
    3, 0, 11, 0,
    0, 5, 0, 1,
    8, 0, 13, 0,
    4, 1, 0, 0,
    99, 101, 103, 107
  ), nrow = 6L, byrow = TRUE,
  dimnames = list(hd_reader_features()$id, barcodes))
  methods::as(Matrix::Matrix(value, sparse = TRUE), "dgCMatrix")
}

hd_reader_write_h5 <- function(path, counts) {
  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all())
  mat <- h5$create_group("matrix")
  # 10x uses zero-based row indices and CSC column pointers.
  mat[["data"]] <- counts@x
  mat[["indices"]] <- counts@i
  mat[["indptr"]] <- counts@p
  mat[["shape"]] <- as.integer(dim(counts))
  mat[["barcodes"]] <- colnames(counts)
  features <- mat$create_group("features")
  info <- hd_reader_features()
  for (key in names(info)) features[[key]] <- info[[key]]
  invisible(path)
}

hd_reader_positions <- function(barcodes, rows, cols, tissue) {
  data.frame(
    barcode = barcodes, in_tissue = as.integer(tissue),
    array_row = as.integer(rows), array_col = as.integer(cols),
    # These image coordinates must not determine physical neighbourhoods.
    pxl_row_in_fullres = c(90123, 101, 42000, 18888),
    pxl_col_in_fullres = c(12, 81777, 302, 99999),
    stringsAsFactors = FALSE
  )
}

test_that("10x H5 CSC and shuffled Parquet positions preserve genes and grid units", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("arrow")
  root <- tempfile("visium-hd-h5-parquet-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  square <- file.path(root, "outs", "binned_outputs", "square_008um")
  spatial <- file.path(square, "spatial")
  dir.create(spatial, recursive = TRUE)
  ids <- c("s_008um_00001_00001-1", "s_008um_00000_00000-1",
           "s_008um_00002_00000-1", "s_008um_00000_00001-1")
  counts <- hd_reader_counts(ids)
  hd_reader_write_h5(file.path(square, "filtered_feature_bc_matrix.h5"), counts)
  positions <- hd_reader_positions(ids, c(1, 0, 2, 0), c(1, 0, 0, 1), c(1, 1, 0, 1))
  arrow::write_parquet(positions[c(3, 4, 1, 2), ],
                       file.path(spatial, "tissue_positions.parquet"))

  hd <- read_visium_hd(root, bin_size = 8, gene_column = "name")
  expected <- rbind(
    "HLA-DRA" = c(3, 1, 2), "L" = c(0, 5, 1),
    "R" = c(8, 0, 0), "MT-CO1" = c(4, 1, 0)
  )
  colnames(expected) <- ids[c(1, 2, 4)]
  expect_s4_class(hd$counts, "dgCMatrix")
  expect_equal(as.matrix(hd$counts), expected)
  expect_identical(rownames(hd$counts), c("HLA-DRA", "L", "R", "MT-CO1"))
  expect_identical(colnames(hd$counts), ids[c(1, 2, 4)])
  expect_equal(unname(hd$coordinates), cbind(c(12, 4, 12), c(12, 4, 4)))
  expect_identical(rownames(hd$coordinates), colnames(hd$counts))
  expect_identical(rownames(hd$metadata), colnames(hd$counts))
  expect_equal(hd$metadata$pxl_row_in_fullres, positions$pxl_row_in_fullres[c(1, 2, 4)])
  expect_equal(as.numeric(Matrix::colSums(hd$counts)), c(15, 7, 3))
  expect_identical(hd$provenance$matrix_kind, "filtered_feature_bc_matrix")
  expect_identical(hd$provenance$feature_policy, "filtered")
  expect_match(hd$provenance$positions_file, "tissue_positions[.]parquet$")

  by_id <- read_visium_hd(file.path(root, "outs"), bin_size = 8,
                           gene_column = "id", in_tissue = FALSE)
  expect_identical(rownames(by_id$counts), hd_reader_features()$id[1:5])
  expect_equal(as.matrix(by_id$counts), as.matrix(counts[1:5, , drop = FALSE]))
  expect_equal(sum(by_id$counts), sum(counts[1:5, ]))
  expect_identical(colnames(by_id$counts), ids)
  expect_equal(unname(by_id$coordinates), cbind(c(12, 4, 4, 12), c(12, 4, 20, 4)))
})

test_that("bin32 H5 fallback uses raw counts before Parquet tissue filtering", {
  skip_if_not_installed("hdf5r")
  skip_if_not_installed("arrow")
  root <- tempfile("visium-hd-raw-rebin-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  square <- file.path(root, "square_016um")
  spatial <- file.path(square, "spatial")
  dir.create(spatial, recursive = TRUE)
  ids <- c("s_016um_00000_00000-1", "s_016um_00000_00001-1",
           "s_016um_00001_00000-1", "s_016um_00002_00002-1")
  raw <- hd_reader_counts(ids)
  hd_reader_write_h5(file.path(square, "raw_feature_bc_matrix.h5"), raw)
  hd_reader_write_h5(file.path(square, "filtered_feature_bc_matrix.h5"),
                      raw[, 1L, drop = FALSE])
  positions <- hd_reader_positions(ids, c(0, 0, 1, 2), c(0, 1, 0, 2), c(1, 0, 0, 0))
  arrow::write_parquet(positions[c(4, 2, 1, 3), ],
                       file.path(spatial, "tissue_positions.parquet"))

  native <- read_visium_hd(root, bin_size = 16)
  expect_identical(colnames(native$counts), ids[1L])
  expect_equal(sum(native$counts), 15)
  expect_identical(native$provenance$matrix_kind, "filtered_feature_bc_matrix")

  binned <- read_visium_hd(root, bin_size = 32)
  expected <- matrix(c(22, 5, 21, 5), ncol = 1L,
                     dimnames = list(c("HLA-DRA", "L", "R", "MT-CO1"),
                                     "s_032um_00000_00000-1"))
  expect_equal(as.matrix(binned$counts), expected)
  expect_equal(sum(binned$counts), 53)
  expect_equal(unname(binned$coordinates), matrix(c(16, 16), nrow = 1L))
  expect_equal(binned$bin_size, 32)
  expect_equal(binned$source_bin_size, 16)
  expect_equal(binned$metadata$source_bin_count, 3L)
  expect_equal(binned$metadata$in_tissue, 1L)
  expect_identical(binned$provenance$matrix_kind, "raw_feature_bc_matrix")
  expect_match(binned$provenance$matrix_source, "raw_feature_bc_matrix[.]h5$")

  all_bins <- read_visium_hd(root, bin_size = 32, in_tissue = FALSE)
  expect_identical(colnames(all_bins$counts),
                   c("s_032um_00000_00000-1", "s_032um_00001_00001-1"))
  expect_equal(sum(all_bins$counts), sum(raw[1:5, ]))
  expect_equal(all_bins$metadata$in_tissue, c(1L, 0L))
  expect_equal(unname(all_bins$coordinates), rbind(c(16, 16), c(48, 48)))
  expect_identical(binned$provenance$feature_policy, "filtered")
  expect_false(binned$provenance$incomplete_count_source)
  # This fixture has the same feature IDs in raw and filtered files. The raw
  # policy therefore changes provenance, without changing its count values.
  raw_policy <- read_visium_hd(root, bin_size = 32, in_tissue = FALSE,
                               feature_policy = "raw")
  expect_equal(raw_policy$counts, all_bins$counts)
  expect_identical(raw_policy$provenance$feature_policy, "raw")
})

test_that("headerless six-column Space Ranger CSV works with plain MEX files", {
  root <- tempfile("visium-hd-headerless-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  square <- file.path(root, "square_008um")
  spatial <- file.path(square, "spatial")
  mex <- file.path(square, "filtered_feature_bc_matrix")
  dir.create(spatial, recursive = TRUE)
  dir.create(mex)
  ids <- c("s_008um_00001_00001-1", "s_008um_00000_00000-1",
           "s_008um_00002_00000-1", "s_008um_00000_00001-1")
  counts <- hd_reader_counts(ids)
  Matrix::writeMM(counts, file.path(mex, "matrix.mtx"))
  utils::write.table(hd_reader_features(), file.path(mex, "features.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(ids, file.path(mex, "barcodes.tsv"))
  positions <- hd_reader_positions(ids, c(1, 0, 2, 0), c(1, 0, 0, 1), c(1, 1, 0, 1))
  utils::write.table(positions[c(4, 3, 2, 1), ],
                     file.path(spatial, "tissue_positions_list.csv"),
                     sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)

  hd <- read_visium_hd(root, bin_size = 8)
  expect_identical(colnames(hd$counts), ids[c(1, 2, 4)])
  expect_identical(rownames(hd$counts), c("HLA-DRA", "L", "R", "MT-CO1"))
  expect_equal(as.numeric(hd$counts["HLA-DRA", ]), c(3, 1, 2))
  expect_equal(as.numeric(Matrix::colSums(hd$counts)), c(15, 7, 3))
  expect_equal(unname(hd$coordinates), cbind(c(12, 4, 12), c(12, 4, 4)))
  expect_match(hd$provenance$positions_file, "tissue_positions_list[.]csv$")
})
