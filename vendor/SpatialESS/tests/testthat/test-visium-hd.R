.hd_fixture <- function(size = 8, n = 4L) {
  grid <- expand.grid(array_col = seq_len(n) - 1L, array_row = seq_len(n) - 1L)
  ids <- sprintf("s_%03dum_%05d_%05d-1", size, grid$array_row, grid$array_col)
  counts <- Matrix::Matrix(matrix(
    seq_len(2L * nrow(grid)), nrow = 2L,
    dimnames = list(c("L", "R"), ids)
  ), sparse = TRUE)
  positions <- data.frame(barcode = ids, in_tissue = 1L, grid,
                          pxl_row_in_fullres = 900 + 5 * grid$array_row,
                          pxl_col_in_fullres = 400 + 5 * grid$array_col)
  list(counts = counts, positions = positions)
}

.hd_write_mex <- function(directory, counts, positions, gzip = FALSE,
                           matrix_kind = "filtered_feature_bc_matrix") {
  matrix_dir <- file.path(directory, matrix_kind)
  dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(directory, "spatial"), recursive = TRUE, showWarnings = FALSE)
  matrix_path <- file.path(matrix_dir, "matrix.mtx")
  Matrix::writeMM(counts, matrix_path)
  utils::write.table(data.frame(paste0("id", seq_len(nrow(counts))), rownames(counts), "Gene Expression"),
                     file.path(matrix_dir, "features.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  writeLines(colnames(counts), file.path(matrix_dir, "barcodes.tsv"))
  utils::write.csv(positions, file.path(directory, "spatial", "tissue_positions.csv"), row.names = FALSE)
  if (gzip) {
    for (file in list.files(matrix_dir, full.names = TRUE)) {
      bytes <- readBin(file, "raw", n = file.info(file)$size)
      connection <- gzfile(paste0(file, ".gz"), "wb")
      writeBin(bytes, connection)
      close(connection)
      unlink(file)
    }
  }
  invisible(directory)
}

test_that("Visium HD preparation aligns by barcode and uses grid micrometres", {
  fixture <- .hd_fixture()
  set.seed(3)
  positions <- fixture$positions[sample(nrow(fixture$positions)), ]
  hd <- prepare_visium_hd(fixture$counts, positions, 8)
  expect_s3_class(hd, "SpatialESSVisiumHD")
  expect_s4_class(hd$counts, "dgCMatrix")
  expect_identical(rownames(hd$metadata), colnames(fixture$counts))
  expect_identical(rownames(hd$coordinates), colnames(fixture$counts))
  expect_equal(unname(hd$coordinates[1:4, ]), cbind(c(4, 12, 20, 28), rep(4, 4)))
  expect_equal(hd$counts, methods::as(fixture$counts, "dgCMatrix"))
})

test_that("Sparse bin aggregation conserves counts and uses true target centres", {
  fixture <- .hd_fixture()
  hd8 <- prepare_visium_hd(fixture$counts, fixture$positions, 8)
  hd16 <- aggregate_visium_hd(hd8, 16)
  hd32 <- aggregate_visium_hd(hd8, 32)
  expect_equal(dim(hd16$counts), c(2L, 4L))
  expect_equal(as.numeric(Matrix::rowSums(hd8$counts)), as.numeric(Matrix::rowSums(hd16$counts)))
  expect_equal(as.numeric(Matrix::rowSums(hd8$counts)), as.numeric(Matrix::rowSums(hd32$counts)))
  expect_equal(hd32$counts, aggregate_visium_hd(hd16, 32)$counts)
  expect_equal(unname(hd32$coordinates), matrix(c(16, 16), nrow = 1))
  expect_identical(colnames(hd32$counts), "s_032um_00000_00000-1")
  expect_equal(hd16$counts[, 1], Matrix::rowSums(fixture$counts[, c(1, 2, 5, 6)]))
  expect_equal(hd16$metadata$source_bin_count, rep(4L, 4))
  # A partial target at a boundary still has the physical centre of its grid bin.
  one <- prepare_visium_hd(fixture$counts[, 1, drop = FALSE], fixture$positions[1, ], 8)
  expect_equal(unname(aggregate_visium_hd(one, 32)$coordinates), matrix(c(16, 16), nrow = 1))
  expect_error(aggregate_visium_hd(hd16, 8), "finer bins")
  expect_error(aggregate_visium_hd(hd16, 24), "integer multiple")
})

test_that("Aggregation is order independent and keeps barcode suffixes separate", {
  fixture <- .hd_fixture(n = 2)
  counts <- cbind(fixture$counts, fixture$counts)
  colnames(counts)[5:8] <- sub("-1$", "-2", colnames(counts)[5:8])
  positions <- rbind(fixture$positions, fixture$positions)
  positions$barcode <- colnames(counts)
  hd <- prepare_visium_hd(counts, positions, 8)
  mixed <- prepare_visium_hd(counts[, 8:1], positions, 8)
  original <- aggregate_visium_hd(hd, 16)
  shuffled <- aggregate_visium_hd(mixed, 16)
  expect_equal(original$counts, shuffled$counts)
  expect_identical(colnames(original$counts), c("s_016um_00000_00000-1", "s_016um_00000_00000-2"))
  expect_equal(original$counts[, 1], original$counts[, 2])
})

test_that("Native resolution is preferred and missing bin32 is aggregated", {
  root <- tempfile("visium_hd_")
  on.exit(unlink(root, recursive = TRUE))
  fixture <- .hd_fixture()
  dir8 <- file.path(root, "outs", "binned_outputs", "square_008um")
  .hd_write_mex(dir8, fixture$counts, fixture$positions, gzip = TRUE)
  .hd_write_mex(dir8, fixture$counts, fixture$positions, gzip = TRUE,
                  matrix_kind = "raw_feature_bc_matrix")
  hd8 <- read_visium_hd(root, bin_size = 8)
  hd16 <- aggregate_visium_hd(hd8, 16)
  dir16 <- file.path(root, "outs", "binned_outputs", "square_016um")
  # Make native bin16 deliberately distinct to show it is read rather than rebuilt.
  .hd_write_mex(dir16, hd16$counts * 2, hd16$metadata)
  .hd_write_mex(dir16, hd16$counts * 2, hd16$metadata,
                  matrix_kind = "raw_feature_bc_matrix")
  native <- read_visium_hd(root, bin_size = 16)
  expect_equal(native$counts, hd16$counts * 2)
  expect_equal(native$source_bin_size, 16)
  hd32 <- read_visium_hd(root, bin_size = 32)
  expect_equal(hd32$counts, aggregate_visium_hd(native, 32)$counts)
  expect_equal(hd32$source_bin_size, 16)
  pinned <- read_visium_hd(root, bin_size = 32, source_bin_size = 8)
  expect_equal(pinned$counts, aggregate_visium_hd(hd8, 32)$counts)
  expect_equal(pinned$source_bin_size, 8)
  expect_error(read_visium_hd(dir8, bin_size = 32), "source_bin_size explicitly")
  expect_equal(read_visium_hd(dir8, bin_size = 32, source_bin_size = 8)$counts, pinned$counts)
  expect_error(read_visium_hd(root, bin_size = 6), "divides")
})

test_that("Raw source bins are summed before target tissue filtering", {
  root <- tempfile("visium_hd_")
  on.exit(unlink(root, recursive = TRUE))
  fixture <- .hd_fixture(n = 2)
  fixture$positions$in_tissue <- c(1L, 0L, 0L, 0L)
  dir8 <- file.path(root, "binned_outputs", "square_008um")
  .hd_write_mex(dir8, fixture$counts[, 1, drop = FALSE], fixture$positions)
  .hd_write_mex(dir8, fixture$counts, fixture$positions, matrix_kind = "raw_feature_bc_matrix")
  native <- read_visium_hd(root, bin_size = 8)
  expect_equal(ncol(native$counts), 1)
  hd16 <- read_visium_hd(root, bin_size = 16)
  expect_equal(as.numeric(hd16$counts[, 1]), as.numeric(Matrix::rowSums(fixture$counts)))
  expect_equal(hd16$metadata$in_tissue, 1L)
  expect_identical(hd16$provenance$matrix_kind, "raw_feature_bc_matrix")
  expect_equal(hd16$metadata$source_bin_count, 4L)
})

test_that("Gene symbols are preserved and duplicates summed", {
  fixture <- .hd_fixture(n = 2)
  rownames(fixture$counts) <- c("HLA-DRA", "HLA-DRA")
  hd <- prepare_visium_hd(fixture$counts, fixture$positions, 8)
  expect_identical(rownames(hd$counts), "HLA-DRA")
  expect_equal(as.numeric(hd$counts[1, ]), as.numeric(Matrix::colSums(fixture$counts)))
})

test_that("Bad counts, identifiers and positions fail before fitting", {
  fixture <- .hd_fixture(n = 2)
  expect_error(prepare_visium_hd(fixture$counts / 2, fixture$positions, 8), "integer raw counts")
  negative <- fixture$counts
  negative[1, 1] <- -1
  expect_error(prepare_visium_hd(negative, fixture$positions, 8), "non-negative integer")
  expect_error(prepare_visium_hd(fixture$counts, fixture$positions[-1, ], 8), "requires a spatial position")
  duplicate <- fixture$positions
  duplicate$barcode[2] <- duplicate$barcode[1]
  expect_error(prepare_visium_hd(fixture$counts, duplicate, 8), "unique")
  duplicate <- fixture$positions
  duplicate$array_row[2] <- duplicate$array_row[1]
  duplicate$array_col[2] <- duplicate$array_col[1]
  expect_error(prepare_visium_hd(fixture$counts, duplicate, 8), "duplicate array-grid")
  expect_error(prepare_visium_hd(fixture$counts, fixture$positions, 16), "resolution encoded")
  expect_error(prepare_visium_hd(fixture$counts, fixture$positions, 3), "positive even integer")
  fixture$positions$in_tissue <- 0L
  expect_error(prepare_visium_hd(fixture$counts, fixture$positions, 8), "No bins remain")
})
