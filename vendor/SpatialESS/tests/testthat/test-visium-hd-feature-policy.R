.hd_feature_policy_positions <- function(size = 8) {
  grid <- expand.grid(array_col = 0:1, array_row = 0:1)
  data.frame(barcode = sprintf("s_%03dum_%05d_%05d-1", size,
                               grid$array_row, grid$array_col),
              in_tissue = 1L, grid, stringsAsFactors = FALSE)
}

.hd_feature_policy_write <- function(square, counts, features, positions,
                                      kind = "raw_feature_bc_matrix",
                                      metadata_only = FALSE) {
  directory <- file.path(square, kind)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(features, file.path(directory, "features.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  if (!metadata_only) {
    Matrix::writeMM(counts, file.path(directory, "matrix.mtx"))
    writeLines(colnames(counts), file.path(directory, "barcodes.tsv"))
  }
  dir.create(file.path(square, "spatial"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(positions, file.path(square, "spatial", "tissue_positions.csv"),
                   row.names = FALSE)
}

.hd_feature_policy_fixture <- function() {
  positions <- .hd_feature_policy_positions()
  features <- data.frame(id = c("id_keep", "id_drop", "id_L"),
    name = c("HLA-DRA", "HLA-DRA", "L"), type = "Gene Expression")
  counts <- Matrix::Matrix(matrix(c(1:4, seq(10, 40, 10), 5:8),
    nrow = 3L, byrow = TRUE,
    dimnames = list(features$id, positions$barcode)), sparse = TRUE)
  list(counts = counts, features = features, positions = positions)
}

test_that("filtered feature IDs are selected before duplicate symbols are summed", {
  root <- tempfile("hd-probe-features-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  x <- .hd_feature_policy_fixture()
  square <- file.path(root, "square_008um")
  .hd_feature_policy_write(square, x$counts, x$features, x$positions)
  .hd_feature_policy_write(square, x$counts[c(1L, 3L), ],
    x$features[c(1L, 3L), ], x$positions,
    kind = "filtered_feature_bc_matrix", metadata_only = TRUE)

  filtered <- read_visium_hd(root, bin_size = 32)
  expect_identical(rownames(filtered$counts), c("HLA-DRA", "L"))
  expect_equal(as.numeric(filtered$counts[, 1L]), c(10, 26))
  expect_identical(filtered$provenance$feature_ids_kept, c("id_keep", "id_L"))
  expect_identical(filtered$provenance$feature_policy, "filtered")
  expect_identical(filtered$provenance$incomplete_count_source, FALSE)
  expect_identical(filtered$provenance$potential_incomplete_source, FALSE)
  expect_match(filtered$provenance$source_feature_id_source,
                "filtered_feature_bc_matrix$")
  expect_equal(filtered$provenance$feature_id_count, 2L)
  by_id <- read_visium_hd(root, 32, gene_column = "id")
  expect_identical(rownames(by_id$counts), c("id_keep", "id_L"))
  expect_equal(as.numeric(by_id$counts[, 1L]), c(10, 26))

  raw <- read_visium_hd(root, 32, feature_policy = "raw")
  expect_equal(as.numeric(raw$counts[, 1L]), c(110, 26))
  expect_identical(raw$provenance$feature_ids_kept, x$features$id)
  expect_identical(raw$provenance$feature_policy, "raw")
  native_raw <- read_visium_hd(root, 8, feature_policy = "raw")
  expect_equal(sum(native_raw$counts), sum(raw$counts))
  native_filtered <- read_visium_hd(root, 8)
  expect_equal(sum(native_filtered$counts), sum(filtered$counts))
})

test_that("missing filtered metadata never silently introduces raw probe features", {
  root <- tempfile("hd-missing-feature-policy-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  x <- .hd_feature_policy_fixture()
  square <- file.path(root, "square_008um")
  .hd_feature_policy_write(square, x$counts, x$features, x$positions)
  expect_error(read_visium_hd(root, 32), "filtered feature metadata")
  expect_error(read_visium_hd(root, 8), "filtered feature metadata|feature policy")
  expect_equal(sum(read_visium_hd(root, 32, feature_policy = "raw")$counts),
               sum(x$counts))
  expect_error(read_visium_hd(square, 32, source_bin_size = 8),
                "filtered feature metadata")
})

test_that("auto coarsening selects the largest usable raw divisor", {
  root <- tempfile("hd-raw-source-selection-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  x <- .hd_feature_policy_fixture()
  square8 <- file.path(root, "square_008um")
  .hd_feature_policy_write(square8, x$counts, x$features, x$positions)
  .hd_feature_policy_write(square8, x$counts[c(1L, 3L), ],
    x$features[c(1L, 3L), ], x$positions, kind = "filtered_feature_bc_matrix")
  pos16 <- data.frame(barcode = "s_016um_00000_00000-1", in_tissue = 1L,
                       array_row = 0L, array_col = 0L)
  filtered16 <- Matrix::Matrix(matrix(c(999, 888), ncol = 1L,
    dimnames = list(c("id_keep", "id_L"), pos16$barcode)), sparse = TRUE)
  square16 <- file.path(root, "square_016um")
  .hd_feature_policy_write(square16, filtered16,
    x$features[c(1L, 3L), ], pos16, kind = "filtered_feature_bc_matrix")

  coarse <- read_visium_hd(root, 32)
  expect_equal(coarse$source_bin_size, 8)
  expect_equal(as.numeric(coarse$counts[, 1L]), c(10, 26))
  expect_identical(coarse$provenance$matrix_kind, "raw_feature_bc_matrix")
  native <- read_visium_hd(root, 16)
  expect_equal(as.numeric(native$counts[, 1L]), c(999, 888))
  expect_equal(native$source_bin_size, 16)
  expect_error(read_visium_hd(root, 32, source_bin_size = 16), "Raw source counts")
  expect_warning(partial <- read_visium_hd(root, 32, source_bin_size = 16,
    allow_filtered_aggregation = TRUE), "incomplete")
  expect_equal(as.numeric(partial$counts[, 1L]), c(999, 888))
  expect_true(partial$provenance$incomplete_count_source)
  expect_true(partial$provenance$potential_incomplete_source)
  # Raw policy does not read filtered native16; it can build raw16 from raw8.
  native_raw <- read_visium_hd(root, 16, feature_policy = "raw")
  expect_equal(native_raw$source_bin_size, 8)
  expect_equal(as.numeric(native_raw$counts[, 1L]), c(110, 26))
})

test_that("all supplied positions are validated even without count columns", {
  x <- .hd_feature_policy_fixture()
  positions <- x$positions
  positions$array_row[4L] <- 7L
  expect_error(prepare_visium_hd(x$counts[, 1L, drop = FALSE], positions, 8),
                "grid indices encoded")
  positions <- x$positions
  positions$in_tissue <- NULL
  hd <- prepare_visium_hd(x$counts, positions, 8, in_tissue = FALSE)
  expect_false(hd$provenance$tissue_mask_available)
  expect_identical(hd$provenance$tissue_rule, "assumed_all_tissue_no_mask")
  coarse <- aggregate_visium_hd(hd, 32)
  expect_false(coarse$provenance$tissue_mask_available)
  expect_identical(coarse$provenance$tissue_rule, "assumed_all_tissue_no_mask")
})

test_that("direct aggregation records potentially truncated source count scope", {
  x <- .hd_feature_policy_fixture()
  positions <- x$positions
  positions$in_tissue[2:4] <- 0L
  hd <- prepare_visium_hd(x$counts, positions, 8, in_tissue = TRUE)
  expect_equal(hd$provenance$counts_removed_by_tissue_filter, 3L)
  coarse <- aggregate_visium_hd(hd, 32)
  expect_true(coarse$provenance$potential_incomplete_source)
  expect_true(coarse$provenance$aggregation[[1L]]$potential_incomplete_source)
  expect_equal(coarse$metadata$source_position_count, 4L)
  expect_equal(coarse$metadata$source_bin_count, 1L)
  expect_equal(as.numeric(coarse$counts[, 1L]), c(1, 10, 5))
  expect_identical(coarse$provenance$aggregation[[1L]]$native_output_equivalence,
                   "not empirically verified")
})
