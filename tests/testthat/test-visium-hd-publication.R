# Publication-facing binning regressions: tissue positions and observed count
# columns are distinct universes. A raw matrix can omit zero-read barcodes.
.hd_publication_positions <- function(size = 8, n = 2L) {
  grid <- expand.grid(array_col = seq_len(n) - 1L,
                      array_row = seq_len(n) - 1L)
  data.frame(
    barcode = sprintf("s_%03dum_%05d_%05d-1", size,
                       grid$array_row, grid$array_col),
    in_tissue = integer(nrow(grid)), grid, stringsAsFactors = FALSE
  )
}

.hd_publication_mex <- function(square, counts, positions,
                                kind = "raw_feature_bc_matrix") {
  mex <- file.path(square, kind)
  spatial <- file.path(square, "spatial")
  dir.create(mex, recursive = TRUE, showWarnings = FALSE)
  dir.create(spatial, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(counts, file.path(mex, "matrix.mtx"))
  utils::write.table(
    data.frame(paste0("id", seq_len(nrow(counts))), rownames(counts),
               "Gene Expression"),
    file.path(mex, "features.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  writeLines(colnames(counts), file.path(mex, "barcodes.tsv"))
  utils::write.csv(positions, file.path(spatial, "tissue_positions.csv"),
                   row.names = FALSE)
}

test_that("zero-read tissue positions affect coarser tissue calls", {
  root <- tempfile("hd-zero-read-tissue-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  positions <- .hd_publication_positions()
  positions$in_tissue[1L] <- 1L
  # The only observed count column is outside tissue. An omitted zero-read
  # source barcode lies in tissue in the same target bin.
  counts <- Matrix::Matrix(matrix(c(3, 7), ncol = 1L,
    dimnames = list(c("L", "R"), positions$barcode[4L])), sparse = TRUE)
  square <- file.path(root, "square_008um")
  .hd_publication_mex(square, counts, positions)
  hd <- read_visium_hd(root, bin_size = 16, feature_policy = "raw")
  expect_identical(colnames(hd$counts), "s_016um_00000_00000-1")
  expect_equal(as.numeric(hd$counts[, 1L]), c(3, 7))
  expect_identical(hd$metadata$in_tissue, 1L)
  expect_equal(unname(hd$coordinates), matrix(c(8, 8), nrow = 1L))
})

test_that("filtered-only coarsening requires explicit partial-count opt-in", {
  root <- tempfile("hd-filtered-coarsening-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  positions <- .hd_publication_positions()
  positions$in_tissue[1L] <- 1L
  counts <- Matrix::Matrix(matrix(c(2, 5), ncol = 1L,
    dimnames = list(c("L", "R"), positions$barcode[1L])), sparse = TRUE)
  square <- file.path(root, "square_008um")
  .hd_publication_mex(square, counts, positions,
                       kind = "filtered_feature_bc_matrix")
  # Native filtered matrices remain valid inputs. The explicit policy concerns
  # reconstructing a different resolution from a truncated barcode universe.
  expect_equal(sum(read_visium_hd(root, bin_size = 8)$counts), 7)
  expect_error(read_visium_hd(root, bin_size = 16), "[Rr]aw|filtered")
  expect_warning(
    partial <- read_visium_hd(root, bin_size = 16,
                              allow_filtered_aggregation = TRUE),
    "[Ff]iltered|partial|incomplete"
  )
  expect_equal(as.numeric(partial$counts[, 1L]), c(2, 5))
  expect_identical(partial$provenance$matrix_kind,
                   "filtered_feature_bc_matrix")
  expect_true(partial$provenance$incomplete_count_source)
})

test_that("nested raw-count summation preserves genes, counts and tissue calls", {
  positions <- .hd_publication_positions(size = 2, n = 8L)
  positions$in_tissue[1L] <- 1L
  # The observed source bin and tissue source bin occupy different 8 um parents
  # but one 16 um parent. This checks preservation through successive binning.
  counts <- Matrix::Matrix(matrix(c(11, 0, 19), ncol = 1L,
    dimnames = list(c("HLA-DRA", "zero_gene", "R"),
                    positions$barcode[nrow(positions)])), sparse = TRUE)
  hd2 <- prepare_visium_hd(counts, positions, bin_size = 2,
                            in_tissue = FALSE)
  direct <- aggregate_visium_hd(hd2, 16)
  via8 <- aggregate_visium_hd(aggregate_visium_hd(hd2, 8), 16)
  via4_8 <- aggregate_visium_hd(
    aggregate_visium_hd(aggregate_visium_hd(hd2, 4), 8), 16
  )
  expect_equal(direct$counts, via8$counts)
  expect_equal(direct$counts, via4_8$counts)
  expect_identical(rownames(direct$counts), rownames(counts))
  expect_equal(as.numeric(Matrix::rowSums(direct$counts)),
               as.numeric(Matrix::rowSums(counts)))
  expect_identical(direct$metadata$in_tissue, 1L)
  expect_identical(via8$metadata$in_tissue, 1L)
  expect_identical(via4_8$metadata$in_tissue, 1L)
})

test_that("missing tissue annotations cannot silently become a tissue mask", {
  positions <- .hd_publication_positions()
  positions$in_tissue <- NULL
  counts <- Matrix::Matrix(matrix(1, nrow = 2L, ncol = nrow(positions),
    dimnames = list(c("L", "R"), positions$barcode)), sparse = TRUE)
  expect_error(prepare_visium_hd(counts, positions, 8, in_tissue = TRUE),
                "in_tissue|tissue")
})
