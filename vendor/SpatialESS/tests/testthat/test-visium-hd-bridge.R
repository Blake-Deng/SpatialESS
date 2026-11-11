hd_bridge_fixture <- function() {
  ids <- paste0("b", seq_len(6))
  counts <- Matrix::Matrix(matrix(
    c(10, 5, 85, 6, 9, 85, 2, 8, 90,
      8, 2, 90, 3, 7, 90, 9, 1, 90),
    nrow = 3L, dimnames = list(c("L", "R", "housekeeping"), ids)
  ), sparse = TRUE)
  positions <- data.frame(barcode = ids, in_tissue = 1L,
                          array_row = 0L, array_col = 0:5)
  hd <- prepare_visium_hd(counts, positions, bin_size = 8)
  weights <- Matrix::Matrix(matrix(
    c(.8, .2, .2, .8, .6, .4, 0, 0, .5, .5, .9, .1),
    ncol = 2, byrow = TRUE, dimnames = list(ids, c("A", "B"))
  ), sparse = TRUE)
  annotations <- data.frame(
    barcode = ids, status = c("ok", "ok", "ok", "filtered", "ambiguous", "reject"),
    dominant_type = c("A", "B", "A", NA, NA, NA),
    max_weight = c(.8, .8, .6, 0, .5, .9), row.names = ids
  )
  list(hd = hd, deconvolution = list(weights = weights, annotations = annotations))
}

test_that("HD bridge normalizes by all genes and audits every bin", {
  f <- hd_bridge_fixture()
  x <- prepare_visium_hd_spatialess(f$hd, f$deconvolution, min_weight = .7)
  expect_identical(colnames(x$expression), c("b1", "b2"))
  expect_equal(as.numeric(x$expression["L", ]), log1p(c(1000, 600)))
  expect_identical(rownames(x$coordinates), colnames(x$expression))
  expect_identical(names(x$group), colnames(x$expression))
  expect_identical(as.character(x$group), c("A", "B"))
  expect_equal(nrow(x$audit), 6L)
  expect_identical(x$audit$exclusion_reason,
                   c("retained", "retained", "below_min_weight",
                     rep("rctd_not_accepted", 3)))
  expect_equal(f$hd$counts["L", "b1"], 10)
  raw <- prepare_visium_hd_spatialess(f$hd, f$deconvolution,
                                     normalization = "none")
  expect_equal(raw$expression, f$hd$counts[, 1:3])
})

test_that("HD bridge realigns keyed deconvolution and rejects stale labels", {
  f <- hd_bridge_fixture()
  original <- prepare_visium_hd_spatialess(f$hd, f$deconvolution)
  f$deconvolution$annotations <- f$deconvolution$annotations[6:1, ]
  f$deconvolution$weights <- f$deconvolution$weights[c(2:6, 1), ]
  expect_equal(prepare_visium_hd_spatialess(f$hd, f$deconvolution), original)
  f$deconvolution$annotations$max_weight[6] <- .9
  expect_error(prepare_visium_hd_spatialess(f$hd, f$deconvolution), "disagrees")
  f <- hd_bridge_fixture()
  f$deconvolution$annotations$barcode[2] <- "b1"
  expect_error(prepare_visium_hd_spatialess(f$hd, f$deconvolution), "exactly once")
  f <- hd_bridge_fixture()
  f$deconvolution$weights[1, ] <- c(.8, .8)
  expect_error(prepare_visium_hd_spatialess(f$hd, f$deconvolution), "sum to one")
  f <- hd_bridge_fixture()
  expect_error(prepare_visium_hd_spatialess(f$hd, f$deconvolution, 1.1), "min_weight")
  expect_error(prepare_visium_hd_spatialess(f$hd, f$deconvolution, .99), "No accepted")
})

test_that("HD communication agrees with direct bin-level SpatialESS", {
  f <- hd_bridge_fixture()
  lr <- data.frame(ligand = "L", receptor = "R", interaction_name = "L_R",
                   annotation = "Secreted Signaling")
  complex <- data.frame(subunit_1 = character())
  cofactor <- data.frame(cofactor1 = character())
  opts <- list(tol = 0, interaction.range = 40, scale.distance = .1,
               min.percent = 0, min.cells.sr = 1L, nboot = 2L, seed.use = 7L)
  prepared <- prepare_visium_hd_spatialess(f$hd, f$deconvolution)
  direct <- do.call(spatialess, c(list(
    expression = prepared$expression, coordinates = prepared$coordinates,
    group = prepared$group, lr = lr, complex = complex, cofactor = cofactor
  ), opts))
  bridged <- do.call(spatialess_from_visium_hd, c(list(
    hd = f$hd, deconvolution = f$deconvolution,
    lr = lr, complex = complex, cofactor = cofactor
  ), opts))
  expect_equal(bridged$records, direct$records)
  expect_identical(bridged$visium_hd$weights, f$deconvolution$weights)
  expect_equal(nrow(bridged$visium_hd$audit), 6L)
  expect_error(spatialess_from_visium_hd(
    f$hd, f$deconvolution, lr, complex, cofactor, ratio = .5
  ), "ratio must equal 1")
})
