# Mock the upstream boundary, exercising alignment, proportions and global-fit
# batching without making spacexr a mandatory package dependency.
.test_hd_rctd_fixture <- function() {
  genes <- paste0("g", seq_len(12L))
  bins <- paste0("b", seq_len(5L))
  counts <- Matrix::Matrix(matrix(10, 12L, 5L, dimnames = list(genes, bins)), sparse = TRUE)
  counts[, 5L] <- 0
  reference <- Matrix::Matrix(matrix(20, 12L, 4L,
    dimnames = list(genes, paste0("r", seq_len(4L)))), sparse = TRUE)
  coordinates <- matrix(seq_len(10L), 5L, 2L, dimnames = list(bins, c("x", "y")))
  hd <- structure(list(counts = counts, coordinates = coordinates,
    metadata = data.frame(barcode = bins), bin_size = 32, source_bin_size = 8,
    provenance = list(source = "synthetic")), class = "SpatialESSVisiumHD")
  list(hd = hd, reference = reference,
       types = setNames(c("A", "A", "B", "B"), colnames(reference)))
}

.test_hd_rctd_backend <- function(drop_bin = "b3") {
  if (!methods::isClass("SpatialESSTestReference")) {
    methods::setClass("SpatialESSTestReference",
      slots = c(counts = "dgCMatrix", cell_types = "factor", nUMI = "numeric"))
  }
  if (!methods::isClass("SpatialESSTestSpatialRNA")) {
    methods::setClass("SpatialESSTestSpatialRNA",
      slots = c(counts = "dgCMatrix", coords = "data.frame", nUMI = "numeric"))
  }
  if (!methods::isClass("SpatialESSTestRCTD")) {
    methods::setClass("SpatialESSTestRCTD", slots = c(
      spatialRNA = "SpatialESSTestSpatialRNA", originalSpatialRNA = "SpatialESSTestSpatialRNA",
      config = "list", results = "list", cell_type_info = "list", internal_vars = "list",
      reference = "SpatialESSTestReference"))
  }
  calls <- new.env(parent = emptyenv())
  calls$create <- calls$bulk <- calls$sigma <- calls$pixels <- 0L
  backend <- list(
    Reference = function(counts, cell_types, nUMI) {
      calls$reference_nUMI <- nUMI
      calls$reference_types <- cell_types
      keep <- nUMI >= 100
      methods::new("SpatialESSTestReference",
                   counts = methods::as(counts[, keep, drop = FALSE], "dgCMatrix"),
                   cell_types = cell_types[keep], nUMI = nUMI[keep])
    },
    SpatialRNA = function(coords, counts, nUMI) {
      calls$spatial_nUMI <- nUMI
      methods::new("SpatialESSTestSpatialRNA", counts = methods::as(counts, "dgCMatrix"),
                   coords = coords, nUMI = nUMI)
    },
    create.RCTD = function(spatialRNA, reference, max_cores, ...) {
      calls$create <- calls$create + 1L
      calls$create_args <- list(...)
      # Mimic upstream's plotting-reference subsampling, which must not be
      # reported as the number of cells used to construct cell-type profiles.
      reduced <- reference
      keep <- match(levels(reference@cell_types), as.character(reference@cell_types))
      keep <- keep[!is.na(keep)]
      reduced@counts <- reduced@counts[, keep, drop = FALSE]
      reduced@cell_types <- reduced@cell_types[keep]
      reduced@nUMI <- reduced@nUMI[keep]
      methods::new("SpatialESSTestRCTD", spatialRNA = spatialRNA,
                   originalSpatialRNA = spatialRNA,
                   config = list(UMI_min = calls$create_args$UMI_min,
                                 counts_MIN = 10, max_cores = max_cores),
                   results = list(), reference = reduced,
                   cell_type_info = list(renorm = list(NULL, c("A", "B"), 2L)),
                   internal_vars = list(cell_types_assigned = FALSE))
    },
    fitBulk = function(model) {
      calls$bulk <- calls$bulk + 1L
      model@config$global_bulk <- TRUE
      model
    },
    choose_sigma_c = function(model) {
      calls$sigma <- calls$sigma + 1L
      expect_true(model@config$global_bulk)
      ids <- setdiff(colnames(model@spatialRNA@counts), drop_bin)
      model@spatialRNA <- .hd_rctd_subset_spatial(model@spatialRNA, ids)
      model@config$global_sigma <- stats::runif(1)
      model
    },
    fitPixels = function(model, doublet_mode) {
      calls$pixels <- calls$pixels + 1L
      expect_true(model@config$global_bulk)
      expect_true(is.finite(model@config$global_sigma))
      ids <- colnames(model@spatialRNA@counts)
      w <- matrix(rep(c(1, 3), each = length(ids)), nrow = length(ids),
                  dimnames = list(ids, c("A", "B")))
      if ("b2" %in% ids) w["b2", ] <- 2
      if ("b4" %in% ids) w["b4", ] <- 0
      if (doublet_mode == "doublet") {
        spot_class <- ifelse(ids == "b4", "reject",
                             ifelse(ids == "b2", "doublet_uncertain", "singlet"))
        results_df <- data.frame(
          spot_class = factor(spot_class, levels = c(
            "reject", "singlet", "doublet_certain", "doublet_uncertain")),
          first_type = factor(rep("A", length(ids)), levels = c("A", "B")),
          second_type = factor(rep("B", length(ids)), levels = c("A", "B")),
          min_score = as.numeric(match(ids, paste0("b", 1:5))),
          singlet_score = rep(7, length(ids)),
          conv_all = ids != "b4", conv_doublet = ids != "b2",
          row.names = ids
        )
        first_weight <- unname(c(b1 = .3, b2 = .6, b3 = .4, b4 = .35)[ids])
        weights_doublet <- cbind(first_type = first_weight, second_type = 1 - first_weight)
        rownames(weights_doublet) <- ids
        model@results <- list(weights = w, results_df = results_df,
                              weights_doublet = weights_doublet)
      } else if (doublet_mode == "multi") {
        model@results <- stats::setNames(lapply(ids, function(id) {
          list(all_weights = c(A = 100, B = 200),
               cell_type_list = c("A", "B"), sub_weights = c(1, 3),
               conf_list = c(TRUE, id != "b2"),
               min_score = match(id, paste0("b", 1:5)),
               conv_all = TRUE, conv_sub = id != "b2")
        }), ids)
      } else model@results <- list(weights = w)
      model
    },
    version = "test-original-api"
  )
  list(backend = backend, calls = calls)
}

test_that("Visium HD RCTD validates raw counts and named alignment before backend loading", {
  x <- .test_hd_rctd_fixture()
  expect_error(run_visium_hd_rctd(list(), x$reference, x$types), "SpatialESSVisiumHD")
  bad <- x$hd
  bad$counts@x[1L] <- -1
  expect_error(run_visium_hd_rctd(bad, x$reference, x$types), "raw integer counts")
  bad$counts@x[1L] <- 0.2
  expect_error(run_visium_hd_rctd(bad, x$reference, x$types), "raw integer counts")
  bad$counts@x[1L] <- Inf
  expect_error(run_visium_hd_rctd(bad, x$reference, x$types), "raw integer counts")
  expect_error(run_visium_hd_rctd(x$hd, x$reference, unname(x$types)), "uniquely named")
  bad_types <- x$types
  bad_types[1L] <- NA
  expect_error(run_visium_hd_rctd(x$hd, x$reference, bad_types), "nonmissing")
  bad_types[1L] <- "A/B"
  expect_error(run_visium_hd_rctd(x$hd, x$reference, bad_types), "without '/'")
  bad <- x$hd
  bad$coordinates <- bad$coordinates[5:1, ]
  expect_error(run_visium_hd_rctd(bad, x$reference, x$types), "coordinates")
  expect_error(run_visium_hd_rctd(x$hd, x$reference, x$types, min_umi = 1000), "No bins")
  expect_error(run_visium_hd_rctd(x$hd, x$reference[1:5, ], x$types), "10 shared genes")
  expect_error(run_visium_hd_rctd(x$hd, x$reference, x$types, batch_size = 0), "batch_size")
  expect_error(run_visium_hd_rctd(x$hd, x$reference, x$types, create_args = list(reference = 1)), "cannot override")
  expect_error(run_visium_hd_rctd(x$hd, x$reference, x$types, run_args = list(doublet_mode = "doublet")), "only supports")
})

test_that("RCTD preserves all-bin alignment and performs one global fit before batching", {
  x <- .test_hd_rctd_fixture()
  mock <- .test_hd_rctd_backend()
  local_mocked_bindings(.hd_rctd_backend = function() mock$backend)
  set.seed(777)
  rng <- .Random.seed
  fit <- run_visium_hd_rctd(x$hd, x$reference, x$types[c(4, 2, 1, 3)],
                            batch_size = 1L, keep_fit = TRUE)
  expect_identical(.Random.seed, rng)
  expect_s3_class(fit, "SpatialESSVisiumHDRCTD")
  expect_s4_class(fit$weights, "dgCMatrix")
  expect_identical(rownames(fit$weights), colnames(x$hd$counts))
  expect_identical(fit$annotations$barcode, colnames(x$hd$counts))
  expect_equal(mock$calls$create, 1L)
  expect_equal(mock$calls$bulk, 1L)
  expect_equal(mock$calls$sigma, 1L)
  expect_equal(mock$calls$pixels, 3L)
  expect_identical(as.character(mock$calls$reference_types), unname(x$types))
  expect_equal(fit$weights["b1", "B"], 0.75)
  expect_identical(fit$annotations$dominant_type[1L], "B")
  expect_identical(fit$annotations$status, c("ok", "ambiguous", "not_fit", "not_fit", "filtered"))
  expect_true(all(is.na(fit$annotations$dominant_type[2:5])))
  expect_equal(as.numeric(Matrix::rowSums(fit$weights)), c(1, 1, 0, 0, 0))
  expect_equal(nrow(fit$fit@results$weights), 3L)
  expect_true(fit$fit@internal_vars$cell_types_assigned)
  expect_equal(fit$provenance$bin_size, 32)
  # The same global fit and pixel estimates must be independent of chunk size.
  other <- run_visium_hd_rctd(x$hd, x$reference, x$types, batch_size = NULL)
  expect_equal(other$weights, fit$weights)
  expect_identical(other$annotations, fit$annotations)
})

test_that("RCTD keeps complete-library UMI totals before intersecting genes", {
  x <- .test_hd_rctd_fixture()
  mock <- .test_hd_rctd_backend(drop_bin = character())
  local_mocked_bindings(.hd_rctd_backend = function() mock$backend)
  extra <- Matrix::Matrix(matrix(c(23, 29, 31, 37), 1L,
    dimnames = list("reference_only", colnames(x$reference))), sparse = TRUE)
  reference <- rbind(x$reference, extra)
  extra <- Matrix::Matrix(matrix(c(7, 11, 13, 17, 0), 1L,
    dimnames = list("spatial_only", colnames(x$hd$counts))), sparse = TRUE)
  x$hd$counts <- rbind(x$hd$counts, extra)
  fit <- run_visium_hd_rctd(x$hd, reference, x$types)
  expect_equal(mock$calls$reference_nUMI, Matrix::colSums(reference))
  expect_equal(mock$calls$spatial_nUMI, Matrix::colSums(x$hd$counts)[1:4])
  expect_equal(fit$diagnostics$shared_genes, 12L)
})

test_that("RCTD output types follow the actual fitted profile subset", {
  x <- .test_hd_rctd_fixture()
  x$reference <- cbind(x$reference, r5 = x$reference[, 1L])
  x$types <- c(x$types, r5 = "C")
  mock <- .test_hd_rctd_backend()
  local_mocked_bindings(.hd_rctd_backend = function() mock$backend)
  fit <- run_visium_hd_rctd(
    x$hd, x$reference, x$types,
    create_args = list(cell_type_names = c("A", "B"))
  )
  expect_identical(colnames(fit$weights), c("A", "B"))
  expect_equal(mock$calls$create_args$cell_type_names, c("A", "B"))
})

test_that("Native doublet pairs and convergence survive without retaining the fit", {
  x <- .test_hd_rctd_fixture()
  mock <- .test_hd_rctd_backend()
  local_mocked_bindings(.hd_rctd_backend = function() mock$backend)
  fit <- run_visium_hd_rctd(
    x$hd, x$reference, x$types, doublet_mode = "doublet",
    batch_size = 1L, keep_fit = FALSE,
    run_args = list(CONFIDENCE_THRESHOLD = 9)
  )
  expect_null(fit$fit)
  expect_identical(names(fit$classification), c("results_df", "weights_doublet"))
  expect_identical(rownames(fit$classification$results_df), c("b1", "b2", "b4"))
  expect_identical(rownames(fit$classification$weights_doublet), c("b1", "b2", "b4"))
  expect_equal(as.numeric(fit$classification$weights_doublet["b1", ]), c(.3, .7))
  expect_equal(as.numeric(fit$weights["b1", ]), c(1, 0))
  expect_equal(as.numeric(fit$classification$weights_doublet["b4", ]), c(.35, .65))
  expect_equal(as.numeric(fit$weights["b4", ]), c(0, 0))
  expect_identical(as.character(fit$classification$results_df["b4", "first_type"]), "A")
  expect_true(is.na(fit$annotations["b4", "first_type"]))
  expect_identical(fit$classification$results_df$conv_all, c(TRUE, TRUE, FALSE))
  expect_identical(fit$classification$results_df$conv_doublet, c(TRUE, FALSE, TRUE))
  expect_identical(fit$classification$results_df$min_score, c(1, 2, 4))
  expect_true(is.factor(fit$classification$results_df$first_type))
  expect_identical(fit$parameters$weight_units,
                    "relative RNA/UMI contribution; not cell-number fraction")
  expect_match(fit$parameters$weight_policy, "derived_annotation_weights")
  expect_identical(fit$parameters$upstream_config$RCTDmode, "doublet")
  expect_equal(fit$parameters$upstream_config$UMI_min, 100)
  expect_equal(fit$parameters$upstream_config$CONFIDENCE_THRESHOLD, 9)
  expect_false("weights" %in% names(fit$classification))
  unbatched <- run_visium_hd_rctd(
    x$hd, x$reference, x$types, doublet_mode = "doublet", batch_size = NULL
  )
  expect_equal(fit$classification, unbatched$classification)
})

test_that("Reference diagnostics precede the plotting-reference subsample", {
  x <- .test_hd_rctd_fixture()
  x$reference[, "r4"] <- 1
  mock <- .test_hd_rctd_backend()
  local_mocked_bindings(.hd_rctd_backend = function() mock$backend)
  fit <- run_visium_hd_rctd(x$hd, x$reference, x$types, keep_fit = TRUE)
  reference <- fit$diagnostics$reference
  expect_equal(reference$input_cells, 4L)
  expect_equal(reference$retained_cells, 3L)
  expect_equal(reference$removed_cells, 1L)
  expect_equal(as.integer(reference$input_cells_per_type), c(2L, 2L))
  expect_equal(as.integer(reference$retained_cells_per_type), c(2L, 1L))
  expect_identical(names(reference$retained_cells_per_type), c("A", "B"))
  expect_equal(reference$retained_umi_range, c(240, 240))
  expect_match(reference$stage, "before create.RCTD")
  expect_equal(ncol(fit$fit@reference@counts), 2L)
  expect_null(fit$classification)
  expect_match(fit$parameters$weight_policy, "RNA_UMI_contributions")
})

test_that("Multi-mode native mixture confidence is retained without unrestricted weights", {
  x <- .test_hd_rctd_fixture()
  mock <- .test_hd_rctd_backend()
  local_mocked_bindings(.hd_rctd_backend = function() mock$backend)
  fit <- run_visium_hd_rctd(
    x$hd, x$reference, x$types, doublet_mode = "multi", batch_size = 1L
  )
  native <- fit$classification$pixel_results
  expect_identical(names(native), c("b1", "b2", "b4"))
  expect_identical(native$b2$conf_list, c(TRUE, FALSE))
  expect_false(native$b2$conv_sub)
  expect_equal(native$b2$sub_weights, c(1, 3))
  expect_equal(native$b2$min_score, 2L)
  expect_false("all_weights" %in% names(native$b2))
  expect_equal(as.numeric(fit$weights["b2", ]), c(.25, .75))
  expect_identical(fit$annotations["b2", "status"], "ambiguous")
  expect_null(fit$fit)
})

test_that("Doublet parsing uses pair weights and preserves rejected and uncertain status", {
  ids <- paste0("b", 1:5)
  results <- list(
    results_df = data.frame(
      spot_class = c("singlet", "doublet_certain", "doublet_uncertain", "reject", "doublet_certain"),
      first_type = c("A", "B", "A", "A", "A"),
      second_type = c("B", "A", "B", "B", "B"), row.names = ids),
    weights_doublet = matrix(c(0.3, 0.8, 0.6, NA, 0.5, 0.7, 0.2, 0.4, NA, 0.5),
      5L, dimnames = list(ids, c("first_type", "second_type"))),
    weights = matrix(99, 5L, 2L, dimnames = list(ids, c("A", "B")))
  )
  parsed <- .hd_rctd_parse(results, rev(ids), c("A", "B"), "doublet")
  ann <- parsed$annotations[ids, ]
  expect_equal(as.numeric(parsed$weights["b1", ]), c(1, 0))
  expect_equal(as.numeric(parsed$weights["b2", ]), c(0.2, 0.8))
  expect_equal(as.numeric(parsed$weights["b3", ]), c(0.6, 0.4))
  expect_equal(as.numeric(parsed$weights["b4", ]), c(0, 0))
  expect_identical(ann$status, c("ok", "ok", "ambiguous", "reject", "ambiguous"))
  expect_identical(ann$dominant_type[1:2], c("A", "B"))
  expect_true(all(is.na(ann$dominant_type[3:5])))
  expect_true(is.na(ann$first_type[4L]))
  expect_true(is.na(ann$second_type[1L]))
})

test_that("Multi mode keeps estimated mixtures while auditing uncertain assignments", {
  results <- list(
    list(cell_type_list = c("B", "A"), sub_weights = c(0.7, 0.3), conf_list = c(TRUE, TRUE)),
    list(cell_type_list = c("A", "B"), sub_weights = c(0.8, 0.2), conf_list = c(TRUE, FALSE))
  )
  parsed <- .hd_rctd_parse(results, c("b1", "b2"), c("A", "B"), "multi")
  expect_equal(as.matrix(parsed$weights), matrix(c(0.3, 0.8, 0.7, 0.2), 2L,
    dimnames = list(c("b1", "b2"), c("A", "B"))))
  expect_identical(parsed$annotations$status, c("ok", "ambiguous"))
  expect_identical(parsed$annotations$dominant_type, c("B", NA_character_))
})

test_that("Invalid upstream weights fail rather than being silently relabelled", {
  w <- matrix(c(1, 2, 3, 4), 2L,
    dimnames = list(c("b1", "b2"), c("A", "B")))
  expect_error(.hd_rctd_normalize(w, c("b1", "b3"), c("A", "B")), "mismatched")
  w[1L] <- NaN
  expect_error(.hd_rctd_normalize(w, c("b1", "b2"), c("A", "B")), "invalid weights")
})

test_that("Installed original spacexr can fit synthetic Visium HD bins", {
  skip_if_not_installed("spacexr")
  skip_if_not(all(c("Reference", "create.RCTD", "fitPixels") %in% getNamespaceExports("spacexr")),
               "The original dmcable/spacexr API is required.")
  set.seed(5)
  genes <- paste0("marker", seq_len(60L))
  profiles <- cbind(A = c(rep(25, 30), rep(1, 30)), B = c(rep(1, 30), rep(25, 30)))
  reference <- matrix(stats::rpois(60L * 100L,
    as.vector(profiles[, rep(1:2, each = 50L)])), 60L,
    dimnames = list(genes, paste0("ref", seq_len(100L))))
  spatial <- matrix(stats::rpois(60L * 20L,
    as.vector(profiles[, rep(1:2, each = 10L)])), 60L,
    dimnames = list(genes, paste0("bin", seq_len(20L))))
  hd <- structure(list(counts = Matrix::Matrix(spatial, sparse = TRUE),
    coordinates = matrix(seq_len(40L), 20L, 2L,
      dimnames = list(colnames(spatial), c("x", "y"))), bin_size = 32),
    class = "SpatialESSVisiumHD")
  fit <- run_visium_hd_rctd(hd, reference,
    setNames(rep(c("A", "B"), each = 50L), colnames(reference)),
    batch_size = 7L, create_args = list(UMI_min_sigma = 100), seed = 5L)
  expect_equal(dim(fit$weights), c(20L, 2L))
  expect_equal(as.numeric(Matrix::rowSums(fit$weights)), rep(1, 20L), tolerance = 1e-8)
  expect_equal(fit$diagnostics$batches, 3L)
})
