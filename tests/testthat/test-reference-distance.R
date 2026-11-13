test_that("CSR distances round squares before summation like the reference", {
  coordinates <- matrix(c(
    1727.2685852050781, 373.58755493164062,
    1713.0665710449218, 343.14104156494142
  ), ncol = 2, byrow = TRUE)
  graph <- build_radius_graph_csr(coordinates, radius = 40, store_distance = TRUE)
  edges <- materialize_csr_edges(graph)
  dx <- coordinates[edges$sender, 1L] - coordinates[edges$receiver, 1L]
  dy <- coordinates[edges$sender, 2L] - coordinates[edges$receiver, 2L]
  expect_identical(edges$distance, sqrt(dx * dx + dy * dy))
})

test_that("the primary graph scales distances after strict raw-coordinate filtering", {
  coordinates <- matrix(c(
    1727.2685852050781, 373.58755493164062,
    1713.0665710449218, 343.14104156494142,
    1727.2685852050781, 373.58755493164062,
    1720.41455078125, 350.30490112304688
  ), ncol = 2, byrow = TRUE, dimnames = list(paste0("c", 1:4), c("x", "y")))
  expression <- Matrix::Matrix(rbind(L = c(.2, .6, .8, .4), R = c(.5, 1, .3, .7)),
    sparse = TRUE, dimnames = list(c("L", "R"), rownames(coordinates)))
  lr <- data.frame(ligand = "L", receptor = "R", interaction_name = "L_R",
    annotation = "Secreted Signaling")
  complex <- data.frame(subunit_1 = character(), row.names = character())
  cofactor <- data.frame(cofactor1 = character(), row.names = character())
  run <- function(ratio, interaction.range = 30) spatialess(
    expression, coordinates, factor(LETTERS[1:4]), lr, complex, cofactor,
    ratio = ratio, tol = 0, interaction.range = interaction.range,
    scale.distance = 1, min.percent = 0, min.cells.sr = 1,
    use.AGAN = FALSE, nboot = 0)
  ratio <- .7
  result <- run(ratio)
  edges <- materialize_csr_edges(result$graph)
  dx <- outer(coordinates[, 1L], coordinates[, 1L], `-`)
  dy <- outer(coordinates[, 2L], coordinates[, 2L], `-`)
  raw <- sqrt(dx * dx + dy * dy)
  expected <- raw[cbind(edges$sender, edges$receiver)] * ratio
  expect_identical(edges$distance, as.double(expected))
  expect_false(any(edges$sender == edges$receiver))
  expect_false(any((edges$sender == 1L & edges$receiver == 3L) |
                  (edges$sender == 3L & edges$receiver == 1L)))
  expect_identical(result$graph$parameters$radius, 30)
  expect_identical(result$graph$parameters$scale, 15)
  expect_true(attr(validate_spatial_csr_graph(result$graph), "valid"))
  expect_identical(result$diagnostics$minimum_nonzero_distance,
    min(raw[raw > 0 & raw <= 30 / ratio]) * ratio)

  # The generic graph's numerical candidate tolerance must not add an edge
  # whose raw Euclidean distance exceeds the official threshold.
  threshold <- raw[1L, 2L] - 1e-14
  boundary <- run(1, interaction.range = threshold)
  boundary_edges <- materialize_csr_edges(boundary$graph)
  expect_false(any((boundary_edges$sender == 1L & boundary_edges$receiver == 2L) |
                  (boundary_edges$sender == 2L & boundary_edges$receiver == 1L)))
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(run(bad), "ratio must")
  }
  expect_error(run(.Machine$double.xmin, 100), "converted interaction radius")
})
