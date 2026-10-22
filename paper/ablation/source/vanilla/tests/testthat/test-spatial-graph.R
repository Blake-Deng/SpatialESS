edge_key <- function(graph) {
  e <- graph$edges
  paste(e$sender_id, e$receiver_id, sep = "->")
}

test_that("spatial hash radius graph matches brute-force oracle", {
  set.seed(11)
  coords <- matrix(runif(160, -5, 5), ncol = 2,
                   dimnames = list(sprintf("c%03d", 1:80), c("x", "y")))
  sample_id <- rep(c("s1", "s2"), each = 40)
  compartment <- rep(rep(c("a", "b"), each = 20), 2)

  for (barrier in c(FALSE, TRUE)) {
    fast <- build_radius_graph(coords, radius = 1.7, sample_id = sample_id,
                               compartment = compartment,
                               same_compartment = barrier,
                               weight = "gaussian", scale = 0.8)
    oracle <- build_radius_graph_bruteforce(
      coords, radius = 1.7, sample_id = sample_id,
      compartment = compartment, same_compartment = barrier,
      weight = "gaussian", scale = 0.8
    )
    expect_identical(edge_key(fast), edge_key(oracle))
    expect_equal(fast$edges$distance, oracle$edges$distance, tolerance = 1e-14)
    expect_equal(fast$edges$weight, oracle$edges$weight, tolerance = 1e-14)
    expect_true(isTRUE(attr(validate_spatial_graph(fast), "valid")))
  }
})

test_that("radius graph respects samples, barriers and self-edge policy", {
  coords <- rbind(a = c(0, 0), b = c(0.1, 0), c = c(0.2, 0), d = c(0.3, 0))
  sample_id <- c("s1", "s2", "s1", "s1")
  compartment <- c("x", "x", "y", "x")
  graph <- build_radius_graph(coords, 1, sample_id, compartment,
                              same_compartment = TRUE)
  expect_false(any(graph$edges$sender == graph$edges$receiver))
  expect_true(all(graph$cells$sample_id[graph$edges$sender] ==
                  graph$cells$sample_id[graph$edges$receiver]))
  expect_true(all(graph$cells$compartment[graph$edges$sender] ==
                  graph$cells$compartment[graph$edges$receiver]))
})

test_that("cell reordering preserves the graph after mapping identifiers", {
  set.seed(29)
  coords <- matrix(runif(100), ncol = 2,
                   dimnames = list(sprintf("cell%02d", 1:50), c("x", "y")))
  g1 <- build_radius_graph(coords, 0.2)
  permutation <- sample.int(nrow(coords))
  g2 <- build_radius_graph(coords[permutation, , drop = FALSE], 0.2)
  expect_setequal(edge_key(g1), edge_key(g2))
})

test_that("kNN graph does not create cross-sample edges", {
  skip_if_not_installed("RANN")
  set.seed(31)
  coords <- matrix(runif(120), ncol = 2,
                   dimnames = list(sprintf("cell%02d", 1:60), c("x", "y")))
  sample_id <- rep(c("left", "right"), each = 30)
  graph <- build_knn_graph(coords, k = 5, sample_id = sample_id)
  checks <- validate_spatial_graph(graph)
  expect_true(unname(checks[["no_cross_sample_edges"]]))
  expect_true(isTRUE(attr(checks, "valid")))
})

