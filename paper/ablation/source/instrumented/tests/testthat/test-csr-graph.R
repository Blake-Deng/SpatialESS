test_that("CSR radius graph matches the existing graph and brute-force oracle", {
  set.seed(101)
  coords <- matrix(runif(240, -4, 4), ncol = 2,
                   dimnames = list(sprintf("c%03d", 1:120), c("x", "y")))
  sample_id <- rep(c("s1", "s2"), each = 60)
  compartment <- rep(rep(c("a", "b"), each = 30), 2)

  csr <- build_radius_graph_csr(
    coords, radius = 1.25, sample_id = sample_id,
    compartment = compartment, same_compartment = TRUE,
    weight = "gaussian", scale = 0.7, store_distance = TRUE
  )
  reference <- build_radius_graph_bruteforce(
    coords, radius = 1.25, sample_id = sample_id,
    compartment = compartment, same_compartment = TRUE,
    weight = "gaussian", scale = 0.7
  )
  compact_edges <- materialize_csr_edges(csr)
  compact_edges <- compact_edges[order(compact_edges$sender, compact_edges$receiver), ]
  reference_edges <- reference$edges[order(reference$edges$sender,
                                           reference$edges$receiver), ]

  expect_identical(compact_edges$sender, reference_edges$sender)
  expect_identical(compact_edges$receiver, reference_edges$receiver)
  expect_equal(compact_edges$distance, reference_edges$distance, tolerance = 1e-14)
  expect_equal(compact_edges$weight, reference_edges$weight, tolerance = 1e-14)
  expect_true(isTRUE(attr(validate_spatial_csr_graph(csr), "valid")))
})

test_that("CSR graph is symmetric and has matching degree offsets", {
  coords <- rbind(a = c(0, 0), b = c(1, 0), c = c(2, 0), d = c(5, 0))
  graph <- build_radius_graph_csr(coords, radius = 1.1, store_distance = FALSE)
  edges <- materialize_csr_edges(graph)
  keys <- paste(edges$sender, edges$receiver, sep = ":")
  reverse_keys <- paste(edges$receiver, edges$sender, sep = ":")
  expect_setequal(keys, reverse_keys)
  expect_equal(diff(graph$offsets), graph$degree)
  expect_length(graph$distances, 0)
})

test_that("CSR graph respects samples, barriers and allocation guard", {
  coords <- rbind(a = c(0, 0), b = c(0.1, 0), c = c(0.2, 0), d = c(0.3, 0))
  graph <- build_radius_graph_csr(
    coords, radius = 1,
    sample_id = c("s1", "s2", "s1", "s1"),
    compartment = c("x", "x", "y", "x"),
    same_compartment = TRUE
  )
  checks <- validate_spatial_csr_graph(graph)
  expect_true(isTRUE(attr(checks, "valid")))
  expect_error(
    build_radius_graph_csr(coords, radius = 1, max_edges = 1),
    "exceeded max_edges"
  )
})

test_that("CSR representation is smaller than materialized edge graph", {
  set.seed(102)
  coords <- matrix(runif(1000), ncol = 2,
                   dimnames = list(sprintf("c%03d", 1:500), c("x", "y")))
  csr <- build_radius_graph_csr(coords, radius = 0.1)
  expanded <- build_radius_graph(coords, radius = 0.1)
  expect_lt(as.numeric(object.size(csr)), as.numeric(object.size(expanded)))
})
