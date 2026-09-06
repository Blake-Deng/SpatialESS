test_that("CSR field solver matches the analytic two-cell solution", {
  coords <- rbind(a = c(0, 0), b = c(1, 0))
  graph <- build_radius_graph_csr(coords, radius = 2, weight = "binary")
  field <- solve_ligand_field_csr(
    graph, source = c(a = 1, b = 0),
    diffusion = 1, decay = 1, tolerance = 1e-12
  )
  expect_true(field$converged)
  expect_equal(unname(field$concentration), c(2 / 3, 1 / 3), tolerance = 1e-11)
  expect_lt(field$relative_residual, 1e-11)
  expect_lt(field$mass_balance_relative_error, 1e-11)
})

test_that("CSR and materialized graph field solvers agree", {
  set.seed(151)
  coords <- matrix(runif(160), ncol = 2,
                   dimnames = list(sprintf("c%02d", 1:80), c("x", "y")))
  csr <- build_radius_graph_csr(coords, radius = 0.2,
                                weight = "gaussian", scale = 0.1)
  expanded <- build_radius_graph(coords, radius = 0.2,
                                 weight = "gaussian", scale = 0.1)
  source <- setNames(runif(80), rownames(coords))
  receptor <- setNames(runif(80), rownames(coords))
  compact_field <- solve_ligand_field_csr(
    csr, source, receptor, diffusion = 0.7, decay = 0.9,
    uptake = 0.2, tolerance = 1e-11
  )
  expanded_field <- solve_ligand_field(
    expanded, source, receptor, diffusion = 0.7, decay = 0.9,
    uptake = 0.2, tolerance = 1e-11
  )
  expect_true(compact_field$converged)
  expect_equal(compact_field$concentration, expanded_field$concentration,
               tolerance = 1e-11)
  expect_equal(compact_field$mass_balance_relative_error,
               expanded_field$mass_balance_relative_error,
               tolerance = 1e-11)
})

test_that("CSR field solver handles zero source exactly", {
  coords <- rbind(a = c(0, 0), b = c(1, 0))
  graph <- build_radius_graph_csr(coords, radius = 2)
  field <- solve_ligand_field_csr(graph, source = c(0, 0))
  expect_true(field$converged)
  expect_equal(field$iterations, 0)
  expect_equal(field$concentration, c(a = 0, b = 0))
})
