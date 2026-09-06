test_that("two-cell steady-state field matches the analytic solution", {
  coords <- rbind(a = c(0, 0), b = c(1, 0))
  graph <- build_radius_graph(coords, radius = 2, weight = "binary")
  field <- solve_ligand_field(graph, source = c(a = 1, b = 0),
                              diffusion = 1, decay = 1,
                              tolerance = 1e-12)
  expect_true(field$converged)
  expect_equal(unname(field$concentration), c(2 / 3, 1 / 3), tolerance = 1e-11)
  expect_lt(field$relative_residual, 1e-11)
  expect_lt(field$mass_balance_relative_error, 1e-11)
  expect_gte(field$minimum_concentration, -1e-12)
})

test_that("zero diffusion reduces to independent reaction equations", {
  coords <- rbind(a = c(0, 0), b = c(1, 0))
  graph <- build_radius_graph(coords, radius = 2)
  field <- solve_ligand_field(
    graph,
    source = c(a = 3, b = 4),
    receptor = c(a = 1, b = 2),
    diffusion = 0, decay = 1, uptake = 2, production_rate = 2,
    tolerance = 1e-13
  )
  expect_equal(unname(field$concentration), c(2, 1.6), tolerance = 1e-12)
  expect_true(field$converged)
})

test_that("compartment barriers block ligand propagation", {
  coords <- rbind(a = c(0, 0), b = c(1, 0), c = c(1.5, 0))
  graph <- build_radius_graph(coords, radius = 2,
                              compartment = c("left", "left", "right"),
                              same_compartment = TRUE)
  field <- solve_ligand_field(graph, source = c(a = 1, b = 0, c = 0),
                              diffusion = 1, decay = 1)
  expect_gt(field$concentration[["b"]], 0)
  expect_equal(field$concentration[["c"]], 0, tolerance = 1e-14)
})

test_that("named field inputs are aligned to graph cells", {
  coords <- rbind(a = c(0, 0), b = c(1, 0), c = c(2, 0))
  graph <- build_radius_graph(coords, radius = 1.1)
  ordered <- solve_ligand_field(graph, source = c(a = 1, b = 0, c = 0))
  shuffled <- solve_ligand_field(graph, source = c(c = 0, a = 1, b = 0))
  expect_equal(ordered$concentration, shuffled$concentration)
})

test_that("zero source returns an exact zero field", {
  coords <- rbind(a = c(0, 0), b = c(1, 0))
  graph <- build_radius_graph(coords, radius = 2)
  field <- solve_ligand_field(graph, source = c(0, 0))
  expect_true(field$converged)
  expect_equal(field$iterations, 0)
  expect_equal(field$concentration, c(a = 0, b = 0))
})
