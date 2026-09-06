test_that("simple LR edge score follows the declared Hill formula", {
  coords <- rbind(c1 = c(0, 0), c2 = c(1, 0), c3 = c(2, 0), c4 = c(3, 0))
  graph <- build_radius_graph(coords, radius = 1.01, weight = "binary")
  expression <- Matrix::Matrix(
    rbind(L = c(2, 0, 4, 0), R = c(0, 3, 0, 5)), sparse = TRUE,
    dimnames = list(c("L", "R"), rownames(coords))
  )
  lr <- data.frame(interaction_name = "L_R", ligand = "L", receptor = "R")
  result <- score_simple_lr_edges(expression, graph, lr, Kh = 0.5,
                                  n = 1, normalize = FALSE)
  expect_identical(paste(result$edges$sender, result$edges$receiver),
                   c("1 2", "3 2", "3 4"))
  expected_product <- c(6, 12, 20)
  expect_equal(result$edges$molecular_product, expected_product)
  expect_equal(result$edges$hill, expected_product / (0.5 + expected_product))
  expect_equal(result$edges$spatial_score, result$edges$hill)
  expect_equal(result$diagnostics$supported_edges, 3)
})

test_that("expression is aligned to graph cell identifiers", {
  coords <- rbind(c1 = c(0, 0), c2 = c(1, 0), c3 = c(2, 0))
  graph <- build_radius_graph(coords, radius = 1.1)
  expression <- Matrix::Matrix(
    rbind(L = c(1, 0, 1), R = c(0, 1, 0)), sparse = TRUE,
    dimnames = list(c("L", "R"), rownames(coords))
  )
  lr <- data.frame(ligand = "L", receptor = "R")
  reference <- score_simple_lr_edges(expression, graph, lr, normalize = FALSE)
  reordered <- score_simple_lr_edges(expression[, c(3, 1, 2)], graph, lr,
                                     normalize = FALSE)
  expect_equal(reference$edges, reordered$edges)
})

test_that("unsupported and non-simple pairs are excluded", {
  coords <- rbind(c1 = c(0, 0), c2 = c(1, 0))
  graph <- build_radius_graph(coords, radius = 2)
  expression <- Matrix::Matrix(
    rbind(L = c(0, 1), R = c(0, 0)), sparse = TRUE,
    dimnames = list(c("L", "R"), rownames(coords))
  )
  lr <- data.frame(interaction_name = c("L_R", "complex_R"),
                   ligand = c("L", "L_complex"), receptor = c("R", "R"))
  result <- score_simple_lr_edges(expression, graph, lr, normalize = FALSE)
  expect_equal(nrow(result$lr), 1)
  expect_equal(nrow(result$edges), 0)
  expect_equal(result$diagnostics$supported_edges, 0)
})
