test_that("streamed CSR LR aggregation matches materialized edge scoring", {
  coords <- rbind(c1 = c(0, 0), c2 = c(1, 0), c3 = c(2, 0), c4 = c(3, 0))
  csr <- build_radius_graph_csr(coords, radius = 1.01, weight = "gaussian",
                                scale = 0.8, store_distance = TRUE)
  expanded <- build_radius_graph(coords, radius = 1.01, weight = "gaussian",
                                 scale = 0.8)
  expression <- Matrix::Matrix(
    rbind(L = c(2, 0, 4, 0), R = c(0, 3, 0, 5)), sparse = TRUE,
    dimnames = list(c("L", "R"), rownames(coords))
  )
  lr <- data.frame(interaction_name = "L_R", ligand = "L", receptor = "R")
  group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))

  streamed <- aggregate_simple_lr_csr(
    expression, csr, lr, group = group, normalize = FALSE, Kh = 0.5
  )
  materialized <- score_simple_lr_edges(
    expression, expanded, lr, normalize = FALSE, Kh = 0.5
  )$edges
  materialized$sender_group <- as.integer(group[materialized$sender])
  materialized$receiver_group <- as.integer(group[materialized$receiver])
  expected <- aggregate(
    cbind(supported_edges = rep(1, nrow(materialized)),
          sum_hill = materialized$hill,
          sum_spatial_score = materialized$spatial_score),
    by = list(sender_group = materialized$sender_group,
              receiver_group = materialized$receiver_group),
    FUN = sum
  )
  expected$mean_spatial_score <- expected$sum_spatial_score / expected$supported_edges
  observed <- streamed$group_pairs[, c(
    "sender_group", "receiver_group", "supported_edges",
    "sum_hill", "sum_spatial_score", "mean_spatial_score"
  )]
  expected <- expected[order(expected$sender_group, expected$receiver_group), ]
  observed <- observed[order(observed$sender_group, observed$receiver_group), ]
  rownames(expected) <- rownames(observed) <- NULL

  expect_equal(observed, expected, tolerance = 1e-14)
  expect_equal(streamed$diagnostics$supported_edges, nrow(materialized))
  expect_lte(streamed$diagnostics$visited_out_edges, length(csr$neighbors))
})

test_that("streamed CSR aggregation handles zero support without dense output", {
  coords <- rbind(c1 = c(0, 0), c2 = c(1, 0))
  csr <- build_radius_graph_csr(coords, radius = 2)
  expression <- Matrix::Matrix(
    rbind(L = c(1, 0), R = c(0, 0)), sparse = TRUE,
    dimnames = list(c("L", "R"), rownames(coords))
  )
  lr <- data.frame(ligand = "L", receptor = "R")
  result <- aggregate_simple_lr_csr(expression, csr, lr,
                                    group = factor(c("A", "B")),
                                    normalize = FALSE)
  expect_equal(nrow(result$group_pairs), 0)
  expect_equal(result$diagnostics$supported_edges, 0)
  expect_equal(result$diagnostics$sparse_group_pairs, 0)
})
