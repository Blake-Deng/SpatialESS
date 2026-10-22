make_tied_permutation_fixture <- function() {
  expression <- Matrix::Matrix(
    matrix(
      1,
      nrow = 2L, ncol = 8L,
      dimnames = list(c("L", "R"), paste0("cell", 1:8))
    ),
    sparse = TRUE
  )
  group <- factor(rep(c("A", "B"), each = 4L), levels = c("A", "B"))
  coords <- cbind(x = seq_len(8), y = 0)
  rownames(coords) <- colnames(expression)
  graph <- build_radius_graph_csr(coords, radius = 20, weight = "binary")
  support <- build_group_support_csr(graph, group)
  blocks <- build_spatial_blocks(coords, block_size = 100, origin = c(0, 0))
  complex <- data.frame(subunit_1 = character(), row.names = character())
  cofactor <- data.frame(cofactor1 = character(), row.names = character())
  lr <- data.frame(interaction_name = "L_R", ligand = "L", receptor = "R")
  components <- prepare_cellchat_lr_components(
    lr, rownames(expression), complex, cofactor
  )
  prepared <- prepare_sparse_trimean(
    expression, genes = components$genes, normalize = FALSE
  )
  list(group = group, support = support, blocks = blocks,
       components = components, prepared = prepared)
}

test_that("calibrated upper tail counts ties while CellChat mode is explicit", {
  fixture <- make_tied_permutation_fixture()
  calibrated <- permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    nperm = 20, seed = 41, finite_correction = TRUE,
    tail = "greater_equal", verbose = FALSE
  )
  compatibility <- permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    nperm = 20, seed = 41, finite_correction = TRUE,
    tail = "strict_greater", verbose = FALSE
  )

  expect_true(all(calibrated$group_pairs$pvalue == 1))
  expect_true(all(compatibility$group_pairs$pvalue == 1 / 21))
  expect_identical(calibrated$parameters$tail, "greater_equal")
  expect_identical(compatibility$parameters$tail, "strict_greater")
})

