test_that("spatial blocks preserve sample and compartment boundaries", {
  coords <- cbind(x = c(0, 1, 2, 3, 0, 1, 2, 3), y = 0)
  rownames(coords) <- paste0("c", 1:8)
  sample_id <- rep(c("s1", "s2"), each = 4)
  compartment <- rep(c("x", "x", "y", "y"), 2)
  blocks <- build_spatial_blocks(
    coords, block_size = 100, sample_id = sample_id,
    compartment = compartment, origin = c(0, 0)
  )
  permutation <- permute_within_spatial_blocks(blocks, seed = 19)

  expect_setequal(unname(permutation), seq_len(nrow(coords)))
  expect_equal(blocks$block_code[permutation], blocks$block_code)
  expect_equal(sample_id[permutation], sample_id)
  expect_equal(compartment[permutation], compartment)
  expect_identical(
    permutation, permute_within_spatial_blocks(blocks, seed = 19)
  )
})

test_that("singleton spatial blocks remain fixed", {
  coords <- cbind(x = c(0, 100, 200), y = 0)
  rownames(coords) <- c("a", "b", "c")
  blocks <- build_spatial_blocks(coords, block_size = 10, origin = c(0, 0))
  expect_identical(
    unname(permute_within_spatial_blocks(blocks, seed = 1)), 1:3
  )
})

test_that("permuted sparse triMean equals explicit profile reassignment", {
  expression <- Matrix::Matrix(
    rbind(L = c(0, 1, 2, 4, 8, 0, 3, 6),
          R = c(5, 4, 3, 2, 1, 0, 7, 8)),
    sparse = TRUE,
    dimnames = list(c("L", "R"), paste0("c", 1:8))
  )
  group <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
  coords <- cbind(x = seq_len(8), y = 0)
  rownames(coords) <- colnames(expression)
  blocks <- build_spatial_blocks(coords, block_size = 100, origin = c(0, 0))
  permutation <- permute_within_spatial_blocks(blocks, seed = 31)
  prepared <- prepare_sparse_trimean(expression, normalize = FALSE)
  observed <- summarize_prepared_trimean(
    prepared, group, source_for_target = permutation
  )$average
  expected <- summarize_cellchat_trimean(
    expression[, permutation, drop = FALSE], group, normalize = FALSE
  )$average
  expect_equal(observed, expected, tolerance = 1e-15)
})

make_permutation_fixture <- function() {
  expression <- Matrix::Matrix(
    rbind(L = c(4, 5, 6, 7, 1, 2, 3, 4),
          R = c(1, 2, 3, 4, 5, 6, 7, 8)),
    sparse = TRUE,
    dimnames = list(c("L", "R"), paste0("c", 1:8))
  )
  group <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
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
  list(expression = expression, group = group, coords = coords,
       support = support, blocks = blocks, components = components,
       prepared = prepared)
}

test_that("spatial permutation inference is reproducible and score-compatible", {
  fixture <- make_permutation_fixture()
  first <- permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    nperm = 20, seed = 101, verbose = FALSE
  )
  second <- permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    nperm = 20, seed = 101, verbose = FALSE
  )
  average <- summarize_prepared_trimean(fixture$prepared, fixture$group)
  direct <- score_cellchat_group_support(
    average, fixture$support, fixture$components
  )$group_pairs

  expect_equal(first$group_pairs$pvalue, second$group_pairs$pvalue)
  expect_equal(first$group_pairs$n_reject, second$group_pairs$n_reject)
  expect_true(all(first$group_pairs$pvalue >= 0 &
                  first$group_pairs$pvalue <= 1))
  expect_gt(first$diagnostics$exchangeable_cell_fraction, 0)
  expect_gt(mean(first$diagnostics$moved_group_fraction), 0)
  key_first <- paste(first$group_pairs$lr_index,
                     first$group_pairs$sender_group,
                     first$group_pairs$receiver_group)
  key_direct <- paste(direct$lr_index, direct$sender_group,
                      direct$receiver_group)
  expect_equal(
    first$group_pairs$probability,
    direct$molecular_probability[match(key_first, key_direct)],
    tolerance = 1e-15
  )
})

test_that("permutation rejects blocks nested entirely inside target groups", {
  fixture <- make_permutation_fixture()
  fixture$blocks <- build_spatial_blocks(
    fixture$coords, block_size = 4, origin = c(1, 0)
  )
  expect_error(
    permutation_cellchat_group_support(
      fixture$prepared, fixture$group, fixture$support,
      fixture$components, fixture$blocks, nperm = 2, verbose = FALSE
    ),
    "No spatial block spans multiple target groups"
  )
})
