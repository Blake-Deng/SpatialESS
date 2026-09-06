make_analytic_fixture <- function(seed = 11L, tied = FALSE,
                                  correlated = FALSE, mild_zero = FALSE) {
  set.seed(seed)
  cell_count <- 240L
  cell_id <- paste0("cell", seq_len(cell_count))
  group <- factor(rep(c("A", "B"), each = cell_count / 2L),
                  levels = c("A", "B"))
  ligand <- if (tied) {
    rep(c(0, 0, 0, 1), length.out = cell_count)
  } else {
    stats::rgamma(cell_count, shape = 5, rate = 2)
  }
  receptor <- if (correlated) {
    ligand * exp(stats::rnorm(cell_count, sd = 0.01))
  } else if (tied) {
    rep(c(0, 0, 1, 1), length.out = cell_count)
  } else {
    stats::rgamma(cell_count, shape = 6, rate = 2)
  }
  if (mild_zero) {
    zero_index <- seq_len(cell_count) %% 5L == 0L
    ligand[zero_index] <- 0
    receptor[zero_index] <- 0
  }
  expression <- Matrix::Matrix(
    rbind(L = ligand, R = receptor), sparse = TRUE,
    dimnames = list(c("L", "R"), cell_id)
  )
  coords <- cbind(x = seq_len(cell_count), y = 0)
  rownames(coords) <- cell_id
  graph <- build_radius_graph_csr(coords, radius = cell_count,
                                  weight = "binary")
  support <- build_group_support_csr(graph, group)
  blocks <- build_spatial_blocks(
    coords, block_size = cell_count + 1, origin = c(0, 0)
  )
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
       components = components, prepared = prepared, coords = coords)
}

test_that("analytic prototype evaluates only supported paracrine simple LR", {
  fixture <- make_analytic_fixture()
  result <- experimental_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks
  )
  pairs <- result$group_pairs
  paracrine <- pairs$sender_group != pairs$receiver_group
  autocrine <- !paracrine

  expect_true(all(pairs$inference_mode[paracrine] ==
                  "experimental_analytic"))
  expect_true(all(is.finite(pairs$pvalue[paracrine])))
  expect_true(all(pairs$pvalue[paracrine] >= 0 &
                  pairs$pvalue[paracrine] <= 1))
  expect_true(all(pairs$fallback_reason[autocrine] == "autocrine_pair"))
  expect_true(all(is.na(pairs$pvalue[autocrine])))
  expect_identical(
    result$parameters$calibration_status, "experimental_unvalidated"
  )
})

test_that("analytic prototype guards tied zero quartiles", {
  fixture <- make_analytic_fixture(tied = TRUE)
  result <- experimental_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks
  )
  paracrine <- result$group_pairs$sender_group !=
    result$group_pairs$receiver_group
  expect_true(all(result$group_pairs$inference_mode[paracrine] ==
                  "permutation_fallback_required"))
  expect_true(all(grepl(
    "zero_fraction|quartile",
    result$group_pairs$fallback_reason[paracrine]
  )))
})


test_that("analytic prototype guards moderate zero inflation", {
  fixture <- make_analytic_fixture(mild_zero = TRUE)
  result <- experimental_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks
  )
  paracrine <- result$group_pairs$sender_group !=
    result$group_pairs$receiver_group
  expect_true(all(
    result$group_pairs$fallback_reason[paracrine] %in%
      c("ligand_excessive_zero_fraction",
        "receptor_excessive_zero_fraction")
  ))
  expect_equal(result$diagnostics$analytic_records, 0)
})

test_that("analytic prototype guards strong ligand-receptor dependence", {
  fixture <- make_analytic_fixture(correlated = TRUE)
  result <- experimental_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    max_abs_influence_correlation = 0.5
  )
  paracrine <- result$group_pairs$sender_group !=
    result$group_pairs$receiver_group
  expect_true(all(result$group_pairs$fallback_reason[paracrine] ==
                  "strong_ligand_receptor_dependence"))
})

test_that("analytic prototype rejects spatially stratified nulls", {
  fixture <- make_analytic_fixture()
  fixture$blocks <- build_spatial_blocks(
    fixture$coords, block_size = 40, origin = c(0, 0)
  )
  result <- experimental_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks
  )
  expect_true(all(result$group_pairs$fallback_reason ==
                  "spatially_stratified_null"))
  expect_equal(result$diagnostics$analytic_records, 0)
})
