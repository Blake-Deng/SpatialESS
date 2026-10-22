test_that("adaptive permutation records per-record stopping and intervals", {
  fixture <- make_stratified_analytic_fixture()
  result <- experimental_adaptive_permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    min_permutations = 20L, batch_size = 20L,
    max_permutations = 60L, seed = 313L,
    confidence_level = 0.95, tail = "greater_equal"
  )
  pairs <- result$group_pairs
  expect_true(all(pairs$permutations_used >= 20L &
                  pairs$permutations_used <= 60L))
  expect_true(all(is.finite(pairs$pvalue) & pairs$pvalue >= 0 &
                  pairs$pvalue <= 1))
  expect_true(all(is.finite(pairs$pvalue_lower) &
                  pairs$pvalue_lower <= pairs$pvalue_upper))
  expect_true(all(pairs$decision %in% c("significant", "not_significant")))
  expect_true(result$diagnostics$batches >= 1L)
})

test_that("adaptive permutation is reproducible and respects max bound", {
  fixture <- make_stratified_analytic_fixture()
  first <- experimental_adaptive_permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    min_permutations = 10L, batch_size = 10L,
    max_permutations = 30L, seed = 417L
  )
  second <- experimental_adaptive_permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    min_permutations = 10L, batch_size = 10L,
    max_permutations = 30L, seed = 417L
  )
  expect_identical(first$group_pairs$pvalue, second$group_pairs$pvalue)
  expect_identical(first$group_pairs$permutations_used,
                   second$group_pairs$permutations_used)
  expect_true(all(first$group_pairs$permutations_used <= 30L))
})
