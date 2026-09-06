test_that("record-sliced permutation exactly matches full matrix reference", {
  fixture <- make_stratified_analytic_fixture()
  full <- permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    nperm = 30L, seed = 719L, finite_correction = TRUE,
    tail = "greater_equal", verbose = FALSE
  )
  selected_rows <- c(2L, 4L)
  selected <- full$group_pairs[
    selected_rows, c("support_pair_index", "lr_index"), drop = FALSE
  ]
  sliced <- permutation_cellchat_records(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks, selected,
    nperm = 30L, seed = 719L, finite_correction = TRUE,
    tail = "greater_equal", verbose = FALSE
  )
  expect_equal(
    sliced$group_pairs$probability,
    full$group_pairs$probability[selected_rows], tolerance = 0
  )
  expect_identical(
    sliced$group_pairs$n_reject,
    full$group_pairs$n_reject[selected_rows]
  )
  expect_equal(
    sliced$group_pairs$pvalue,
    full$group_pairs$pvalue[selected_rows], tolerance = 0
  )
})

test_that("adaptive record slicing never exceeds full matrix work", {
  fixture <- make_stratified_analytic_fixture()
  result <- experimental_adaptive_permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    min_permutations = 20L, batch_size = 20L,
    max_permutations = 80L, seed = 821L,
    confidence_level = 0.95
  )
  expect_lte(
    result$diagnostics$record_permutation_evaluations,
    result$diagnostics$full_matrix_record_permutation_evaluations
  )
  expect_true(result$parameters$multiple_testing != "")
})
