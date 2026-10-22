test_that("hybrid uses analytic inference only to screen non-significance", {
  fixture <- make_stratified_analytic_fixture()
  fixture$support$group_pairs <- subset(
    fixture$support$group_pairs,
    sender_group != receiver_group
  )
  result <- experimental_hybrid_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    alpha = 1e-12, screening_margin = 0, nperm = 20, seed = 91
  )
  expect_false(result$parameters$permutation_called)
  expect_true(all(result$group_pairs$inference_mode ==
                  "experimental_hybrid_screened_nonsignificant"))
  expect_true(all(result$group_pairs$pvalue == 1))
  expect_true(all(is.finite(result$group_pairs$analytic_pvalue)))
  expect_false(any(result$group_pairs$permutation_confirmed))
})

test_that("hybrid replaces unsupported records with permutation reference", {
  fixture <- make_stratified_analytic_fixture()
  result <- experimental_hybrid_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    alpha = 1e-12, screening_margin = 0, nperm = 20, seed = 91
  )
  direct <- permutation_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    nperm = 20, seed = 91, finite_correction = TRUE,
    tail = "greater_equal", verbose = FALSE
  )
  key_result <- paste(result$group_pairs$lr_index,
                      result$group_pairs$support_pair_index)
  key_direct <- paste(direct$group_pairs$lr_index,
                      direct$group_pairs$support_pair_index)
  match_index <- match(key_result, key_direct)
  fallback <- result$group_pairs$inference_mode !=
    "experimental_hybrid_screened_nonsignificant"
  expect_true(any(fallback))
  expect_equal(
    result$group_pairs$pvalue[fallback],
    direct$group_pairs$pvalue[match_index[fallback]]
  )
  expect_true(all(result$group_pairs$n_reject[fallback] ==
                  direct$group_pairs$n_reject[match_index[fallback]]))
  expect_true(all(result$group_pairs$permutation_confirmed[fallback]))
})

test_that("every potentially significant hybrid result is permutation-confirmed", {
  fixture <- make_stratified_analytic_fixture()
  result <- experimental_hybrid_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    alpha = 0.05, screening_margin = 0.01, nperm = 99, seed = 91
  )
  significant <- result$group_pairs$pvalue <= result$parameters$alpha
  expect_true(all(result$group_pairs$permutation_confirmed[significant]))
  expect_false(any(
    result$group_pairs$inference_mode[significant] ==
      "experimental_hybrid_screened_nonsignificant"
  ))
})


test_that("analytic candidates always enter fixed permutation", {
  fixture <- make_stratified_analytic_fixture()
  result <- experimental_hybrid_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks,
    alpha = 0.2, screening_margin = 0, nperm = 99, seed = 91
  )
  candidate <- result$group_pairs$inference_mode ==
    "experimental_hybrid_permutation_candidate"
  expect_true(any(candidate))
  expect_true(all(result$group_pairs$permutation_confirmed[candidate]))
  expect_true(all(is.finite(result$group_pairs$n_reject[candidate])))
})

test_that("hybrid validates the one-way screening boundary", {
  fixture <- make_stratified_analytic_fixture()
  expect_error(
    experimental_hybrid_cellchat_group_support(
      fixture$prepared, fixture$group, fixture$support,
      fixture$components, fixture$blocks,
      alpha = 0.99, screening_margin = 0.01
    ),
    "alpha \\+ screening_margin"
  )
})
