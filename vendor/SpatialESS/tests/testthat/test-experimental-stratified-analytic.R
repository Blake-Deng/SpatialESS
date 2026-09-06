test_that("stratified analytic evaluates spatially shared paracrine LR", {
  fixture <- make_stratified_analytic_fixture()
  result <- experimental_stratified_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks
  )
  pairs <- result$group_pairs
  paracrine <- pairs$sender_group != pairs$receiver_group
  expect_true(all(pairs$inference_mode[paracrine] ==
                  "experimental_stratified_analytic"))
  expect_true(all(is.finite(pairs$pvalue[paracrine])))
  expect_true(all(pairs$shared_blocks[paracrine] == 8L))
  expect_true(all(pairs$fallback_reason[!paracrine] == "autocrine_pair"))
})

test_that("stratified analytic rejects blocks nested inside groups", {
  fixture <- make_stratified_analytic_fixture(nested = TRUE)
  result <- experimental_stratified_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks
  )
  paracrine <- result$group_pairs$sender_group !=
    result$group_pairs$receiver_group
  expect_true(all(result$group_pairs$fallback_reason[paracrine] ==
                  "insufficient_shared_blocks"))
  expect_equal(result$diagnostics$analytic_records, 0)
})

test_that("stratified analytic guards zero-inflated quartiles", {
  fixture <- make_stratified_analytic_fixture(tied = TRUE)
  result <- experimental_stratified_analytic_cellchat_group_support(
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


test_that("stratified analytic guards moderate zero inflation", {
  fixture <- make_stratified_analytic_fixture(mild_zero = TRUE)
  result <- experimental_stratified_analytic_cellchat_group_support(
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

test_that("stratified analytic is reproducible", {
  fixture <- make_stratified_analytic_fixture()
  first <- experimental_stratified_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks
  )
  second <- experimental_stratified_analytic_cellchat_group_support(
    fixture$prepared, fixture$group, fixture$support,
    fixture$components, fixture$blocks
  )
  expect_identical(first$group_pairs$pvalue, second$group_pairs$pvalue)
  expect_identical(first$group_pairs$fallback_reason,
                   second$group_pairs$fallback_reason)
})
