test_that("equal-weight spatial null quantiles exactly reproduce R type 7", {
  set.seed(71)
  probabilities <- seq(0, 1, by = 0.05)
  for (values in list(
      stats::runif(101),
      sample(c(0, 1, 2, 3, 5), 101, replace = TRUE),
      rep(4, 101))) {
    observed <- SpatialESS:::.spatialess_weighted_quantile(
      values, rep(1, length(values)), probabilities
    )
    expected <- as.numeric(stats::quantile(
      values, probabilities, names = FALSE, type = 7
    ))
    expect_equal(observed, expected, tolerance = 0)
  }
})

test_that("weighted spatial null quantiles are finite and monotone", {
  values <- c(0, 1, 2, 3, 5, 8)
  weights <- c(1, 1, 2, 3, 5, 8)
  observed <- SpatialESS:::.spatialess_weighted_quantile(
    values, weights, seq(0, 1, by = 0.01)
  )
  expect_true(all(is.finite(observed)))
  expect_true(all(diff(observed) >= 0))
  expect_equal(observed[c(1L, length(observed))], range(values))
})
