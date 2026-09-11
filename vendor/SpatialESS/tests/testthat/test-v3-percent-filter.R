test_that("V3 formatting preserves the 8/82 boundary and column context", {
  sizes <- c(90, 82, 63)
  positive <- cbind(c(11, 8, 3), c(11, 8, 3))
  mask <- SpatialESS:::.spatialess_v3_percent_mask(positive, sizes, .1)
  expect_true(mask[2, 1])
  expect_false(positive[2, 1] / sizes[2] >= .1)
  expected <- format(data.frame(a = positive[, 1]/sizes,
                                b = positive[, 2]/sizes), digits = 1) >= .1
  expect_identical(mask, unname(expected))
})

test_that("fraction masks match the official aggregate and format expression", {
  old <- options(scipen = 0, OutDec = ".")
  on.exit(options(old), add = TRUE)
  set.seed(217)
  for (trial in seq_len(80)) {
    sizes <- sample(1:120, 7, replace = TRUE)
    positive <- cbind(vapply(sizes, function(n) sample.int(n + 1, 1) - 1L, integer(1)),
                      vapply(sizes, function(n) sample.int(n + 1, 1) - 1L, integer(1)))
    group <- factor(rep(seq_along(sizes), sizes))
    dataLR <- do.call(rbind, lapply(seq_along(sizes), function(i) {
      cbind(seq_len(sizes[i]) <= positive[i, 1], seq_len(sizes[i]) <= positive[i, 2])
    }))
    aggregated <- aggregate(1 * (dataLR > 0), list(group), mean)
    for (threshold in c(0, .01, .1, .5, 1)) {
      expected <- format(aggregated[, -1], digits = 1) >= threshold
      actual <- SpatialESS:::.spatialess_v3_percent_mask(positive, sizes, threshold)
      expect_identical(actual, unname(expected))
    }
  }
})

test_that("format options and very small fractions are not hand rounded", {
  old <- options()
  on.exit(options(old), add = TRUE)
  sizes <- c(1000000, 82, 1, 10000)
  positive <- cbind(c(1, 8, 1, 0), c(0, 9, 1, 1))
  for (decimal in c(".", ",")) for (scipen in c(-9, 0, 999)) {
    options(OutDec = decimal, scipen = scipen)
    fractions <- data.frame(a = positive[, 1]/sizes, b = positive[, 2]/sizes)
    for (threshold in c(0, .00001, .1, 1)) {
      expect_identical(SpatialESS:::.spatialess_v3_percent_mask(positive, sizes, threshold),
                       unname(format(fractions, digits = 1) >= threshold))
    }
    expect_true(all(SpatialESS:::.spatialess_v3_percent_mask(positive, sizes, 0)))
  }
})
