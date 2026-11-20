test_that("explicit V3 filtering matches frozen official observed and permutation results", {
  skip_if_not_installed("future.apply")
  f <- readRDS(test_path("fixtures", "v3_percent_boundary.rds"))
  for (mode in c("default", "relaxed")) {
    filters <- if (mode == "default") list(min.percent = .1, min.cells.sr = 5L) else list(min.percent = 0, min.cells.sr = 1L)
    result <- do.call(spatialess, c(list(expression = f$expression,
      coordinates = f$coordinates, group = f$group, lr = f$lr,
      complex = f$db$complex, cofactor = f$db$cofactor,
      percent_filter = "spatialcellchat_v3"), f$parameters, filters))
    expected <- f[[mode]]
    key <- function(x) paste(x$interaction_name, x$sender_group_name, x$receiver_group_name, sep = "|")
    expect_setequal(key(result$records), key(expected))
    observed <- result$records[match(key(expected), key(result$records)), , drop = FALSE]
    expect_equal(nrow(result$records), nrow(expected))
    expect_lte(max(abs(observed$probability - expected$probability)), 1e-15)
    expect_identical(observed$pvalue, expected$pvalue)
    expect_identical(observed$pvalue < .05, expected$pvalue < .05)
    expect_identical(result$parameters$percentage_filter, "SpatialCellChat_V3_format_digits1")
    expect_identical(result$parameters$percent_filter, "spatialcellchat_v3")
  }
})
