test_that("sparse group triMeans match the R type-7 arithmetic exactly", {
  reference <- function(x) mean(stats::quantile(
    x, c(.25, .5, .5, .75), type = 7, names = FALSE
  ))
  set.seed(1706)
  for (size in c(1L, 2L, 3L, 4L, 5L, 6L, 9L, 21L, 100L)) {
    for (iteration in seq_len(20L)) {
      values <- c(0, 10^runif(size, -15, 2))
      expression <- Matrix::Matrix(matrix(values, 1L), sparse = TRUE,
        dimnames = list("L", paste0("c", seq_along(values))))
      group <- factor(rep("A", length(values)), levels = "A")
      expected <- reference(values)
      direct <- summarize_cellchat_trimean(
        expression, group, genes = "L", normalize = FALSE
      )$average
      prepared <- prepare_sparse_trimean(expression, normalize = FALSE)
      permuted <- summarize_prepared_trimean(
        prepared, group, source_for_target = rev(seq_along(values))
      )$average
      expect_identical(as.numeric(direct), expected)
      expect_identical(as.numeric(permuted), expected)
    }
  }
})
