test_that("multi-select sparse triMean covers zero and small-group boundaries", {
  set.seed(20260731)
  cases <- list()
  case_index <- 1L
  for (cell_count in 1:10) {
    for (positive_count in 0:cell_count) {
      positive <- if (positive_count) {
        sample(c(0.5, 1, 1, 2, 4), positive_count, replace = TRUE)
      } else {
        numeric()
      }
      values <- sample(c(rep(0, cell_count - positive_count), positive))
      cases[[case_index]] <- values
      case_index <- case_index + 1L
    }
  }
  cases <- c(
    cases,
    list(
      c(rep(0, 24), 8),
      c(rep(0, 25), 8),
      c(rep(0, 49), 3),
      c(rep(0, 50), 3),
      c(rep(0, 74), 2),
      c(rep(0, 75), 2),
      rep(3, 101)
    )
  )

  for (i in seq_along(cases)) {
    values <- cases[[i]]
    expression <- Matrix::Matrix(
      matrix(values, nrow = 1L,
             dimnames = list("gene", paste0("cell", seq_along(values)))),
      sparse = TRUE
    )
    observed <- summarize_cellchat_trimean(
      expression, rep("group", length(values)),
      genes = "gene", normalize = FALSE
    )$average[1L, 1L]
    expected <- mean(stats::quantile(
      values, c(0.25, 0.5, 0.5, 0.75), type = 7, names = FALSE
    ))
    expect_equal(observed, expected, tolerance = 1e-15,
                 info = sprintf("boundary case %d", i))
  }
})

test_that("sort and multi-select type-7 kernels are numerically identical", {
  set.seed(42)
  for (cell_count in c(8L, 31L, 100L, 1001L)) {
    for (zero_fraction in c(0, 0.24, 0.49, 0.74, 0.9)) {
      positive_count <- max(1L, round(cell_count * (1 - zero_fraction)))
      positive <- sample(c(0.25, 0.5, 1, 2, 4, 8),
                         positive_count, replace = TRUE)
      result <- SpatialESS:::benchmark_sparse_type7_cpp(
        positive, cell_count, iterations = 2L
      )
      expect_equal(result$select_value, result$sort_value, tolerance = 1e-15)
    }
  }
})

