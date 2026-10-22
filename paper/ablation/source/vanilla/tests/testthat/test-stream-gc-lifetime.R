test_that("percentage filtering preserves group sizes across R allocations", {
  n_cells <- 88L
  expression <- methods::as(Matrix::Matrix(matrix(1, n_cells, 2L), sparse = TRUE), "dgCMatrix")
  graph <- build_radius_graph_csr(cbind(seq_len(n_cells), 0), radius = 2,
                                  store_distance = TRUE)
  empty <- matrix(0L, 1L, 0L)
  args <- list(expression, graph$offsets, graph$neighbors, graph$distances,
               rep(seq_len(44L), 2L), 44L,
               matrix(1L), matrix(2L), empty, empty, empty, empty,
               FALSE, FALSE, FALSE,
               cbind(seq_len(n_cells), rev(seq_len(n_cells))),
               .5, 1, 1, 1, .1, 1L, TRUE, 1e8, 5e8)
  native <- getFromNamespace("spatialess_stream_cpp", "SpatialESS")
  expected <- do.call(native, args)
  old <- gctorture2(0L)
  on.exit(gctorture2(old), add = TRUE)
  for (iteration in seq_len(3L)) {
    gctorture2(1L)
    observed <- do.call(native, args)
    gctorture2(0L)
    expect_identical(observed, expected)
  }
})
