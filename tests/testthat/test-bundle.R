test_that("bundle reader preserves dimensions and cell order", {
  path <- tempfile("bundle_")
  dir.create(path)
  con <- gzfile(file.path(path, "expression.mtx.gz"), "wt")
  writeLines(
    c(
      "%%MatrixMarket matrix coordinate real general",
      "%",
      "2 2 2",
      "1 1 1",
      "2 2 2"
    ),
    con
  )
  close(con)
  writeLines(c("L", "R"), file.path(path, "genes.tsv"))
  writeLines(c("c1", "c2"), file.path(path, "cells.tsv"))
  utils::write.table(
    data.frame(
      cell = c("c1", "c2"), group = c("A", "B"),
      coordinate_1 = c(0, 1), coordinate_2 = c(0, 1)
    ),
    file.path(path, "metadata.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE
  )
  bundle <- read_spatialess_bundle(path)
  expect_s4_class(bundle$expression, "dgCMatrix")
  expect_identical(dim(bundle$expression), c(2L, 2L))
  expect_identical(colnames(bundle$expression), c("c1", "c2"))
  expect_identical(names(bundle$group), c("c1", "c2"))
})
