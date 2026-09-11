coverage_fixture <- function() {
  list(
    genes = c("L", "R", "R2", "ZERO"),
    complex = data.frame(subunit_1 = c("R", "R", ""),
                         subunit_2 = c("R2", "missing", ""),
                         row.names = c("complete", "partial", "empty")),
    cofactor = data.frame(cofactor1 = "missing_cofactor", row.names = "optional"),
    lr = data.frame(interaction_name = paste0("lr", 1:6),
                    ligand = c("L", "L", "L", "absent", "L", "ZERO"),
                    receptor = c("R", "complete", "partial", "R", "empty", "R"),
                    agonist = "optional", row.names = paste0("row", 1:6))
  )
}

test_that("coverage requires every subunit and reports exact excluded genes", {
  f <- coverage_fixture()
  expect_warning(a <- filter_cellchat_lr(f$lr, f$genes, f$complex, f$cofactor), "Excluded 3 of 6")
  expect_identical(a$valid, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(a$lr, f$lr[c(1, 2, 6), , drop = FALSE])
  expect_identical(a$excluded, f$lr[c(3, 4, 5), , drop = FALSE])
  expect_identical(a$audit$missing_receptor_genes[3], "missing")
  expect_identical(a$audit$reason[5], "unresolved_receptor_definition")
  expect_setequal(a$missing_genes$gene[a$missing_genes$required], c("missing", "absent"))
  expect_true(all(!a$missing_genes$required[a$missing_genes$role == "agonist"]))
  # A measured zero-expression row is retained; there is no abundance filter.
  expect_true(a$valid[6])
})

test_that("empty input, missing identifiers and malformed tables are diagnosed", {
  f <- coverage_fixture()
  expect_error(filter_cellchat_lr(f$lr, c("L", "L"), f$complex), "unique")
  expect_error(filter_cellchat_lr(f$lr, c("L", ""), f$complex), "non-empty")
  expect_error(filter_cellchat_lr(f$lr, f$genes, data.frame(a = "R")), "subunit")
  blank <- filter_cellchat_lr(f$lr[FALSE, ], f$genes, f$complex)
  expect_length(blank$valid, 0)
  expect_equal(nrow(blank$missing_genes), 0)
  malformed <- f$lr[1:2, ]; malformed$ligand <- c(NA, "")
  expect_warning(a <- filter_cellchat_lr(malformed, f$genes, f$complex), "retained 0")
  expect_true(all(a$audit$reason == "unresolved_ligand_definition"))
  expect_error(prepare_cellchat_lr_components(f$lr, f$genes, f$complex, f$cofactor),
               "filter_cellchat_lr")
})

test_that("explicit and in-call filters preserve records, RNG and full-matrix normalization", {
  f <- coverage_fixture()
  # The unused row sets the normalization maximum: filtering must not remove it.
  expr <- Matrix::Matrix(rbind(L = c(1, 2, 3, 2, 0, 1), R = c(3, 2, 1, 2, 1, 0),
                              R2 = c(1, 1, 2, 0, 1, 3), ZERO = rep(0, 6), UNUSED = rep(100, 6)), sparse = TRUE)
  colnames(expr) <- paste0("c", 1:6)
  coords <- cbind(x = 1:6, y = 0); rownames(coords) <- colnames(expr)
  g <- factor(rep(c("A", "B"), each = 3))
  f$lr$annotation <- rep(c("Cell-Cell Contact", "Secreted Signaling"), 3)
  args <- list(expression = expr, coordinates = coords, group = g,
               complex = f$complex, cofactor = f$cofactor,
               interaction.range = 10, contact.range = 5, tol = 0,
               min.percent = 0, min.cells.sr = 1, nboot = 19, seed.use = 7)
  expect_error(do.call(spatialess, c(args, list(lr = f$lr))), "missing genes")
  before <- do.call(spatialess, c(args, list(lr = f$lr[c(1, 2, 6), ])))
  rng <- .Random.seed
  expect_warning(after <- do.call(spatialess, c(args, list(lr = f$lr, lr_missing = "filter"))), "Excluded 3")
  expect_identical(after$records, before$records)
  expect_identical(after$graph, before$graph)
  expect_identical(.Random.seed, rng)
  expect_equal(nrow(after$lr_filter$excluded), 3)
  expect_identical(after$lr_filter$audit$input_row, 1:6)
  expect_error(suppressWarnings(do.call(spatialess, c(args, list(lr = f$lr[3:5, ], lr_missing = "filter")))), "No complete LR")
})

test_that("factor columns and custom component column names retain input order", {
  f <- coverage_fixture(); lr <- f$lr[c(2, 1), ]
  lr$ligand <- factor(lr$ligand)
  names(lr)[names(lr) == "ligand"] <- "L"
  names(lr)[names(lr) == "receptor"] <- "R"
  a <- filter_cellchat_lr(lr, f$genes, f$complex, ligand_col = "L", receptor_col = "R")
  expect_identical(a$lr, lr)
  expect_true(all(a$valid))
})
