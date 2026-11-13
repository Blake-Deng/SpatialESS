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

test_that("measured LR rows win name collisions without changing cofactor resolution", {
  genes <- c("L", "R")
  complex <- data.frame(subunit_1 = c("missing_L", "missing_R"),
                         row.names = c("L", "R"))
  cofactor <- data.frame(cofactor1 = "missing_cofactor", row.names = "L")
  lr <- data.frame(ligand = "L", receptor = "R", agonist = "L")

  audit <- filter_cellchat_lr(lr, genes, complex, cofactor)
  components <- prepare_cellchat_lr_components(lr, genes, complex, cofactor)
  expect_identical(audit$valid, TRUE)
  expect_identical(audit$lr, lr)
  expect_identical(audit$audit$ligand_genes, "L")
  expect_identical(audit$audit$receptor_genes, "R")
  expect_identical(audit$audit$missing_ligand_genes, "")
  expect_identical(audit$audit$missing_receptor_genes, "")
  expect_identical(audit$missing_genes$gene, "missing_cofactor")
  expect_identical(audit$missing_genes$role, "agonist")
  expect_identical(audit$missing_genes$required, FALSE)
  expect_identical(components$ligand[[1L]], "L")
  expect_identical(components$receptor[[1L]], "R")
  expect_identical(attr(components$ligand, "spatialess_single_gene"), TRUE)
  expect_identical(attr(components$receptor, "spatialess_single_gene"), TRUE)
  expect_length(components$agonist[[1L]], 0L)
})

test_that("in-call LR filtering retains direct genes that collide with incomplete complexes", {
  expression <- Matrix::Matrix(rbind(L=c(.2,.2,.1,.1),R=c(.5,.5,1,1)),sparse=TRUE,
    dimnames=list(c("L","R"),paste0("c",1:4)))
  coordinates <- cbind(x=1:4,y=0);rownames(coordinates)<-colnames(expression)
  lr <- data.frame(ligand="L",receptor="R",interaction_name="L_R")
  collision <- data.frame(subunit_1=c("missing_L","missing_R"),row.names=c("L","R"))
  empty <- data.frame(subunit_1=character(),row.names=character())
  cofactor <- data.frame(cofactor1=character(),row.names=character())
  args <- list(expression=expression,coordinates=coordinates,
    group=factor(c("A","A","B","B")),lr=lr,cofactor=cofactor,
    tol=0,interaction.range=10,scale.distance=1,min.percent=0,min.cells.sr=1,
    nboot=7L,seed.use=19L)
  filtered <- do.call(spatialess,c(args,list(complex=collision,lr_missing="filter")))
  direct <- do.call(spatialess,c(args,list(complex=empty)))
  expect_identical(filtered$records,direct$records)
  expect_identical(filtered$lr_filter$valid,TRUE)
  expect_equal(nrow(filtered$lr_filter$excluded),0)
})
