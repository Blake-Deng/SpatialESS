single_boundary_indices <- function(ligand, receptor, genes,
                                     complex = data.frame(
                                       subunit_1 = character(),
                                       row.names = character())) {
  cofactor <- data.frame(cofactor1 = character(), row.names = character())
  lr <- data.frame(interaction_name = "test", ligand = ligand, receptor = receptor)
  components <- prepare_cellchat_lr_components(lr, genes, complex, cofactor)
  list(components = components, indices = lapply(components[c(
    "ligand", "receptor", "co_a", "co_i", "agonist", "antagonist"
  )], .spatialess_index_matrix, genes = genes))
}

single_boundary_native_scores <- function(average, fixture, Kh = .5, n = 1) {
  indices <- fixture$indices
  components <- fixture$components
  extra <- list(components$has_agonist, components$has_antagonist, Kh, n)
  support <- do.call(score_cellchat_group_support_cpp, c(
    list(average, 1L, 1L, 1, 1, 1), indices, extra, list(100)
  ))$group_pairs$molecular_probability
  matrix_score <- as.numeric(do.call(score_cellchat_group_matrix_cpp, c(
    list(average, 1L, 1L), indices, extra
  )))
  record <- do.call(score_cellchat_group_records_cpp, c(
    list(average, 1L, 1L, 1L), indices, extra
  ))
  c(support = support, matrix = matrix_score, records = record)
}

test_that("all group scorers preserve a true single-gene probability tie", {
  genes <- c("L", "R")
  fixture <- single_boundary_indices("L", "R", genes)
  score <- function(l, r) single_boundary_native_scores(
    matrix(c(l, r), 2L, 1L, dimnames = list(genes, "A")), fixture
  )
  expected <- (.2 * .5) / (.5 + .2 * .5)
  observed <- score(.2, .5)
  bootstrap <- score(.1, 1)
  expect_identical(unname(observed), rep(expected, 3L))
  expect_identical(bootstrap, observed)
  expect_identical(unname(bootstrap > observed), rep(FALSE, 3L))
  # Sparse support omits a zero score; dense and record scorers retain zero.
  zero_scores <- score(0, 1)
  expect_identical(names(zero_scores), c("matrix", "records"))
  expect_identical(unname(zero_scores), rep(0, 2L))
})

test_that("named one-subunit complexes retain the official geometric mean", {
  genes <- c("L", "R")
  complex <- data.frame(subunit_1 = c("L", "R"),
                        row.names = c("LC1", "L"))
  fixture <- single_boundary_indices("LC1", "R", genes, complex)
  expect_identical(attr(fixture$indices$ligand, "spatialess_single_gene"), FALSE)
  expect_identical(attr(fixture$indices$receptor, "spatialess_single_gene"), TRUE)
  average <- matrix(c(.1, 1), 2L, 1L, dimnames = list(genes, "A"))
  ligand <- exp(mean(log(.1)))
  expected <- ligand / (.5 + ligand)
  scores <- single_boundary_native_scores(average, fixture)
  expect_identical(unname(scores), rep(expected, 3L))

  # The upstream computeExpr_LR first checks expression rownames, even when the
  # database contains a complex with the same name and different subunits.
  direct <- single_boundary_indices("L", "R", genes, complex)
  expect_identical(direct$components$ligand[[1]], "L")
  expect_identical(attr(direct$indices$ligand, "spatialess_single_gene"), TRUE)
  direct_expected <- .1 / (.5 + .1)
  direct_scores <- single_boundary_native_scores(average, direct)
  expect_identical(unname(direct_scores), rep(direct_expected, 3L))
  expect_false(identical(scores, direct_scores))
})

test_that("old component objects are recovered only when molecular kind is known", {
  genes <- c("L", "R")
  direct <- single_boundary_indices("L", "R", genes)$components
  attr(direct$ligand, "spatialess_single_gene") <- NULL
  attr(direct$receptor, "spatialess_single_gene") <- NULL
  recovered <- .spatialess_component_index_matrices(direct, genes)
  expect_identical(attr(recovered$ligand, "spatialess_single_gene"), TRUE)
  expect_identical(attr(recovered$receptor, "spatialess_single_gene"), TRUE)

  complex <- data.frame(subunit_1 = "L", row.names = "LC1")
  ambiguous <- single_boundary_indices("LC1", "R", genes, complex)$components
  attr(ambiguous$ligand, "spatialess_single_gene") <- NULL
  attr(ambiguous$receptor, "spatialess_single_gene") <- NULL
  expect_error(.spatialess_component_index_matrices(ambiguous, genes),
               "re-run prepare_cellchat_lr_components")
})

test_that("complex scorers use R mean rounding for all subunit lengths", {
  set.seed(418)
  for (width in c(2L, 3L, 4L, 7L)) {
    genes <- c(paste0("L", seq_len(width)), "R")
    complex <- as.data.frame(as.list(setNames(
      genes[seq_len(width)], paste0("subunit_", seq_len(width))
    )), stringsAsFactors = FALSE)
    rownames(complex) <- "LC"
    fixture <- single_boundary_indices("LC", "R", genes, complex)
    for (iteration in seq_len(30L)) {
      ligand <- 10^runif(width, -12, 0)
      average <- matrix(c(ligand, 1), width + 1L, 1L,
                        dimnames = list(genes, "A"))
      upstream_ligand <- exp(mean(log(ligand)))
      expected <- upstream_ligand / (.5 + upstream_ligand)
      expect_identical(unname(single_boundary_native_scores(average, fixture)),
                       rep(expected, 3L))
    }
  }
})

test_that("cofactor and regulator arithmetic follows the reference expression order", {
  genes <- c(paste0("L", 1:3), paste0("R", 1:2), paste0("CA", 1:3),
             paste0("CI", 1:3), paste0("AG", 1:4), paste0("AN", 1:4))
  complex <- data.frame(subunit_1 = c("L1", "R1"),
                        subunit_2 = c("L2", "R2"), subunit_3 = c("L3", ""),
                        row.names = c("LC", "RC"))
  cofactor <- data.frame(cofactor1 = c("CA1", "CI1", "AG1", "AN1"),
                         cofactor2 = c("CA2", "CI2", "AG2", "AN2"),
                         cofactor3 = c("CA3", "CI3", "AG3", "AN3"),
                         cofactor4 = c("", "", "AG4", "AN4"),
                         row.names = c("CA", "CI", "AG", "AN"))
  lr <- data.frame(interaction_name = "LC_RC", ligand = "LC", receptor = "RC",
                   co_A_receptor = "CA", co_I_receptor = "CI",
                   agonist = "AG", antagonist = "AN")
  components <- prepare_cellchat_lr_components(lr, genes, complex, cofactor)
  fixture <- list(components = components,
    indices = .spatialess_component_index_matrices(components, genes))
  set.seed(607)
  for (power in c(1, 2, .7, 3.9)) {
    for (iteration in seq_len(30L)) {
      values <- setNames(10^runif(length(genes), -5, 0), genes)
      average <- matrix(values, length(genes), 1L,
                        dimnames = list(genes, "A"))
      ligand <- exp(mean(log(values[paste0("L", 1:3)])))
      receptor <- exp(mean(log(values[paste0("R", 1:2)])))
      receptor <- receptor * prod(1 + values[paste0("CA", 1:3)]) /
        prod(1 + values[paste0("CI", 1:3)])
      Kh_n <- .5^power
      product_n <- (ligand * receptor)^power
      expected <- product_n / (Kh_n + product_n)
      ag <- values[paste0("AG", 1:4)]^power
      agonist <- prod(1 + ag / (Kh_n + ag))
      an <- values[paste0("AN", 1:4)]^power
      antagonist <- prod(Kh_n / (Kh_n + an))
      expected <- expected * (agonist * agonist)
      expected <- expected * (antagonist * antagonist)
      expect_identical(unname(single_boundary_native_scores(
        average, fixture, n = power)), rep(unname(expected), 3L))
    }
  }
})

test_that("single-gene block-null tail counts match the direct upstream rule", {
  expression <- Matrix::Matrix(
    rbind(L = c(.2, .2, .1, .1), R = c(.5, .5, 1, 1)), sparse = TRUE,
    dimnames = list(c("L", "R"), paste0("c", 1:4))
  )
  group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  coordinates <- cbind(x = 1:4, y = 0)
  rownames(coordinates) <- colnames(expression)
  graph <- build_radius_graph_csr(coordinates, radius = 10, weight = "binary")
  support <- build_group_support_csr(graph, group)
  blocks <- build_spatial_blocks(coordinates, block_size = 100, origin = c(0, 0))
  components <- single_boundary_indices("L", "R", rownames(expression))$components
  prepared <- prepare_sparse_trimean(expression, genes = components$genes,
                                     normalize = FALSE)
  observed <- summarize_prepared_trimean(prepared, group)$average
  direct_score <- function(avg) {
    product <- avg["L", "A"] * avg["R", "A"]
    product / (.5 + product)
  }
  observed_direct <- direct_score(observed)
  direct_reject <- vapply(seq_len(100L), function(index) {
    permutation <- permute_within_spatial_blocks(blocks, seed = index)
    average <- summarize_prepared_trimean(
      prepared, group, source_for_target = permutation
    )$average
    direct_score(average) > observed_direct
  }, logical(1))
  fit <- permutation_cellchat_group_support(
    prepared, group, support, components, blocks, nperm = 100L, seed = 1L,
    finite_correction = FALSE, tail = "strict_greater", verbose = FALSE
  )$group_pairs
  row <- fit[fit$sender_group == 1L & fit$receiver_group == 1L, ]
  expect_equal(nrow(row), 1L)
  expect_identical(as.integer(row$n_reject), as.integer(sum(direct_reject)))
  expect_identical(row$pvalue, sum(direct_reject) / 100)
  expect_identical(row$probability, observed_direct)
  expect_identical(as.integer(row$n_reject), 67L)
  expect_identical(row$pvalue, .67)
})
