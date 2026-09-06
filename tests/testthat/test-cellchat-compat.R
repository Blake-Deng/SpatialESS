test_that("sparse triMean matches CellChat type-7 definition", {
  expression <- Matrix::Matrix(
    rbind(
      g1 = c(0, 1, 2, 8, 0, 0, 4, 12),
      g2 = c(1, 3, 5, 7, 2, 4, 6, 8)
    ),
    sparse = TRUE,
    dimnames = list(c("g1", "g2"), paste0("c", 1:8))
  )
  group <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
  observed <- summarize_cellchat_trimean(
    expression, group, genes = rownames(expression), normalize = FALSE
  )$average
  tri_mean <- function(x) mean(quantile(x, c(0.25, 0.5, 0.5, 0.75), type = 7))
  expected <- vapply(levels(group), function(level) {
    apply(as.matrix(expression[, group == level, drop = FALSE]), 1, tri_mean)
  }, numeric(nrow(expression)))
  expect_equal(observed, expected, tolerance = 1e-15)
})

test_that("CSR group support matches a materialized edge aggregation", {
  coords <- rbind(c1 = c(0, 0), c2 = c(1, 0), c3 = c(2, 0), c4 = c(3, 0))
  graph <- build_radius_graph_csr(
    coords, radius = 1.01, weight = "gaussian", scale = 0.8
  )
  group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  observed <- build_group_support_csr(graph, group)$group_pairs
  edges <- materialize_csr_edges(graph)
  expected <- aggregate(
    cbind(supported_edges = rep(1, nrow(edges)),
          sum_spatial_weight = edges$weight),
    by = list(sender_group = as.integer(group[edges$sender]),
              receiver_group = as.integer(group[edges$receiver])),
    FUN = sum
  )
  expected$mean_spatial_weight <- expected$sum_spatial_weight /
    expected$supported_edges
  columns <- c("sender_group", "receiver_group", "supported_edges",
               "sum_spatial_weight", "mean_spatial_weight")
  expected <- expected[order(expected$sender_group, expected$receiver_group), ]
  rownames(expected) <- NULL
  expect_equal(observed[, columns], expected[, columns], tolerance = 1e-15)
})

test_that("complex and cofactor group score follows CellChat formulas", {
  genes <- c("L1", "L2", "R1", "R2", "CA", "CI", "AGO", "ANT")
  values <- rbind(
    L1 = c(2, 3, 4, 5, 4, 5, 6, 7),
    L2 = c(1, 2, 3, 4, 3, 4, 5, 6),
    R1 = c(3, 4, 5, 6, 2, 3, 4, 5),
    R2 = c(2, 3, 4, 5, 1, 2, 3, 4),
    CA = c(1, 2, 1, 2, 2, 3, 2, 3),
    CI = c(2, 1, 2, 1, 1, 2, 1, 2),
    AGO = c(1, 1, 2, 2, 2, 2, 3, 3),
    ANT = c(2, 2, 1, 1, 3, 3, 2, 2)
  )
  expression <- Matrix::Matrix(
    values, sparse = TRUE,
    dimnames = list(genes, paste0("c", seq_len(ncol(values))))
  )
  group <- factor(rep(c("A", "B"), each = 4), levels = c("A", "B"))
  coords <- cbind(x = seq_len(ncol(values)), y = 0)
  rownames(coords) <- colnames(expression)
  graph <- build_radius_graph_csr(coords, radius = 20, weight = "binary")

  complex <- data.frame(subunit_1 = c("L1", "R1"),
                        subunit_2 = c("L2", "R2"),
                        row.names = c("LC", "RC"))
  cofactor <- data.frame(cofactor1 = c("CA", "CI", "AGO", "ANT"),
                         cofactor2 = "",
                         row.names = c("coa", "coi", "ago", "ant"))
  lr <- data.frame(
    interaction_name = "LC_RC", ligand = "LC", receptor = "RC",
    agonist = "ago", antagonist = "ant",
    co_A_receptor = "coa", co_I_receptor = "coi"
  )
  components <- prepare_cellchat_lr_components(
    lr, rownames(expression), complex, cofactor
  )
  averages <- summarize_cellchat_trimean(
    expression, group, components, normalize = FALSE
  )
  support <- build_group_support_csr(graph, group)
  observed <- score_cellchat_group_support(
    averages, support, components, Kh = 0.5, n = 1
  )$group_pairs

  avg <- averages$average
  ligand <- sqrt(avg["L1", ] * avg["L2", ])
  receptor <- sqrt(avg["R1", ] * avg["R2", ]) *
    (1 + avg["CA", ]) / (1 + avg["CI", ])
  agonist <- 1 + avg["AGO", ] / (0.5 + avg["AGO", ])
  antagonist <- 0.5 / (0.5 + avg["ANT", ])
  product <- outer(ligand, receptor)
  expected <- product / (0.5 + product) *
    outer(agonist, agonist) * outer(antagonist, antagonist)
  expected_at_support <- expected[cbind(observed$sender_group,
                                        observed$receiver_group)]
  expect_equal(observed$molecular_probability, expected_at_support,
               tolerance = 1e-14)
  expect_equal(observed$mean_weighted_probability,
               observed$molecular_probability, tolerance = 1e-15)
})
