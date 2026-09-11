#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialCellChat)
  library(SpatialESS)
})

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

complex <- data.frame(
  subunit_1 = c("L1", "R1"), subunit_2 = c("L2", "R2"),
  row.names = c("LC", "RC")
)
cofactor <- data.frame(
  cofactor1 = c("CA", "CI", "AGO", "ANT"), cofactor2 = "",
  row.names = c("coa", "coi", "ago", "ant")
)
lr <- data.frame(
  interaction_name = "LC_RC", ligand = "LC", receptor = "RC",
  agonist = "ago", antagonist = "ant",
  co_A_receptor = "coa", co_I_receptor = "coi"
)

components <- prepare_cellchat_lr_components(
  lr, rownames(expression), complex, cofactor
)
engine_result <- spatialess(
  expression, coords, group, lr, complex, cofactor,
  ratio = 1, tol = 0, interaction.range = 20, contact.range = 1,
  scale.distance = 1, min.percent = 0, min.cells.sr = 1, nboot = 0
)
stopifnot(inherits(engine_result, "SpatialESSResult"), nrow(engine_result$records) > 0)

average <- summarize_cellchat_trimean(
  expression, group, components, normalize = TRUE
)
graph <- build_radius_graph_csr(coords, radius = 20, weight = "binary")
support <- build_group_support_csr(graph, group)
observed <- score_cellchat_group_support(
  average, support, components, Kh = 0.5, n = 1
)$group_pairs

normalized <- expression / max(expression@x)
reference_average <- aggregate(
  t(as.matrix(normalized)), list(group), FUN = SpatialCellChat::triMean
)
reference_average <- t(reference_average[, -1, drop = FALSE])
colnames(reference_average) <- levels(group)
rownames(reference_average) <- rownames(expression)

ligand <- SpatialCellChat::computeExpr_LR(lr$ligand, reference_average, complex)
receptor <- SpatialCellChat::computeExpr_LR(lr$receptor, reference_average, complex)
co_a <- SpatialCellChat::computeExpr_coreceptor(
  cofactor, reference_average, lr, type = "A"
)
co_i <- SpatialCellChat::computeExpr_coreceptor(
  cofactor, reference_average, lr, type = "I"
)
receptor <- receptor * co_a / co_i
product <- Matrix::crossprod(matrix(ligand, nrow = 1),
                             matrix(receptor, nrow = 1))
p1 <- product / (0.5 + product)
agonist <- SpatialCellChat::computeExpr_agonist(
  reference_average, lr, cofactor, index.agonist = 1, Kh = 0.5, n = 1
)
antagonist <- SpatialCellChat::computeExpr_antagonist(
  reference_average, lr, cofactor, index.antagonist = 1, Kh = 0.5, n = 1
)
reference_probability <- p1 *
  Matrix::crossprod(matrix(agonist, nrow = 1)) *
  Matrix::crossprod(matrix(antagonist, nrow = 1))
expected <- reference_probability[cbind(observed$sender_group,
                                        observed$receiver_group)]

average_difference <- max(abs(as.matrix(summarize_cellchat_trimean(
  expression, group, components, normalize = TRUE)$average) - reference_average))
probability_difference <- max(abs(observed$molecular_probability - expected))
stopifnot(average_difference < 1e-14, probability_difference < 1e-14)

cat(sprintf("CellChat triMean max difference: %.3g\n", average_difference))
cat(sprintf("CellChat complex/cofactor probability max difference: %.3g\n",
            probability_difference))
cat(sprintf("Spatially supported group pairs checked: %d\n", nrow(observed)))
out <- Sys.getenv("SPATIALESS_MINI_OUTPUT", "spatialess_mini_records.tsv")
utils::write.table(observed, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote %s\n", normalizePath(out, mustWork = FALSE)))
