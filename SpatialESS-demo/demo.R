library(Matrix)
library(SpatialESS)
expression <- Matrix(rbind(L = c(1, 2, 1, 0, 2, 1), R = c(2, 1, 0, 2, 1, 3)), sparse = TRUE)
colnames(expression) <- paste0("cell", 1:6)
coordinates <- cbind(x = 1:6, y = 0)
rownames(coordinates) <- colnames(expression)
group <- factor(rep(c("A", "B"), each = 3))
db <- list(
  interaction = data.frame(interaction_name = c("L_R", "L_RC"),
                            ligand = "L", receptor = c("R", "RC")),
  complex = data.frame(subunit_1 = "R", subunit_2 = "unmeasured", row.names = "RC"),
  cofactor = data.frame(cofactor1 = character(), row.names = character())
)
# A warning is expected: L_RC cannot be measured with this panel.
coverage <- filter_cellchat_lr(db$interaction, rownames(expression), db$complex, db$cofactor)
print(coverage$audit)
print(coverage$missing_genes)
result <- spatialess(expression, coordinates, group, coverage$lr, db$complex, db$cofactor,
                     tol = 0, interaction.range = 5, min.cells.sr = 1, min.percent = 0,
                     nboot = 20, seed.use = 1)
stopifnot(nrow(result$records) > 0, nrow(coverage$lr) == 1)
print(result$records)
