root <- "/home/dzf/cellchat_acceleration/releases/SpatialESS_compatibility_20260911"
old <- "/home/dzf/cellchat_acceleration/validation/stratified_V3_patched_20260910"
.libPaths(c(file.path(root, "lib_new"), file.path(old, "lib_patched"), .libPaths()))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SpatialESS))
suppressPackageStartupMessages(library(SpatialCellChat))
input <- file.path(root, "validation/dual_track")
manifest <- jsonlite::fromJSON(file.path(input, "manifest.json"))
genes <- read.delim(file.path(input, "genes.tsv"))$gene
env <- new.env()
data("CellChatDB.human", package = "SpatialCellChat", envir = env)
db <- env$CellChatDB.human
coverage <- filter_cellchat_lr(db$interaction, genes, db$complex, db$cofactor)
full_db <- db
db$interaction <- coverage$lr
saveRDS(full_db, file.path(input, "database.rds"))
for (case in manifest$cases$case) {
  d <- file.path(input, case)
  cells <- read.delim(file.path(d, "metadata.tsv"))
  counts <- as(readMM(file.path(d, "counts.mtx")), "CsparseMatrix")
  dimnames(counts) <- list(genes, cells$cell_id)
  x <- counts %*% Diagonal(x = 10000/colSums(counts))
  x@x <- log1p(x@x)
  dimnames(x) <- dimnames(counts)
  coordinates <- as.matrix(cells[, c("x_centroid", "y_centroid")])
  rownames(coordinates) <- cells$cell_id
  group <- factor(cells$group)
  names(group) <- cells$cell_id
  chat <- createSpatialCellChat(x, meta = data.frame(labels = group, row.names = cells$cell_id),
    group.by = "labels", datatype = "spatial", coordinates = coordinates,
    spatial.factors = list(ratio = 1, tol = 5), do.sparse = TRUE)
  chat@DB <- db
  chat <- subsetData(chat)
  chat@LR$LRsig <- coverage$lr
  saveRDS(chat, file.path(d, "chat_input.rds"))
  parameters <- list(ratio = 1, tol = 5, interaction.range = 30, contact.range = 10,
    scale.distance = manifest$scale_distance, Kh = .5, n = 1, use.AGAN = TRUE,
    contact.dependent = TRUE, contact.dependent.forced = FALSE, nboot = 100L, seed.use = 1L)
  f <- list(expression = chat@data.signaling, coordinates = coordinates, group = group,
    lr = coverage$lr, db = db, full_db = full_db, parameters = parameters)
  saveRDS(f, file.path(d, "fixture.rds"))
  saveRDS(list(cells = colnames(f$expression), genes = rownames(f$expression),
    expression = f$expression, coordinates = coordinates, group = group,
    lr = f$lr, db = db, parameters = parameters), file.path(d, "input_contract.rds"))
  write.table(coverage$audit, file.path(d, "lr_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(case, "prepared", ncol(f$expression), "cells", nrow(f$expression), "genes", nrow(f$lr), "LR\n")
}
writeLines(capture.output(sessionInfo()), file.path(input, "sessionInfo.txt"))
