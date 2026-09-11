root <- "/home/dzf/cellchat_acceleration/releases/SpatialESS_compatibility_20260911"
prior <- "/home/dzf/cellchat_acceleration/validation/collaborator_cervical_20260911"
.libPaths(c("/home/dzf/cellchat_acceleration/validation/stratified_V3_patched_20260910/lib_patched", .libPaths()))
suppressPackageStartupMessages(library(SpatialCellChat))
input <- readRDS(file.path(prior, "inputs/chat_normalized.rds"))
records <- function(chat) {
  p <- chat@net[["prob", exact = TRUE]]
  q <- chat@net[["pval", exact = TRUE]]
  ix <- which(p > 0, arr.ind = TRUE)
  data.frame(sender_group_name = dimnames(p)[[1]][ix[, 1]],
    receiver_group_name = dimnames(p)[[2]][ix[, 2]], interaction_name = dimnames(p)[[3]][ix[, 3]],
    probability = p[ix], pvalue = q[ix])
}
fixture <- list(expression = input@data.signaling, coordinates = as.matrix(input@images$coordinates),
  group = input@idents, lr = input@LR$LRsig, db = input@DB,
  default = records(readRDS(file.path(prior, "results/official_normalized/chat_complete.rds"))),
  relaxed = records(readRDS(file.path(prior, "results/relaxed_normalized/chat_complete.rds"))),
  parameters = list(ratio = 1, tol = 5, interaction.range = 30, contact.range = 10,
    scale.distance = 1, nboot = 100L, seed.use = 1L),
  provenance = list(official_commit = "be88f300464cf970b36a6165ea7d2e41a47a7511",
    official_patch = "dims=c(NC,NC) only; percentage filter unchanged",
    dataset = "10x Xenium Prime FFPE Human Cervical Cancer",
    labels = "23 public graph clusters; not biological cell-type annotation"))
d <- file.path(root, "work/SpatialESS/tests/testthat/fixtures")
dir.create(d, showWarnings = FALSE)
saveRDS(fixture, file.path(d, "v3_percent_boundary.rds"))
mini <- file.path(root, "work/SpatialESS/validation/fixtures")
db_env <- new.env()
data("CellChatDB.human", package = "SpatialCellChat", envir = db_env)
public_fixture <- fixture
public_fixture$db <- db_env$CellChatDB.human
saveRDS(public_fixture, file.path(mini, "cervical_1000.rds"))
saveRDS(fixture$default, file.path(mini, "cervical_1000_expected_records.rds"))
write.table(fixture$default, gzfile(file.path(mini, "cervical_1000_expected_records.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
jsonlite::write_json(c(fixture$provenance, list(cells = 1000, expected_records = 57,
  nboot = 100, seed = 1, min_percent = .1, min_cells_sr = 5,
  expected_source = "frozen official V3 group probabilities and p-values, not ESS-generated")),
  file.path(mini, "cervical_source_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Saved official-anchored regression fixture", nrow(fixture$default), nrow(fixture$relaxed), "\n")
