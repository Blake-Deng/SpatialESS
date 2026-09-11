root <- "/home/dzf/cellchat_acceleration/releases/SpatialESS_compatibility_20260911"
.libPaths(c("/home/dzf/cellchat_acceleration/validation/stratified_V3_patched_20260910/lib_patched", .libPaths()))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SpatialCellChat))
future::plan(future::sequential)
e <- new.env(parent = asNamespace("SpatialCellChat"))
source(file.path(root, "reference/original/R/spatial.R"), local = e)
stopifnot(!grepl("dims = c(NC, NC)", paste(deparse(body(e$computeCellDistance)), collapse = "\n"), fixed = TRUE))
out <- file.path(root, "validation/original_distance")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rows <- list()
for (case in c("center_1000", "center_5000", "stratified_1000", "stratified_5000")) {
  f <- readRDS(file.path(root, "validation/dual_track", case, "fixture.rds"))
  a <- e$computeCellDistance(f$coordinates, interaction.range = 30, contact.range = 10, ratio = 1, tol = 5)
  b <- SpatialCellChat::computeCellDistance(f$coordinates, interaction.range = 30, contact.range = 10, ratio = 1, tol = 5)
  aa <- summary(a$d.spatial)
  bb <- summary(b$d.spatial)
  nc <- ncol(f$expression)
  rows[[case]] <- data.frame(case = case, cells = nc, original_rows = nrow(a$d.spatial),
    original_cols = ncol(a$d.spatial), patched_rows = nrow(b$d.spatial),
    distance_entries_identical = identical(aa$i, bb$i) && identical(aa$j, bb$j) && identical(aa$x, bb$x),
    original_complete_shape = identical(dim(a$d.spatial), c(nc, nc)),
    patched_complete_shape = identical(dim(b$d.spatial), c(nc, nc)))
}
write.table(do.call(rbind, rows), file.path(out, "dimensions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(do.call(rbind, rows))
# Execute the unmodified distance path with all 1076 candidates, not a synthetic multiplication.
original_compute <- SpatialCellChat::computeCommunProb
environment(original_compute) <- e
chat <- readRDS(file.path(root, "validation/dual_track/stratified_1000/chat_input.rds"))
options(error = function() {
  message <- geterrmessage()
  writeLines(message, file.path(out, "original_full_call_error.txt"))
  writeLines(if (grepl("non-conformable matrix dimensions", message, fixed = TRUE))
    "expected_original_dimension_failure" else "unexpected_failure", file.path(out, "status.txt"))
  traceback(20)
  quit(status = 1L, save = "no")
})
chat <- original_compute(chat, LR.use = chat@LR$LRsig, tol = 5, interaction.range = 30,
                         contact.range = 10, scale.distance = 1)
writeLines("unexpected_original_success", file.path(out, "status.txt"))
