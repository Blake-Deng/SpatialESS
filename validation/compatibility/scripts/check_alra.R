root <- "/home/dzf/cellchat_acceleration/releases/SpatialESS_compatibility_20260911"
prior <- "/home/dzf/cellchat_acceleration/validation/collaborator_cervical_20260911"
.libPaths(c(file.path(root, "lib_new"),
  "/home/dzf/cellchat_acceleration/validation/stratified_V3_patched_20260910/lib_patched", .libPaths()))
suppressPackageStartupMessages(library(SpatialESS))
suppressPackageStartupMessages(library(SpatialCellChat))
chat <- readRDS(file.path(prior, "inputs/chat_alra_meringue.rds"))
official <- readRDS(file.path(prior, "results/official_alra_meringue/chat_complete.rds"))
stopifnot(identical(chat@data.signaling, official@data.signaling),
  identical(chat@idents, official@idents), identical(chat@LR$LRsig, official@LR$LRsig))
ess <- spatialess(chat@data.signaling, as.matrix(chat@images$coordinates), chat@idents,
  lr = chat@LR$LRsig, complex = chat@DB$complex, cofactor = chat@DB$cofactor,
  ratio = 1, tol = 5, interaction.range = 30, contact.range = 10, scale.distance = 1,
  min.percent = .1, min.cells.sr = 5L, nboot = 100L, seed.use = 1L)
p <- official@net[["prob", exact = TRUE]]; q <- official@net[["pval", exact = TRUE]]
ix <- which(p > 0, arr.ind = TRUE)
a <- data.frame(sender_group_name = dimnames(p)[[1]][ix[, 1]],
  receiver_group_name = dimnames(p)[[2]][ix[, 2]], interaction_name = dimnames(p)[[3]][ix[, 3]],
  probability = p[ix], pvalue = q[ix])
key <- function(x) paste(x$interaction_name, x$sender_group_name, x$receiver_group_name, sep = "|")
b <- ess$records[match(key(a), key(ess$records)), ]
summary <- data.frame(preprocessing = "shared_ALRA_MERINGUE", cells = ncol(chat@data.signaling),
  lr = nrow(chat@LR$LRsig), official_active = nrow(a), ess_active = nrow(ess$records),
  active_jaccard = length(intersect(key(a), key(ess$records)))/length(union(key(a), key(ess$records))),
  max_abs_probability_difference = max(abs(a$probability - b$probability)),
  pvalue_exact_fraction = mean(a$pvalue == b$pvalue),
  significance_disagreements = sum((a$pvalue < .05) != (b$pvalue < .05)))
summary$pass <- with(summary, active_jaccard == 1 && max_abs_probability_difference <= 1e-15 &&
  pvalue_exact_fraction == 1 && significance_disagreements == 0)
out <- file.path(root, "validation/alra")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(ess, file.path(out, "ess.rds"))
write.table(summary, file.path(out, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
print(summary)
stopifnot(summary$pass)
