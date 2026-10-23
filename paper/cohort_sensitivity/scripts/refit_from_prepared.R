#!/usr/bin/env Rscript
# Refit using the supplied saved matrix without raw Seurat inputs or LR inference.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: refit_from_prepared.R BUNDLE NEW_OUTPUT CORES")
bundle <- normalizePath(args[1]); out <- args[2]; cores <- as.integer(args[3])
stopifnot(cores > 0L, !file.exists(out))
dir.create(out, recursive = TRUE); out <- normalizePath(out)
for (sub in c("scripts", "data", "logs", "reports", "figures")) dir.create(file.path(out, sub))
file.copy(list.files(file.path(bundle, "scripts"), full.names = TRUE), file.path(out, "scripts"))
for (name in c("prepared_slice_scores.rds", "primary_lmm.tsv.gz", "slice_scores.tsv.gz", "sample_metadata.tsv", "feature_mapping.csv")) {
  stopifnot(file.copy(file.path(bundle, "data", name), file.path(out, "data", name)))
}
file.copy(list.files(file.path(bundle, "reports"), pattern = "^(condition_by_|design_rank_audit|input_validation)", full.names = TRUE), file.path(out, "reports"))
engine <- new.env(parent = globalenv())
sys.source(file.path(bundle, "scripts/multisample_source.R"), envir = engine)
prepared <- readRDS(file.path(bundle, "data/prepared_slice_scores.rds"))
prepared$unit_metadata$condition <- factor(prepared$unit_metadata$condition, levels = c("Control", "Disease"))
prepared$unit_metadata$run <- factor(prepared$unit_metadata$run)
prepared$unit_metadata$tma <- factor(prepared$unit_metadata$tma)
tasks <- split(seq_len(ncol(prepared$score)), ceiling(seq_len(ncol(prepared$score))/100L))
for (model in c("run", "tma")) {
  started <- proc.time()[["elapsed"]]
  design <- if (model == "run") ~ condition + run else ~ condition + tma
  fits <- parallel::mclapply(tasks, function(idx) {
    part <- prepared; part$score <- prepared$score[, idx, drop = FALSE]
    part$feature_metadata <- prepared$feature_metadata[idx, , drop = FALSE]
    engine$fit_multisample_communication_lmm(part, design, "conditionDisease", patient_col = "patient_id")
  }, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
  stopifnot(!any(vapply(fits, inherits, logical(1), "try-error")))
  fit <- do.call(rbind, fits); stopifnot(identical(fit$feature_id, colnames(prepared$score)))
  p <- fit$p_value; p[fit$status != "ok"] <- NA_real_
  fit$q_value <- p.adjust(p, "BH"); fit$ci_lower <- fit$estimate - qnorm(.975)*fit$std_error
  fit$ci_upper <- fit$estimate + qnorm(.975)*fit$std_error
  data.table::fwrite(fit, file.path(out, "data", paste0(model,"_adjusted_lmm.tsv.gz")), sep = "\t", na = "NA")
  cat(model, "finished in", proc.time()[["elapsed"]]-started, "seconds\n")
}
writeLines(capture.output(sessionInfo()), file.path(out, "logs/R_session_info.txt"))
cat("Next run run_patient_sensitivity.py, verify_robust_layer.R, then summarize_and_plot.py using this new output.\n")
