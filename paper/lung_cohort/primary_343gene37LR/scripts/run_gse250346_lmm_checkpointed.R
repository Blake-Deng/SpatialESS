#!/usr/bin/env Rscript
# Feature fits are independent; BH correction is recomputed after all chunks.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript run_gse250346_lmm_checkpointed.R REPO INPUT_ROOT OUTPUT CORES")
}
repo <- normalizePath(args[1L])
root <- normalizePath(args[2L])
out <- args[3L]
cores <- as.integer(args[4L])
stopifnot(is.finite(cores), cores >= 1L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
dir.create(file.path(out, "chunks"), showWarnings = FALSE)
engine <- new.env(parent = globalenv())
source_file <- file.path(repo, "R", "multisample.R")
sys.source(source_file, envir = engine)
source(file.path(repo, "benchmarks", "read_legacy_lmm_tsv.R"))
stopifnot(requireNamespace("lme4", quietly = TRUE))
write_tsv <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", row.names = FALSE,
              quote = TRUE, na = "NA")
}
log_event <- function(...) {
  cat(format(Sys.time(), "%FT%T%z"), ..., "\n")
  flush.console()
}
input_files <- sort(list.files(file.path(root, "full_sample_results"),
                              pattern = "\\.rds$", full.names = TRUE))
meta_file <- file.path(root, "sample_metadata_from_rds.tsv")
stopifnot(length(input_files) == 45L, file.exists(meta_file))
input_md5 <- tools::md5sum(c(source_file, meta_file, input_files))
signature <- list(md5 = unname(input_md5), names = basename(names(input_md5)),
                  R = R.version.string, lme4 = as.character(packageVersion("lme4")),
                  chunk_size = 50L)
signature_file <- file.path(out, "checkpoint_signature.rds")
if (file.exists(signature_file)) {
  stopifnot(identical(readRDS(signature_file), signature))
} else {
  saveRDS(signature, signature_file)
}
write_tsv(data.frame(file = signature$names, md5 = signature$md5),
          "input_and_source_checksums.tsv")
file.copy(source_file, file.path(out, "multisample_source.R"), overwrite = TRUE)
log_event("Loading", length(input_files), "sample profiles")
metadata <- read.delim(meta_file, stringsAsFactors = FALSE)
results <- lapply(input_files, readRDS)
names(results) <- sub("\\.rds$", "", basename(input_files))
metadata <- metadata[match(names(results), metadata$sample), , drop = FALSE]
stopifnot(!anyNA(metadata$patient), !anyNA(metadata$sample))
prepared <- engine$prepare_multisample_communication(
  results,
  data.frame(sample_id = metadata$sample, patient_id = metadata$patient,
             condition = metadata$disease_status, sample_affect = metadata$sample_affect,
             tma = metadata$tma, run = metadata$run),
  unit = "sample", covariates = c("sample_affect", "tma", "run")
)
rm(results)
stopifnot(nrow(prepared$score) == 45L, ncol(prepared$score) == 50040L)
write_tsv(prepared$unit_metadata, "sample_metadata.tsv")
tasks <- split(seq_len(ncol(prepared$score)),
               ceiling(seq_len(ncol(prepared$score)) / signature$chunk_size))
checkpoint <- function(id) file.path(out, "chunks", sprintf("chunk_%04d.rds", id))
pending <- which(!vapply(seq_along(tasks), function(i) file.exists(checkpoint(i)),
                        logical(1L)))
log_event("Features", ncol(prepared$score), "pending chunks", length(pending),
          "workers", cores)
started <- proc.time()[["elapsed"]]
done <- parallel::mclapply(pending, function(id) {
  selected <- tasks[[id]]
  part <- prepared
  part$score <- prepared$score[, selected, drop = FALSE]
  part$feature_metadata <- prepared$feature_metadata[selected, , drop = FALSE]
  lmm <- engine$fit_multisample_communication_lmm(
    part, ~ condition, "conditionDisease", patient_col = "patient_id",
    min_units = 6L, min_nonzero = 2L, REML = FALSE,
    max_zero_fraction = 0.98, min_response_sd = 1e-12, rank_tolerance = 1e-10,
    optimizer = "bobyqa", optimizer_fallbacks = c("nloptwrap", "Nelder_Mead")
  )
  glm <- engine$fit_multisample_communication_glm(
    part, ~ condition, "conditionDisease", min_units = 6L, min_nonzero = 2L
  )
  tmp <- paste0(checkpoint(id), ".tmp")
  saveRDS(list(lmm = lmm, glm = glm), tmp)
  stopifnot(file.rename(tmp, checkpoint(id)))
  log_event("Completed chunk", id, "features", length(selected))
  id
}, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
stopifnot(!any(vapply(done, inherits, logical(1L), "try-error")),
          all(file.exists(vapply(seq_along(tasks), checkpoint, character(1L)))))
parts <- lapply(seq_along(tasks), function(i) readRDS(checkpoint(i)))
lmm <- do.call(rbind, lapply(parts, "[[", "lmm"))
glm <- do.call(rbind, lapply(parts, "[[", "glm"))
rownames(lmm) <- rownames(glm) <- NULL
stopifnot(identical(lmm$feature_id, colnames(prepared$score)),
          identical(glm$feature_id, lmm$feature_id), !anyDuplicated(lmm$feature_id))
for (kind in c("lmm", "glm")) {
  model <- get(kind)
  p <- model$p_value
  p[model$status != "ok"] <- NA_real_
  model$q_value <- p.adjust(p, "BH")
  stopifnot(all(is.na(model$q_value[model$status != "ok"])))
  assign(kind, model)
}
write_tsv(lmm, "patient_random_intercept_lmm.tsv")
write_tsv(glm, "slice_level_glm.tsv")
status <- as.data.frame(table(lmm$status), stringsAsFactors = FALSE)
names(status) <- c("status", "features")
write_tsv(status, "gse250346_lmm_fit_status.tsv")
count_messages <- function(values, column, name) {
  values <- values[!is.na(values) & nzchar(values)]
  tab <- as.data.frame(table(values), stringsAsFactors = FALSE)
  names(tab) <- c(column, "features")
  write_tsv(tab, name)
}
count_messages(lmm$fit_error, "fit_error", "gse250346_lmm_failure_summary.tsv")
classify <- function(message) {
  if (is.na(message) || !nzchar(message)) return(NA_character_)
  if (grepl("negative eigen", message, fixed = TRUE)) return("negative_hessian_eigenvalue")
  if (grepl("numerically singular", message, fixed = TRUE)) return("numerically_singular_hessian")
  grad <- grepl("max|grad|", message, fixed = TRUE)
  scale <- grepl("very large eigenvalue", message, fixed = TRUE)
  if (grad && scale) return("gradient_and_scale_warning")
  if (grad) return("gradient_warning")
  if (scale) return("scale_warning")
  "other"
}
count_messages(vapply(lmm$convergence_message, classify, character(1L)),
               "convergence_class", "gse250346_lmm_convergence_summary.tsv")
count_messages(lmm$optimizer_attempts, "optimizer_attempts",
               "gse250346_lmm_optimizer_summary.tsv")
summary <- data.frame(
  cohort = "GSE250346", samples = nrow(prepared$score),
  patients = length(unique(metadata$patient)),
  control_samples = sum(metadata$disease_status == "Control"),
  disease_samples = sum(metadata$disease_status == "Disease"),
  model_glm = "log1p(score) ~ condition",
  model_lmm = "log1p(score) ~ condition + (1 | patient)",
  coefficient = "conditionDisease", features = nrow(lmm),
  elapsed_seconds = if (length(pending) == length(tasks)) {
    proc.time()[["elapsed"]] - started
  } else NA_real_,
  invocation_elapsed_seconds = proc.time()[["elapsed"]] - started,
  chunks_fitted_this_invocation = length(pending),
  total_chunks = length(tasks),
  elapsed_scope = "Excludes input loading; elapsed_seconds is NA on checkpoint resume",
  workers = cores, native_threads_per_worker = 1L,
  glm_ok = sum(glm$status == "ok"), lmm_ok = sum(lmm$status == "ok"),
  lmm_convergence_warning = sum(lmm$status == "convergence_warning"),
  lmm_numerical_failure = sum(lmm$status == "numerical_failure"),
  lmm_fit_failed = sum(lmm$status == "fit_failed"),
  lmm_rank_deficient = sum(lmm$status == "rank_deficient"),
  lmm_insufficient_nonzero = sum(lmm$status == "insufficient_nonzero"),
  lmm_insufficient_zero_support = sum(lmm$status == "insufficient_zero_support"),
  lmm_insufficient_variation = sum(lmm$status == "insufficient_variation"),
  lmm_singular = sum(lmm$status == "ok" & lmm$singular %in% TRUE),
  lmm_non_singular = sum(lmm$status == "ok" & lmm$singular %in% FALSE),
  max_zero_fraction = 0.98, min_response_sd = 1e-12, rank_tolerance = 1e-10,
  optimizer_sequence = "bobyqa,nloptwrap,Nelder_Mead",
  p_method = "Wald_normal", fdr_method = "BH_over_all_status_ok",
  source_md5 = signature$md5[1L]
)
write_tsv(summary, "gse250346_lmm_model_summary.tsv")
sig_g <- glm$status == "ok" & !is.na(glm$q_value) & glm$q_value < 0.05
sig_l <- lmm$status == "ok" & !is.na(lmm$q_value) & lmm$q_value < 0.05
common <- glm$status == "ok" & lmm$status == "ok"
write_tsv(data.frame(
  glm_significant = sum(sig_g), lmm_significant = sum(sig_l),
  overlap_significant = sum(sig_g & sig_l),
  significant_jaccard = sum(sig_g & sig_l) / sum(sig_g | sig_l),
  estimate_spearman_ok = cor(glm$estimate[common], lmm$estimate[common], method = "spearman"),
  direction_agreement_ok = mean(sign(glm$estimate[common]) == sign(lmm$estimate[common])),
  paired_valid_features = sum(common), lmm_ok = sum(lmm$status == "ok"),
  lmm_non_singular = sum(lmm$status == "ok" & lmm$singular %in% FALSE),
  lmm_significant_non_singular = sum(sig_l & lmm$singular %in% FALSE)
), "gse250346_glm_vs_lmm_summary.tsv")
baseline <- read_legacy_lmm_tsv(file.path(root, "slice_glm_vs_patient_lmm",
                                        "patient_random_intercept_lmm.tsv"))
stopifnot(setequal(baseline$feature_id, lmm$feature_id))
baseline <- baseline[match(lmm$feature_id, baseline$feature_id), , drop = FALSE]
transitions <- as.data.frame(table(baseline_status = baseline$status,
                                   robust_status = lmm$status))
names(transitions)[3L] <- "features"
write_tsv(transitions[transitions$features > 0, ], "gse250346_lmm_status_transitions.tsv")
write_tsv(data.frame(feature_id = lmm$feature_id, baseline_status = baseline$status,
                    robust_status = lmm$status, baseline_q = baseline$q_value,
                    robust_q = lmm$q_value),
          "gse250346_lmm_feature_transitions.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
print(summary)
log_event("COMPLETE: all features present; global BH recalculated; no excluded fit has q-value")
writeLines("complete", file.path(out, "COMPLETE"))
