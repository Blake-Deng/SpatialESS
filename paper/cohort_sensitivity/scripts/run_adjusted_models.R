#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: run_adjusted_models.R REPO INPUT_ROOT OUTPUT CORES")
repo <- normalizePath(args[1]); input <- normalizePath(args[2]); out <- normalizePath(args[3])
cores <- as.integer(args[4]); stopifnot(cores >= 1L)
library(data.table)
setDTthreads(1L)
set.seed(20261012)
engine <- new.env(parent = globalenv())
source_file <- file.path(repo, "R/multisample.R")
sys.source(source_file, envir = engine)
dir.create(file.path(out, "data"), showWarnings = FALSE)
write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "NA", quote = "auto")
log_event <- function(...) { cat(format(Sys.time(), "%FT%T%z"), ..., "\n"); flush.console() }
files <- sort(list.files(file.path(input, "full_sample_results"), pattern = "\\.rds$", full.names = TRUE))
meta_file <- file.path(input, "sample_metadata_from_rds.tsv")
stopifnot(length(files) == 45L)
primary_file <- file.path(repo, "results/lmm/patient_random_intercept_lmm.tsv.gz")
primary <- read.delim(gzfile(primary_file), check.names = FALSE)
md5 <- tools::md5sum(c(files, meta_file, source_file, primary_file))
write_tsv(data.frame(file = names(md5), md5 = unname(md5)), file.path(out, "data/input_checksums.tsv"))
file.copy(source_file, file.path(out, "scripts/multisample_source.R"), overwrite = TRUE)
signature <- list(md5 = md5, R = R.version.string, lme4 = as.character(packageVersion("lme4")),
                  models = c("run", "tma"), chunk_size = 100L)
sig_path <- file.path(out, "data/checkpoint_signature.rds")
if (file.exists(sig_path)) stopifnot(identical(signature, readRDS(sig_path))) else saveRDS(signature, sig_path)
prepared_path <- file.path(out, "data/prepared_slice_scores.rds")
if (file.exists(prepared_path)) {
  prepared <- readRDS(prepared_path)
} else {
  log_event("Loading original 45 slice profiles")
  results <- lapply(files, readRDS); names(results) <- sub("\\.rds$", "", basename(files))
  meta <- read.delim(meta_file)
  prepared <- engine$prepare_multisample_communication(results,
    data.frame(sample_id = meta$sample, patient_id = meta$patient, condition = meta$disease_status,
               sample_affect = meta$sample_affect, run = meta$run, tma = meta$tma),
    unit = "sample", covariates = c("sample_affect", "run", "tma"))
  prepared$score <- prepared$score[, match(primary$feature_id, colnames(prepared$score)), drop = FALSE]
  prepared$feature_metadata <- prepared$feature_metadata[match(primary$feature_id, prepared$feature_metadata$feature_id), ]
  stopifnot(identical(colnames(prepared$score), primary$feature_id), !anyNA(prepared$score),
            identical(dim(prepared$score), c(45L, 50040L)))
  saveRDS(prepared, prepared_path)
  rm(results); gc()
}
meta <- prepared$unit_metadata
meta$condition <- factor(meta$condition, levels = c("Control", "Disease"))
meta$run <- factor(meta$run); meta$tma <- factor(meta$tma); meta$sample_affect <- factor(meta$sample_affect)
prepared$unit_metadata <- meta
sd_response <- apply(log1p(prepared$score), 2, sd)
valid_sd <- !is.na(primary$response_sd)
checks <- data.frame(check = c("45_slices", "50040_unique_features", "nonzero_support_matches_primary", "response_sd_matches_primary"),
  passed = c(nrow(prepared$score) == 45L, !anyDuplicated(primary$feature_id) && ncol(prepared$score) == 50040L,
             all(colSums(prepared$score > 0) == primary$nonzero_units),
             max(abs(sd_response[valid_sd] - primary$response_sd[valid_sd])) < 1e-12))
write_tsv(checks, file.path(out, "reports/input_validation.tsv")); stopifnot(all(checks$passed))
write_tsv(meta, file.path(out, "data/sample_metadata.tsv"))
write_tsv(prepared$feature_metadata, file.path(out, "data/feature_metadata.tsv"))
write_tsv(data.frame(sample_id = rownames(prepared$score), prepared$score, check.names = FALSE),
          file.path(out, "data/slice_scores.tsv.gz"))
file.copy(primary_file, file.path(out, "data/primary_lmm.tsv.gz"), overwrite = TRUE)
for (field in c("run", "tma", "sample_affect")) {
  write_tsv(as.data.frame(table(condition = meta$condition, category = meta[[field]])),
            file.path(out, "reports", paste0("condition_by_", field, ".tsv")))
}
designs <- list(primary = ~ condition, run = ~ condition + run, tma = ~ condition + tma,
                affect = ~ condition + sample_affect, redundant_run_tma = ~ condition + run + tma)
rank_audit <- do.call(rbind, lapply(names(designs), function(name) {
  x <- model.matrix(designs[[name]], meta)
  data.frame(design = name, formula = paste(deparse(designs[[name]]), collapse = ""),
             rank = qr(x, tol = 1e-10)$rank, columns = ncol(x), full_rank = qr(x, tol = 1e-10)$rank == ncol(x))
}))
write_tsv(rank_audit, file.path(out, "reports/design_rank_audit.tsv"))
stopifnot(all(rank_audit$full_rank[rank_audit$design %in% c("primary", "run", "tma")]),
          !rank_audit$full_rank[rank_audit$design == "affect"])
tasks <- split(seq_len(ncol(prepared$score)), ceiling(seq_len(ncol(prepared$score)) / 100L))
for (model in c("run", "tma")) {
  chunk_dir <- file.path(out, "data", paste0(model, "_chunks")); dir.create(chunk_dir, showWarnings = FALSE)
  checkpoint <- function(id) file.path(chunk_dir, sprintf("chunk_%04d.rds", id))
  pending <- which(!vapply(seq_along(tasks), function(i) file.exists(checkpoint(i)), logical(1)))
  log_event("Starting", model, "pending chunks", length(pending), "workers", cores)
  started <- proc.time()[["elapsed"]]
  done <- parallel::mclapply(pending, function(id) {
    idx <- tasks[[id]]; part <- prepared
    part$score <- prepared$score[, idx, drop = FALSE]
    part$feature_metadata <- prepared$feature_metadata[idx, , drop = FALSE]
    fit <- engine$fit_multisample_communication_lmm(part, designs[[model]], "conditionDisease",
      patient_col = "patient_id", min_units = 6L, min_nonzero = 2L, REML = FALSE,
      max_zero_fraction = 0.98, min_response_sd = 1e-12, rank_tolerance = 1e-10,
      optimizer = "bobyqa", optimizer_fallbacks = c("nloptwrap", "Nelder_Mead"))
    tmp <- paste0(checkpoint(id), ".tmp"); saveRDS(fit, tmp); stopifnot(file.rename(tmp, checkpoint(id)))
    if (id %% 25L == 0L) log_event("Completed", model, "chunk", id)
    id
  }, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
  stopifnot(!any(vapply(done, inherits, logical(1L), "try-error")),
            all(file.exists(vapply(seq_along(tasks), checkpoint, character(1)))))
  fit <- do.call(rbind, lapply(seq_along(tasks), function(i) readRDS(checkpoint(i))))
  stopifnot(identical(fit$feature_id, primary$feature_id))
  p <- fit$p_value; p[fit$status != "ok"] <- NA_real_
  fit$q_value <- p.adjust(p, "BH"); fit$ci_lower <- fit$estimate - qnorm(.975) * fit$std_error
  fit$ci_upper <- fit$estimate + qnorm(.975) * fit$std_error
  stopifnot(all(is.na(fit$q_value[fit$status != "ok"])))
  write_tsv(fit, file.path(out, "data", paste0(model, "_adjusted_lmm.tsv.gz")))
  write_tsv(data.frame(model = model, formula = paste(deparse(designs[[model]]), collapse = ""),
                       elapsed_seconds_this_invocation = proc.time()[["elapsed"]] - started,
                       workers = cores, ok = sum(fit$status == "ok"),
                       singular_ok = sum(fit$status == "ok" & fit$singular %in% TRUE),
                       hits = sum(fit$status == "ok" & fit$q_value < .05, na.rm = TRUE)),
            file.path(out, "reports", paste0(model, "_run_summary.tsv")))
  log_event("Completed full model", model, "BH hits", sum(fit$q_value < .05, na.rm = TRUE))
}
writeLines(capture.output(sessionInfo()), file.path(out, "logs/R_session_info.txt"))
writeLines("complete", file.path(out, "logs/ADJUSTED_MODELS_COMPLETE"))
log_event("COMPLETE")
