#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SpatialESS)
})

root <- Sys.getenv(
  "GSE250346_MULTISAMPLE_ROOT",
  "/data/dzf/GSE250346/benchmarks/multisample"
)
input_dir <- file.path(root, "full_sample_results")
metadata_path <- file.path(root, "sample_metadata_from_rds.tsv")
out_dir <- Sys.getenv(
  "GSE250346_LMM_OUT",
  file.path(root, "slice_glm_vs_patient_lmm")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.delim(metadata_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
results <- lapply(files, readRDS)
names(results) <- sub("\\.rds$", "", basename(files))
metadata <- metadata[match(names(results), metadata$sample), , drop = FALSE]
if (anyNA(metadata$sample) || anyNA(metadata$patient)) {
  stop("Sample metadata does not cover all result files.", call. = FALSE)
}

prepared <- prepare_multisample_communication(
  results = results,
  sample_metadata = data.frame(
    sample_id = metadata$sample,
    patient_id = metadata$patient,
    condition = metadata$disease_status,
    sample_affect = metadata$sample_affect,
    tma = metadata$tma,
    run = metadata$run,
    stringsAsFactors = FALSE
  ),
  sample_col = "sample_id",
  patient_col = "patient_id",
  condition_col = "condition",
  covariates = c("sample_affect", "tma", "run"),
  unit = "sample"
)

design <- ~ condition
coefficient <- "conditionDisease"
start <- proc.time()[["elapsed"]]
slice_glm <- fit_multisample_communication_glm(
  prepared, design = design, coefficient = coefficient,
  min_units = 6L, min_nonzero = 2L
)
lmm <- fit_multisample_communication_lmm(
  prepared, design = design, coefficient = coefficient,
  patient_col = "patient_id", min_units = 6L, min_nonzero = 2L,
  REML = FALSE
)
elapsed <- proc.time()[["elapsed"]] - start

write.table(slice_glm, file.path(out_dir, "slice_level_glm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lmm, file.path(out_dir, "patient_random_intercept_lmm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
status_summary <- as.data.frame(table(lmm$status), stringsAsFactors = FALSE)
colnames(status_summary) <- c("status", "features")
write.table(status_summary, file.path(out_dir, "lmm_status_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
failure_messages <- lmm$fit_error[!is.na(lmm$fit_error) & nzchar(lmm$fit_error)]
if (length(failure_messages)) {
  failure_counts <- sort(table(failure_messages), decreasing = TRUE)
  failure_summary <- data.frame(
    fit_error = names(failure_counts),
    features = as.integer(failure_counts),
    stringsAsFactors = FALSE
  )
} else {
  failure_summary <- data.frame(
    fit_error = character(), features = integer(),
    stringsAsFactors = FALSE
  )
}
write.table(failure_summary, file.path(out_dir, "lmm_failure_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
convergence_messages <- lmm$convergence_message[
  !is.na(lmm$convergence_message) & nzchar(lmm$convergence_message)
]
if (length(convergence_messages)) {
  convergence_class <- vapply(convergence_messages, function(message) {
    if (grepl("negative eigen", message, fixed = TRUE)) {
      "negative_hessian_eigenvalue"
    } else if (grepl("numerically singular", message, fixed = TRUE)) {
      "numerically_singular_hessian"
    } else if (grepl("max|grad|", message, fixed = TRUE) &
               grepl("very large eigenvalue", message, fixed = TRUE)) {
      "gradient_and_scale_warning"
    } else if (grepl("max|grad|", message, fixed = TRUE)) {
      "gradient_warning"
    } else if (grepl("very large eigenvalue", message, fixed = TRUE)) {
      "scale_warning"
    } else {
      "other"
    }
  }, character(1L))
  convergence_counts <- sort(table(convergence_class), decreasing = TRUE)
  convergence_summary <- data.frame(
    convergence_class = names(convergence_counts),
    features = as.integer(convergence_counts),
    stringsAsFactors = FALSE
  )
} else {
  convergence_summary <- data.frame(
    convergence_class = character(), features = integer(),
    stringsAsFactors = FALSE
  )
}
write.table(
  convergence_summary,
  file.path(out_dir, "lmm_convergence_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    cohort = "GSE250346",
    samples = nrow(prepared$score),
    patients = length(unique(prepared$unit_metadata$patient_id)),
    control_samples = sum(prepared$unit_metadata$condition == "Control"),
    disease_samples = sum(prepared$unit_metadata$condition == "Disease"),
    model_glm = "log1p(score) ~ condition",
    model_lmm = "log1p(score) ~ condition + (1 | patient)",
    coefficient = coefficient,
    elapsed_seconds = elapsed,
    glm_ok = sum(slice_glm$status == "ok"),
    lmm_ok = sum(lmm$status == "ok"),
    lmm_convergence_warning = sum(lmm$status == "convergence_warning"),
    lmm_numerical_failure = sum(lmm$status == "numerical_failure"),
    lmm_fit_failed = sum(lmm$status == "fit_failed"),
    lmm_rank_deficient = sum(lmm$status == "rank_deficient"),
    lmm_insufficient_zero_support = sum(
      lmm$status == "insufficient_zero_support"
    ),
    lmm_insufficient_variation = sum(
      lmm$status == "insufficient_variation"
    ),
    max_zero_fraction = 0.98,
    min_response_sd = 1e-12,
    rank_tolerance = 1e-10,
    optimizer_sequence = "bobyqa,nloptwrap,Nelder_Mead",
    lmm_singular = sum(lmm$status == "ok" & lmm$singular %in% TRUE,
                       na.rm = TRUE),
    lmm_non_singular = sum(lmm$status == "ok" & !lmm$singular %in% TRUE,
                           na.rm = TRUE),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "model_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
print(read.delim(file.path(out_dir, "model_summary.tsv")))
