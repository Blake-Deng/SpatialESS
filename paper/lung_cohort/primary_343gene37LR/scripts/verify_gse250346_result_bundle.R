#!/usr/bin/env Rscript
# Verify published feature tables with base R and without the expression data.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript verify_gse250346_result_bundle.R TABLE_ROOT")
root <- normalizePath(args[1L])
new <- file.path(root, "gse250346_lmm_robust_20260907")
old <- file.path(root, "gse250346_lmm_baseline")
read_gz <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
}
lmm <- read_gz(file.path(new, "patient_random_intercept_lmm.tsv.gz"))
glm <- read_gz(file.path(new, "slice_level_glm.tsv.gz"))
baseline <- read_gz(file.path(old, "patient_random_intercept_lmm.tsv.gz"))
baseline_glm <- read_gz(file.path(old, "slice_level_glm.tsv.gz"))
meta <- read.delim(file.path(new, "sample_metadata.tsv"))
stopifnot(nrow(lmm) == 50040L, !anyDuplicated(lmm$feature_id),
          identical(lmm$feature_id, glm$feature_id),
          setequal(lmm$feature_id, baseline$feature_id),
          setequal(glm$feature_id, baseline_glm$feature_id))
baseline <- baseline[match(lmm$feature_id, baseline$feature_id), ]
baseline_glm <- baseline_glm[match(glm$feature_id, baseline_glm$feature_id), ]
ok <- lmm$status == "ok"
failed <- lmm$status == "numerical_failure"
p <- lmm$p_value
p[!ok] <- NA_real_
stopifnot(nrow(meta) == 45L, length(unique(meta$patient_id)) == 35L,
          all(lmm$total_units == 45L), all(lmm$patient_count == 35L),
          sum(ok) == 39570L, sum(failed) == 292L,
          sum(lmm$status == "convergence_warning") == 1984L,
          sum(lmm$status == "insufficient_nonzero") == 8194L,
          all(is.finite(p[ok])), all(p[ok] >= 0 & p[ok] <= 1),
          all(is.finite(lmm$estimate[ok])), all(is.finite(lmm$std_error[ok])),
          all(lmm$std_error[ok] > 0),
          all(is.na(lmm$q_value[!ok])), all(is.na(lmm$p_value[failed])),
          all(lmm$optimizer_attempts[failed] == "bobyqa,nloptwrap,Nelder_Mead"),
          isTRUE(all.equal(p.adjust(p, "BH"), lmm$q_value, tolerance = 1e-12)),
          sum(lmm$q_value < 0.05, na.rm = TRUE) == 5147L,
          identical(glm$status, baseline_glm$status),
          identical(glm$estimate, baseline_glm$estimate))
common <- ok & baseline$status == "ok"
stopifnot(sum(common) == 39242L,
          max(abs(lmm$estimate[common] - baseline$estimate[common])) == 0,
          sum(baseline$status == "convergence_warning" & ok) == 328L,
          sum(baseline$status == "fit_failed" & ok) == 0L)
audit <- read.delim(file.path(root, "gse250346_lmm_result_audit.tsv"))
stopifnot(nrow(audit) == 8L, all(audit$passed))
print(table(lmm$status))
cat("PASS: 50,040 features; complete cohort; global BH; exclusions; baseline agreement.\n")
