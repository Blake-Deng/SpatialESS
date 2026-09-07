#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript audit_gse250346_lmm_tables.R NEW_OUTPUT BASELINE_OUTPUT")
out <- normalizePath(args[1L])
old <- normalizePath(args[2L])
script_path <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())])
source(file.path(dirname(normalizePath(script_path)), "read_legacy_lmm_tsv.R"))
stopifnot(file.exists(file.path(out, "COMPLETE")))
read_table <- function(name, root = out) {
  read.delim(file.path(root, name), check.names = FALSE, stringsAsFactors = FALSE)
}
write_table <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = TRUE)
}
lmm <- read_table("patient_random_intercept_lmm.tsv")
glm <- read_table("slice_level_glm.tsv")
baseline <- read_legacy_lmm_tsv(file.path(old, "patient_random_intercept_lmm.tsv"))
baseline_glm <- read_table("slice_level_glm.tsv", old)
stopifnot(nrow(lmm) == 50040L, !anyDuplicated(lmm$feature_id),
          setequal(lmm$feature_id, baseline$feature_id),
          identical(lmm$feature_id, glm$feature_id))
baseline <- baseline[match(lmm$feature_id, baseline$feature_id), ]
baseline_glm <- baseline_glm[match(glm$feature_id, baseline_glm$feature_id), ]
ok <- lmm$status == "ok"
stopifnot(all(is.finite(lmm$p_value[ok])), all(is.finite(lmm$estimate[ok])),
          all(is.finite(lmm$std_error[ok])), all(lmm$std_error[ok] > 0),
          all(is.na(lmm$q_value[!ok])),
          all(lmm$p_value[ok] >= 0 & lmm$p_value[ok] <= 1),
          all(lmm$zero_fraction >= 0 & lmm$zero_fraction <= 1))
expected <- rep(NA_real_, nrow(lmm))
expected[ok] <- p.adjust(lmm$p_value[ok], "BH")
stopifnot(isTRUE(all.equal(expected, lmm$q_value, tolerance = 1e-12)))
all_statuses <- c("ok", "convergence_warning", "numerical_failure", "fit_failed",
                  "insufficient_nonzero", "insufficient_zero_support",
                  "insufficient_variation", "rank_deficient", "insufficient_units")
stopifnot(all(lmm$status %in% all_statuses),
          all(lmm$total_units == 45L), all(lmm$patient_count == 35L))
failed <- lmm$status %in% c("numerical_failure", "fit_failed")
stopifnot(all(is.na(lmm$p_value[failed])),
          all(!is.na(lmm$fit_error[failed])),
          all(lmm$optimizer_attempts[failed] == "bobyqa,nloptwrap,Nelder_Mead"))
common <- ok & baseline$status == "ok"
sig <- ok & !is.na(lmm$q_value) & lmm$q_value < 0.05
old_sig <- baseline$status == "ok" & !is.na(baseline$q_value) & baseline$q_value < 0.05
write_table(data.frame(
  baseline_valid = sum(baseline$status == "ok"), robust_valid = sum(ok),
  baseline_warning = sum(baseline$status == "convergence_warning"),
  robust_warning = sum(lmm$status == "convergence_warning"),
  baseline_failed = sum(baseline$status == "fit_failed"),
  robust_failed = sum(failed),
  old_failed_now_ok = sum(baseline$status == "fit_failed" & ok),
  old_warning_now_ok = sum(baseline$status == "convergence_warning" & ok),
  baseline_significant = sum(old_sig), robust_significant = sum(sig),
  overlap_significant = sum(sig & old_sig),
  significant_jaccard = sum(sig & old_sig) / sum(sig | old_sig),
  common_valid = sum(common),
  estimate_spearman_common_valid = cor(lmm$estimate[common], baseline$estimate[common],
                                       method = "spearman"),
  max_abs_estimate_difference_common_valid = max(abs(lmm$estimate[common] -
                                                       baseline$estimate[common])),
  direction_agreement_common_valid = mean(sign(lmm$estimate[common]) ==
                                           sign(baseline$estimate[common]))
), "gse250346_lmm_before_after_summary.tsv")
non_singular <- ok & lmm$singular %in% FALSE
q_subset <- p.adjust(lmm$p_value[non_singular], "BH")
write_table(data.frame(
  valid_non_singular = sum(non_singular),
  significant_using_primary_global_BH = sum(sig & non_singular),
  significant_after_BH_recomputed_in_non_singular_subset = sum(q_subset < 0.05),
  scope = "Sensitivity: restrict to non-singular fits; BH universe is explicitly changed"
), "gse250346_lmm_non_singular_sensitivity.tsv")
support <- aggregate(cbind(zero_fraction, response_sd) ~ status, lmm,
                     FUN = median, na.action = na.pass, na.rm = TRUE)
write_table(support, "gse250346_lmm_support_diagnostics.tsv")
write_table(data.frame(
  check = c("complete_50040_unique_features", "all_45_slices_35_patients",
            "global_BH_matches", "non_ok_excluded_from_BH",
            "failed_have_no_p_value", "failed_keep_three_optimizer_attempts",
            "glm_status_matches_baseline", "glm_estimates_match_baseline"),
  passed = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
             identical(glm$status, baseline_glm$status),
             isTRUE(all.equal(glm$estimate, baseline_glm$estimate, tolerance = 1e-12)))
), "gse250346_lmm_result_audit.tsv")
print(read_table("gse250346_lmm_before_after_summary.tsv"))
print(read_table("gse250346_lmm_result_audit.tsv"))
stopifnot(all(read_table("gse250346_lmm_result_audit.tsv")$passed))
