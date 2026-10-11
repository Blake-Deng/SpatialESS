#!/usr/bin/env Rscript
# Recheck saved records only; no inference, permutations or model fitting.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L))
root <- normalizePath(args[[1L]])
archive <- file.path(root, "paper", "provenance", "ovarian")
output <- if (length(args) == 2L) args[[2L]] else
  file.path(root, "paper", "provenance", "verification")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(digits = 17)

key <- function(x) {
  k <- paste(x$sender_group_name, x$receiver_group_name,
             x$interaction_name, sep = "\r")
  stopifnot(!anyDuplicated(k))
  k
}
compare <- function(a, b, case, comparison) {
  ka <- key(a); kb <- key(b)
  stopifnot(setequal(ka, kb))
  b <- b[match(ka, kb), , drop = FALSE]
  pd <- abs(a$probability - b$probability)
  pv <- abs(a$pvalue - b$pvalue)
  stopifnot(all(is.finite(pd)), all(is.finite(pv)))
  data.frame(case = case, comparison = comparison,
             reference_active = nrow(a), candidate_active = nrow(b),
             active_jaccard = 1,
             probability_pearson = cor(a$probability, b$probability),
             max_abs_probability_difference = max(pd),
             pvalue_exact_fraction = mean(pv == 0),
             max_abs_pvalue_difference = max(pv),
             significant_reference_lt005 = sum(a$pvalue < 0.05),
             significant_candidate_lt005 = sum(b$pvalue < 0.05),
             significance_disagreements = sum((a$pvalue < 0.05) != (b$pvalue < 0.05)))
}

rows <- list()
for (case in c("raw_1000", "raw_5000", "raw_16247",
               "sparkle_1000", "sparkle_5000", "sparkle_16247")) {
  folder <- file.path(archive, "bridge_regression", "runs", case)
  a <- readRDS(file.path(folder, "old_bridge", "records.rds"))
  b <- readRDS(file.path(folder, "new_bridge", "records.rds"))
  c <- readRDS(file.path(folder, "direct_main", "records.rds"))
  rows[[length(rows) + 1L]] <- compare(a, b, case, "old_vs_new")
  rows[[length(rows) + 1L]] <- compare(b, c, case, "bridge_vs_direct")
}
bridge <- do.call(rbind, rows)
stopifnot(all(bridge$max_abs_probability_difference == 0),
          all(bridge$pvalue_exact_fraction == 1),
          all(bridge$significance_disagreements == 0))
write.table(bridge, file.path(output, "bridge_12_rechecked.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

official <- readRDS(file.path(archive, "official_pair", "official", "records.rds"))
ess <- readRDS(file.path(archive, "official_pair", "spatialess", "records.rds"))
historical <- compare(official, ess, "sparkle_16247", "archived_official_pair")
assert_ovarian <- function(x, expected_error) {
  stopifnot(x$reference_active == 37808L, x$candidate_active == 37808L,
            x$significant_reference_lt005 == 8417L,
            x$significant_candidate_lt005 == 8417L,
            x$pvalue_exact_fraction == 1,
            x$significance_disagreements == 0,
            abs(x$max_abs_probability_difference / expected_error - 1) < 1e-12)
}
assert_ovarian(historical, 4.98732999343332e-18)
write.table(historical, file.path(output, "ovarian_archived_pair_rechecked.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

current <- readRDS(file.path(root, "results", "records", "ovarian_sparkle_16247.rds"))
current_check <- compare(official, current, "ovarian_sparkle_16247", "current_0.1.4_vs_archived_official")
assert_ovarian(current_check, 1.32923186346723e-16)
write.table(current_check, file.path(output, "ovarian_current_release_rechecked.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

tab <- read.delim(gzfile(file.path(root, "paper", "figure_package", "source_data",
                                  "archived", "sparkle_probability_records.tsv.gz")),
                  check.names = FALSE)
ko <- key(official); kt <- key(tab)
stopifnot(setequal(ko, kt))
tab <- tab[match(ko, kt), , drop = FALSE]
ke <- key(ess); ess <- ess[match(ko, ke), , drop = FALSE]
# The archived TSV used 15 significant decimal digits; verify that loss
# separately from the original RDS-to-RDS fidelity calculation.
stopifnot(all(abs(tab$official_probability - official$probability) <=
              5e-14 * abs(official$probability)),
          all(abs(tab$spatialess_probability - ess$probability) <=
              5e-14 * abs(ess$probability)),
          identical(tab$official_pvalue, official$pvalue),
          identical(tab$spatialess_pvalue, ess$pvalue))
serialization <- data.frame(
  field = c("official_probability", "spatialess_probability"),
  max_absolute_serialization_difference = c(max(abs(tab$official_probability - official$probability)),
                                           max(abs(tab$spatialess_probability - ess$probability))),
  max_relative_serialization_difference = c(max(abs(tab$official_probability - official$probability) / abs(official$probability)),
                                           max(abs(tab$spatialess_probability - ess$probability) / abs(ess$probability))),
  tolerance_relative = 5e-14)
write.table(serialization, file.path(output, "archived_tsv_serialization_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "verification_sessionInfo.txt"))
print(bridge)
print(historical)
print(current_check)
cat("PASS: 12 archived bridge comparisons, archived ovarian pair, current-release linkage, figure source records.\n")
