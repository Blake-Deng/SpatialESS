#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: summarize_lr_interface_regression.R VALIDATION_ROOT")
root <- normalizePath(args[1]); out <- file.path(root, "work", "SpatialESS", "validation", "lr_interface")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cases <- list.dirs(file.path(root, "runs"), recursive = FALSE)
key <- function(x) paste(x$interaction_name, x$sender_group_name, x$receiver_group_name, sep = "|")
compare <- function(a, b) {
  ka <- key(a); kb <- key(b)
  stopifnot(!anyDuplicated(ka), !anyDuplicated(kb))
  common <- intersect(ka, kb); i <- match(common, ka); j <- match(common, kb)
  sig_a <- ka[a$pvalue <= 0.05]; sig_b <- kb[b$pvalue <= 0.05]
  data.frame(active_a = nrow(a), active_b = nrow(b),
    active_jaccard = if(length(union(ka,kb))) length(common) / length(union(ka,kb)) else 1,
    max_abs_probability_difference = if(length(common)) max(abs(a$probability[i] - b$probability[j])) else NA_real_,
    pvalue_exact_fraction = if(length(common)) mean(a$pvalue[i] == b$pvalue[j]) else NA_real_,
    significant_a = length(sig_a), significant_b = length(sig_b),
    significance_disagreements = length(setdiff(sig_a, sig_b)) + length(setdiff(sig_b, sig_a)))
}
resources <- list(); comparisons <- list(); official <- list()
for (directory in cases) {
  olddir <- file.path(directory, "old_manual")
  old <- readRDS(file.path(olddir, "records.rds"))
  for (mode in c("old_manual", "new_explicit", "new_filter")) {
    d <- file.path(directory, mode)
    s <- read.delim(file.path(d, "summary.tsv")); log <- readLines(file.path(d, "time_verbose.txt"))
    s$peak_rss_gib <- as.numeric(sub(".*: ", "", grep("Maximum resident set size", log, value = TRUE))) / 1024^2
    elapsed <- strsplit(sub(".*\\): ", "", grep("Elapsed ", log, value = TRUE)), ":")[[1]]
    s$process_wall_seconds <- sum(rev(as.numeric(elapsed)) * 60^(seq_along(elapsed) - 1))
    s$exit_code <- as.integer(sub(".*: ", "", grep("Exit status", log, value = TRUE)))
    resources[[length(resources) + 1L]] <- s
    if (mode == "old_manual") next
    new <- readRDS(file.path(d, "records.rds"))
    c <- cbind(s[, c("case", "cells", "mode", "lr", "nboot", "seed")], compare(old, new))
    c$records_identical <- identical(old, new)
    c$input_contract_identical <- identical(readRDS(file.path(olddir, "input_contract.rds")), readRDS(file.path(d, "input_contract.rds")))
    c$result_contract_identical <- identical(readRDS(file.path(olddir, "result_contract.rds")), readRDS(file.path(d, "result_contract.rds")))
    c$rng_identical <- identical(readRDS(file.path(olddir, "rng.rds")), readRDS(file.path(d, "rng.rds")))
    c$pass <- c$records_identical & c$input_contract_identical & c$result_contract_identical & c$rng_identical &
      c$active_jaccard == 1 & c$max_abs_probability_difference == 0 & c$pvalue_exact_fraction == 1 & c$significance_disagreements == 0
    comparisons[[length(comparisons) + 1L]] <- c
    if (mode == "new_filter" && s$case %in% c("cosmx", "cosmx_stratified")) {
      official_root <- if (s$case == "cosmx") "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead" else
        "/data/dzf/SecAct_2026/results/spatialcellchat_v3_stratified_patched_ladder_20260909/cells010000"
      p <- file.path(official_root, sprintf("cells%06d_full_spatialcellchat_v3_perm00020", s$cells))
      if (file.exists(file.path(p, "records.rds"))) {
        b <- readRDS(file.path(p, "records.rds")); comparison <- compare(b,new)
        os <- read.delim(file.path(p,"summary.tsv"))
        parameters_ok <- identical(as.integer(c(os$cells,os$lr,os$nperm)), as.integer(c(s$cells,s$lr,s$nboot))) &&
          os$radius_um == 30 && os$tol_um == 5 && os$contact_um == 10 && os$scale_distance == 20
        official[[length(official)+1L]] <- cbind(data.frame(case=s$case,cells=s$cells,
          comparator=if(s$case=="cosmx") "official_original_V3" else "official_V3_dimension_patched",
          source_records=file.path(p,"records.rds"), source_md5=unname(tools::md5sum(file.path(p,"records.rds"))),
          historical_output=TRUE, parameters_checked=parameters_ok), comparison)
      }
    }
  }
}
r <- do.call(rbind,resources); c <- do.call(rbind,comparisons); o <- do.call(rbind,official)
write.table(r,file.path(out,"resources.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(c,file.path(out,"new_vs_old.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(o,file.path(out,"new_vs_archived_official.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
stopifnot(nrow(r)==33L,nrow(c)==22L,all(r$exit_code==0),all(c$pass),
          all(o$parameters_checked),all(o$active_jaccard==1),
          all(o$max_abs_probability_difference<=1e-15),all(o$pvalue_exact_fraction==1),
          all(o$significance_disagreements==0))
cat("PASS:",nrow(r),"independent processes;",nrow(c),"new/old comparisons;",nrow(o),"archived official comparisons\n")
print(r[,c("case","cells","mode","inference_and_filter_seconds","peak_rss_gib")])
