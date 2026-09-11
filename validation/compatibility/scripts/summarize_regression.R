args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
old_bridge <- "/home/dzf/cellchat_acceleration/releases/SpatialESS_main_SPARKLE_20260910/runs"
compare <- function(a, b) {
  key <- function(x) paste(x$interaction_name, x$sender_group_name, x$receiver_group_name, sep = "|")
  stopifnot(!anyDuplicated(key(a)), !anyDuplicated(key(b)))
  ka <- key(a); kb <- key(b)
  union <- union(ka, kb); common <- intersect(ka, kb)
  aa <- a[match(common, ka), ]; bb <- b[match(common, kb), ]
  ja <- if (length(union)) length(common)/length(union) else 1
  delta <- if (length(common)) max(abs(aa$probability - bb$probability)) else 0
  exact <- if (ja == 1) { if (length(common)) mean(aa$pvalue == bb$pvalue) else 1 } else NA_real_
  sig <- function(x, keys) { q <- x$pvalue[match(union, keys)]; !is.na(q) & q < .05 }
  disagreements <- sum(sig(a, ka) != sig(b, kb))
  data.frame(reference_active = nrow(a), candidate_active = nrow(b), active_jaccard = ja,
    max_abs_probability_difference = delta, pvalue_exact_fraction = exact,
    significance_disagreements = disagreements,
    full_records_identical = identical(a, b),
    pass = ja == 1 && delta <= 1e-15 && exact == 1 && disagreements == 0)
}
rows <- list()
for (n in c(1000, 5000, 10000, 50000)) {
  d <- file.path(root, "validation/legacy", paste0("cosmx_", n))
  a <- file.path(d, "baseline"); b <- file.path(d, "candidate")
  m <- compare(readRDS(file.path(a, "records.rds")), readRDS(file.path(b, "records.rds")))
  m$cells <- n
  m$inputs_identical <- identical(readRDS(file.path(a, "input_contract.rds")), readRDS(file.path(b, "input_contract.rds")))
  m$rng_identical <- identical(readRDS(file.path(a, "rng.rds")), readRDS(file.path(b, "rng.rds")))
  ca <- readRDS(file.path(a, "result_contract.rds")); cb <- readRDS(file.path(b, "result_contract.rds"))
  cb$parameters$percentage_filter <- NULL
  m$contract_except_filter_marker_identical <- identical(ca, cb)
  rows[[as.character(n)]] <- m
}
x <- do.call(rbind, rows)
write.table(x, file.path(root, "validation/legacy/comparisons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(x)
stopifnot(all(x$pass), all(x$full_records_identical), all(x$inputs_identical),
          all(x$rng_identical), all(x$contract_except_filter_marker_identical))
history <- read.delim(file.path(root, "work/SpatialESS/validation/lr_interface/new_vs_archived_official.tsv"))
rows <- list()
for (n in c(1000, 5000, 10000, 50000)) {
  h <- history[history$case == "cosmx" & history$cells == n, ]
  stopifnot(nrow(h) == 1L, unname(tools::md5sum(h$source_records)) == h$source_md5)
  candidate <- readRDS(file.path(root, "validation/legacy", paste0("cosmx_", n), "candidate/records.rds"))
  m <- compare(readRDS(h$source_records), candidate)
  m$cells <- n; m$reference <- h$comparator; m$archived_official <- TRUE
  m$reference_md5_verified <- TRUE
  rows[[as.character(n)]] <- m
}
x <- do.call(rbind, rows)
write.table(x, file.path(root, "validation/legacy/archived_official.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(all(x$pass))
bridge_root <- file.path(root, "validation/bridge/runs")
rows <- list()
for (condition in c("raw", "sparkle")) for (n in c(1000, 5000, 16247)) {
  case <- paste(condition, n, sep = "_")
  a <- file.path(old_bridge, case, "new_bridge")
  b <- file.path(bridge_root, case, "new_bridge")
  c <- file.path(bridge_root, case, "direct_main")
  for (mode in c("old_vs_new", "bridge_vs_direct")) {
    left <- if (mode == "old_vs_new") a else b
    right <- if (mode == "old_vs_new") b else c
    m <- compare(readRDS(file.path(left, "records.rds")), readRDS(file.path(right, "records.rds")))
    m$case <- case; m$comparison <- mode
    m$inputs_identical <- identical(readRDS(file.path(left, "input.rds")), readRDS(file.path(right, "input.rds")))
    m$rng_identical <- identical(readRDS(file.path(left, "rng.rds")), readRDS(file.path(right, "rng.rds")))
    ca <- readRDS(file.path(left, "contract.rds")); cb <- readRDS(file.path(right, "contract.rds"))
    ca$parameters$percentage_filter <- NULL; cb$parameters$percentage_filter <- NULL
    m$contract_except_filter_marker_identical <- identical(ca, cb)
    rows[[paste(case, mode)]] <- m
  }
}
x <- do.call(rbind, rows)
write.table(x, file.path(root, "validation/bridge/comparisons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(x)
stopifnot(all(x$pass), all(x$full_records_identical), all(x$inputs_identical),
          all(x$rng_identical), all(x$contract_except_filter_marker_identical))
