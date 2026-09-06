#!/usr/bin/env Rscript

source_root <- Sys.getenv(
  "SPATIALESS_RNA_BENCHMARK_ROOT",
  "/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902"
)
out_dir <- file.path(
  normalizePath("."),
  "benchmark", "ovarian_20260904", "source_data"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_records <- function(engine, method) {
  readRDS(file.path(source_root, "results", engine, "records.rds"))
}
key <- function(x) paste(
  x$interaction_name, x$sender_group_name, x$receiver_group_name, sep = "\r"
)
align_pair <- function(reference, query, reference_name, query_name) {
  reference$key <- key(reference)
  query$key <- key(query)
  all_keys <- union(reference$key, query$key)
  reference <- reference[match(all_keys, reference$key), , drop = FALSE]
  query <- query[match(all_keys, query$key), , drop = FALSE]
  data.frame(
    interaction_name = ifelse(
      !is.na(reference$interaction_name),
      reference$interaction_name, query$interaction_name
    ),
    sender_group_name = ifelse(
      !is.na(reference$sender_group_name),
      reference$sender_group_name, query$sender_group_name
    ),
    receiver_group_name = ifelse(
      !is.na(reference$receiver_group_name),
      reference$receiver_group_name, query$receiver_group_name
    ),
    setNames(list(replace(reference$probability, is.na(reference$probability), 0)),
             paste0(reference_name, "_probability")),
    setNames(list(replace(query$probability, is.na(query$probability), 0)),
             paste0(query_name, "_probability")),
    setNames(list(replace(reference$pvalue, is.na(reference$pvalue), 1)),
             paste0(reference_name, "_pvalue")),
    setNames(list(replace(query$pvalue, is.na(query$pvalue), 1)),
             paste0(query_name, "_pvalue")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
write_gzip_tsv <- function(value, path) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection))
  write.table(
    value, connection, sep = "\t", quote = FALSE, row.names = FALSE
  )
}

official <- read_records(
  "official_rctd_sparkle_full_20260904", "unused"
)
spatialess <- read_records(
  "spatialess_rctd_sparkle_full_20260904", "unused"
)
official_vs_ess <- align_pair(official, spatialess, "official", "spatialess")
write_gzip_tsv(
  official_vs_ess,
  file.path(out_dir, "official_vs_spatialess_records.tsv.gz")
)

raw <- read_records("spatialess_rctd_raw", "unused")
sparkle <- read_records("spatialess_rctd_sparkle", "unused")
raw_vs_sparkle <- align_pair(raw, sparkle, "raw", "sparkle")
write_gzip_tsv(
  raw_vs_sparkle,
  file.path(out_dir, "raw_vs_sparkle_records.tsv.gz")
)

copy_map <- c(
  comparison_metrics.tsv = "results/comparison_rctd_labels/comparison_metrics.tsv",
  lr_aggregate_comparison.tsv = "results/comparison_rctd_labels/lr_aggregate_comparison.tsv",
  COL1A2_SDC4_sender_class_summary.tsv = "results/comparison_rctd_labels/COL1A2_SDC4_sender_class_summary.tsv",
  marker_preservation_summary.tsv = "results/rctd_sparkle_impact/marker_preservation_summary.tsv",
  marker_preservation_by_gene.tsv = "results/rctd_sparkle_impact/marker_preservation_by_gene.tsv",
  expression_reference_summary.csv = paste0(
    "results/expression_reference_rctd_labels/method_level_new/",
    "ovarian_x1000-1800_y300-1100_summary.csv"
  ),
  official_vs_spatialess.tsv = paste0(
    "results/official_rctd_validation/",
    "official_vs_spatialess_sparkle_full_20260904.tsv"
  )
)
for (destination in names(copy_map)) {
  file.copy(
    file.path(source_root, unname(copy_map[[destination]])),
    file.path(out_dir, destination),
    overwrite = TRUE
  )
}

read_peak_kib <- function(path) {
  value <- readLines(path)
  line <- grep("Maximum resident set size", value, value = TRUE)
  as.numeric(sub(".*: ", "", line))
}
read_wall_seconds <- function(path) {
  value <- readLines(path)
  line <- grep("Elapsed (wall clock)", value, value = TRUE, fixed = TRUE)
  stamp <- sub(".*: ", "", line)
  fields <- as.numeric(strsplit(stamp, ":", fixed = TRUE)[[1L]])
  if (length(fields) == 3L) {
    fields[[1L]] * 3600 + fields[[2L]] * 60 + fields[[3L]]
  } else {
    fields[[1L]] * 60 + fields[[2L]]
  }
}

official_summary <- read.delim(file.path(
  source_root, "results", "official_rctd_sparkle_full_20260904", "summary.tsv"
))
ess_summary <- read.delim(file.path(
  source_root, "results", "spatialess_rctd_sparkle_full_20260904", "summary.tsv"
))
official_time <- file.path(
  source_root, "logs", "official_rctd_sparkle_full_20260904.time.txt"
)
ess_time <- file.path(
  source_root, "logs", "spatialess_rctd_sparkle_full_20260904.time.txt"
)
resource <- data.frame(
  method = c("Official SpatialCellChat V3", "SpatialESS"),
  core_seconds = c(
    official_summary$end_to_end_seconds, ess_summary$inference_seconds
  ),
  process_wall_seconds = c(
    read_wall_seconds(official_time), read_wall_seconds(ess_time)
  ),
  peak_rss_gib = c(
    read_peak_kib(official_time), read_peak_kib(ess_time)
  ) / 1048576,
  active_records = c(
    official_summary$active_records, ess_summary$active_records
  ),
  significant_records = c(
    official_summary$significant_records, ess_summary$significant_records
  ),
  stringsAsFactors = FALSE
)
write.table(
  resource, file.path(out_dir, "resource_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

change <- read.delim(file.path(out_dir, "comparison_metrics.tsv"))
metric <- setNames(change$value, change$metric)
transition <- data.frame(
  category = c("Shared significant", "Lost after SPARKLE", "Gained after SPARKLE"),
  records = c(
    metric[["sparkle_significant_records"]] - metric[["significant_gained"]],
    metric[["significant_lost"]],
    metric[["significant_gained"]]
  )
)
write.table(
  transition, file.path(out_dir, "significant_transition.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(resource)
print(transition)
