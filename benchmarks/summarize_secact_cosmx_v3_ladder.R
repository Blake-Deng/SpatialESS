#!/usr/bin/env Rscript

root <- Sys.getenv(
  "SECACT_V3_OUT",
  "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead"
)
nperm <- as.integer(Sys.getenv("SECACT_V3_NPERM", "20"))
sizes_text <- Sys.getenv(
  "SECACT_V3_SIZES",
  "250 1000 5000 10000 25000 50000 100000 200000 443515"
)
sizes <- as.integer(strsplit(trimws(sizes_text), "[[:space:]]+")[[1L]])
sizes <- sizes[!is.na(sizes)]
methods <- c("spatialess", "spatialcellchat_v3")

time_value <- function(lines, label) {
  hit <- grep(label, lines, fixed = TRUE, value = TRUE)
  if (!length(hit)) return(NA_character_)
  value <- hit[[length(hit)]]
  if (identical(label, "Elapsed (wall clock)")) {
    return(trimws(sub("^.*\\):[[:space:]]*", "", value)))
  }
  trimws(sub("^[^:]+:", "", value))
}
elapsed_seconds <- function(value) {
  if (is.na(value) || !nzchar(value)) return(NA_real_)
  fields <- as.numeric(strsplit(value, ":", fixed = TRUE)[[1L]])
  if (anyNA(fields)) return(NA_real_)
  if (length(fields) == 3L) return(fields[1L] * 3600 + fields[2L] * 60 + fields[3L])
  if (length(fields) == 2L) return(fields[1L] * 60 + fields[2L])
  if (length(fields) == 1L) return(fields)
  NA_real_
}
first_or_na <- function(x, name) {
  if (is.null(x) || !name %in% colnames(x) || !nrow(x)) return(NA)
  x[[name]][[1L]]
}

rows <- list()
index <- 1L
for (cells in sizes) {
  for (method in methods) {
    run_dir <- file.path(root, sprintf(
      "cells%06d_full_%s_perm%05d", cells, method, nperm
    ))
    summary_path <- file.path(run_dir, "summary.tsv")
    time_path <- file.path(run_dir, "time_verbose.txt")
    method_summary <- if (file.exists(summary_path)) {
      read.delim(summary_path, check.names = FALSE)
    } else {
      NULL
    }
    time_lines <- if (file.exists(time_path)) {
      readLines(time_path, warn = FALSE)
    } else {
      character()
    }
    exit_status <- suppressWarnings(as.integer(time_value(time_lines, "Exit status")))
    status <- if (!dir.exists(run_dir)) {
      "pending"
    } else if (file.exists(summary_path) && identical(exit_status, 0L)) {
      "completed"
    } else if (length(time_lines) && !is.na(exit_status)) {
      "failed"
    } else {
      "running_or_interrupted"
    }
    elapsed_text <- time_value(time_lines, "Elapsed (wall clock)")
    rss_kb <- suppressWarnings(as.numeric(time_value(
      time_lines, "Maximum resident set size (kbytes)"
    )))
    rows[[index]] <- data.frame(
      cells = cells,
      panel = "full",
      method = method,
      nperm = nperm,
      status = status,
      exit_status = exit_status,
      process_wall_seconds = elapsed_seconds(elapsed_text),
      peak_rss_kb = rss_kb,
      peak_rss_gib = rss_kb / 1024^2,
      theoretical_one_dense_distance_gib = cells^2 * 8 / 1024^3,
      signaling_genes = first_or_na(method_summary, "signaling_genes"),
      groups = first_or_na(method_summary, "groups"),
      lr = first_or_na(method_summary, "lr"),
      preparation_seconds = first_or_na(method_summary, "preparation_seconds"),
      individual_or_stream_seconds = first_or_na(
        method_summary, "individual_or_stream_seconds"
      ),
      group_permutation_seconds = first_or_na(
        method_summary, "group_permutation_seconds"
      ),
      method_end_to_end_seconds = first_or_na(
        method_summary, "end_to_end_seconds"
      ),
      graph_edges = first_or_na(method_summary, "graph_edges"),
      stored_active_links = first_or_na(method_summary, "stored_active_links"),
      active_group_lr_records = first_or_na(
        method_summary, "active_group_lr_records"
      ),
      significant_p_0_05 = first_or_na(method_summary, "significant_p_0_05"),
      run_directory = run_dir,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
}
performance <- do.call(rbind, rows)
write.table(
  performance, file.path(root, "ladder_method_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

accuracy_rows <- list()
index <- 1L
for (cells in sizes) {
  comparison_dir <- file.path(root, sprintf(
    "cells%06d_full_comparison_perm%05d", cells, nperm
  ))
  comparison_path <- file.path(comparison_dir, "comparison_summary.tsv")
  if (!file.exists(comparison_path)) next
  value <- read.delim(comparison_path, check.names = FALSE)
  value$comparison_directory <- comparison_dir
  accuracy_rows[[index]] <- value
  index <- index + 1L
}
accuracy <- if (length(accuracy_rows)) {
  do.call(rbind, accuracy_rows)
} else {
  data.frame(
    cells = integer(), panel = character(), nperm = integer(),
    stringsAsFactors = FALSE
  )
}
write.table(
  accuracy, file.path(root, "ladder_accuracy_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

ess <- performance[performance$method == "spatialess", ]
v3 <- performance[performance$method == "spatialcellchat_v3", ]
names(ess)[!names(ess) %in% c("cells", "panel", "nperm")] <- paste0(
  names(ess)[!names(ess) %in% c("cells", "panel", "nperm")], "_ess"
)
names(v3)[!names(v3) %in% c("cells", "panel", "nperm")] <- paste0(
  names(v3)[!names(v3) %in% c("cells", "panel", "nperm")], "_v3"
)
pair <- merge(ess, v3, by = c("cells", "panel", "nperm"), all = TRUE)
pair$process_speedup_v3_over_ess <- pair$process_wall_seconds_v3 /
  pair$process_wall_seconds_ess
pair$method_speedup_v3_over_ess <- pair$method_end_to_end_seconds_v3 /
  pair$method_end_to_end_seconds_ess
pair$peak_rss_reduction_fraction <- 1 -
  pair$peak_rss_gib_ess / pair$peak_rss_gib_v3
if (nrow(accuracy)) {
  pair <- merge(pair, accuracy, by = c("cells", "panel", "nperm"), all.x = TRUE)
}
write.table(
  pair, file.path(root, "ladder_pair_overview.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
print(performance[, c(
  "cells", "method", "status", "process_wall_seconds", "peak_rss_gib",
  "active_group_lr_records", "significant_p_0_05"
)])
