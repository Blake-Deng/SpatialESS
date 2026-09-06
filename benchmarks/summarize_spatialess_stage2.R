#!/usr/bin/env Rscript

root <- Sys.getenv(
  "SECACT_INTEGRATED_OUT",
  "/data/dzf/SecAct_2026/results/spatialess_integrated_stage2"
)
parse_elapsed <- function(value) {
  fields <- as.numeric(strsplit(value, ":", fixed = TRUE)[[1L]])
  if (anyNA(fields)) return(NA_real_)
  if (length(fields) == 3L) {
    fields[[1L]] * 3600 + fields[[2L]] * 60 + fields[[3L]]
  } else if (length(fields) == 2L) {
    fields[[1L]] * 60 + fields[[2L]]
  } else fields[[1L]]
}
time_field <- function(lines, label) {
  hit <- grep(label, lines, fixed = TRUE, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub("^.*:[[:space:]]*", "", tail(hit, 1L)))
}
directories <- list.dirs(root, recursive = FALSE, full.names = TRUE)
directories <- directories[grepl(
  "^cells[0-9]+_full_(permutation|guarded_hybrid)_perm[0-9]+$",
  basename(directories)
)]
one <- function(directory) {
  summary_path <- file.path(directory, "summary.tsv")
  time_path <- file.path(directory, "time_verbose.txt")
  summary <- if (file.exists(summary_path)) {
    read.delim(summary_path, check.names = FALSE)
  } else NULL
  lines <- if (file.exists(time_path)) readLines(time_path, warn = FALSE) else character()
  exit <- suppressWarnings(as.integer(time_field(lines, "Exit status")))
  elapsed <- grep("Elapsed (wall clock)", lines, fixed = TRUE, value = TRUE)
  elapsed <- if (length(elapsed)) {
    trimws(sub("^.*\\):[[:space:]]*", "", tail(elapsed, 1L)))
  } else NA_character_
  rss <- suppressWarnings(as.numeric(
    time_field(lines, "Maximum resident set size (kbytes)")
  ))
  get <- function(name) {
    if (!is.null(summary) && name %in% names(summary)) summary[[name]][1L] else NA
  }
  data.frame(
    cells = get("cells"), inference = get("inference"),
    status = if (!is.null(summary) && identical(exit, 0L)) {
      "completed"
    } else if (!is.na(exit)) "failed" else "running_or_incomplete",
    exit_status = exit,
    process_wall_seconds = if (is.na(elapsed)) NA_real_ else parse_elapsed(elapsed),
    peak_rss_gib = rss / 1024^2,
    lr = get("lr"), contact_edges = get("contact_edges"),
    diffusion_edges = get("diffusion_edges"),
    active_records = get("active_records"),
    significant_records_0_05 = get("significant_records_0_05"),
    unconfirmed_significant_records = get("unconfirmed_significant_records"),
    run_directory = directory,
    stringsAsFactors = FALSE
  )
}
result <- if (length(directories)) {
  do.call(rbind, lapply(directories, one))
} else data.frame()
if (nrow(result)) result <- result[order(result$cells, result$inference), ]
write.table(
  result, file.path(root, "stage2_integrated_results.tsv"),
  sep = "	", quote = FALSE, row.names = FALSE, na = "NA"
)
print(result)
