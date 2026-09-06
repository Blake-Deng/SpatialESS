#!/usr/bin/env Rscript

root <- Sys.getenv(
  "SPATIALESS_PUBLICATION_OUT",
  "/data/dzf/SecAct_2026/results/spatialess_publication_validation")
parse_elapsed <- function(value) {
  fields <- as.numeric(strsplit(value, ":", fixed = TRUE)[[1L]])
  if (length(fields) == 3L) fields[1L] * 3600 + fields[2L] * 60 + fields[3L]
  else if (length(fields) == 2L) fields[1L] * 60 + fields[2L] else fields[1L]
}
time_field <- function(lines, label) {
  hit <- grep(label, lines, fixed = TRUE, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub("^.*:[[:space:]]*", "", tail(hit, 1L)))
}
read_run <- function(directory) {
  summary <- read.delim(file.path(directory, "summary.tsv"),
                        check.names = FALSE)
  lines <- readLines(file.path(directory, "time_verbose.txt"), warn = FALSE)
  elapsed <- grep("Elapsed (wall clock)", lines, fixed = TRUE, value = TRUE)
  elapsed <- trimws(sub("^.*\\):[[:space:]]*", "", tail(elapsed, 1L)))
  summary$process_wall_seconds <- parse_elapsed(elapsed)
  summary$peak_rss_kb <- as.numeric(
    time_field(lines, "Maximum resident set size (kbytes)"))
  summary$peak_rss_gib <- summary$peak_rss_kb / 1024^2
  summary$run_directory <- directory
  summary
}
directories <- list.dirs(root, recursive = FALSE, full.names = TRUE)
directories <- directories[file.exists(file.path(directories, "summary.tsv")) &
                           file.exists(file.path(directories, "time_verbose.txt"))]
runs <- do.call(rbind, lapply(directories, read_run))
runs <- runs[order(runs$experiment, runs$nperm, runs$diffusion_radius_um,
                   runs$replicate), ]
write.table(runs, file.path(root, "publication_validation_runs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

record_key <- function(value)
  paste(value$mechanism, value$sender_group_name, value$receiver_group_name,
        value$interaction_name, sep = "|")
compare_records <- function(reference, candidate) {
  ref_key <- record_key(reference); candidate_key <- record_key(candidate)
  common <- intersect(ref_key, candidate_key)
  ref_index <- match(common, ref_key); candidate_index <- match(common, candidate_key)
  ref_sig <- ref_key[reference$pvalue <= 0.05]
  candidate_sig <- candidate_key[candidate$pvalue <= 0.05]
  union_sig <- union(ref_sig, candidate_sig)
  data.frame(
    reference_records = nrow(reference), candidate_records = nrow(candidate),
    active_jaccard = length(common) / length(union(ref_key, candidate_key)),
    max_abs_probability_difference = max(abs(
      reference$probability[ref_index] - candidate$probability[candidate_index])),
    max_abs_pvalue_difference = max(abs(
      reference$pvalue[ref_index] - candidate$pvalue[candidate_index])),
    exact_probability_fraction = mean(
      reference$probability[ref_index] == candidate$probability[candidate_index]),
    exact_pvalue_fraction = mean(
      reference$pvalue[ref_index] == candidate$pvalue[candidate_index]),
    probability_spearman = suppressWarnings(cor(
      reference$probability[ref_index], candidate$probability[candidate_index],
      method = "spearman")),
    pvalue_spearman = suppressWarnings(cor(
      reference$pvalue[ref_index], candidate$pvalue[candidate_index],
      method = "spearman")),
    significance_jaccard = if (length(union_sig))
      length(intersect(ref_sig, candidate_sig)) / length(union_sig) else 1,
    significance_threshold_disagreements =
      length(setdiff(union_sig, intersect(ref_sig, candidate_sig))),
    stringsAsFactors = FALSE)
}
read_records <- function(directory)
  readRDS(file.path(directory, "records.rds"))

repeat_runs <- runs[runs$experiment == "reproducibility", ]
repeat_reference <- read_records(repeat_runs$run_directory[1L])
reproducibility <- do.call(rbind, lapply(seq_len(nrow(repeat_runs)), function(i) {
  cbind(replicate = repeat_runs$replicate[i],
        compare_records(repeat_reference,
                        read_records(repeat_runs$run_directory[i])))
}))
write.table(reproducibility, file.path(root, "reproducibility_accuracy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

repeat_resources <- data.frame(
  runs = nrow(repeat_runs),
  median_method_seconds = median(repeat_runs$method_end_to_end_seconds),
  cv_method_seconds = sd(repeat_runs$method_end_to_end_seconds) /
    mean(repeat_runs$method_end_to_end_seconds),
  median_process_wall_seconds = median(repeat_runs$process_wall_seconds),
  cv_process_wall_seconds = sd(repeat_runs$process_wall_seconds) /
    mean(repeat_runs$process_wall_seconds),
  median_peak_rss_gib = median(repeat_runs$peak_rss_gib),
  cv_peak_rss = sd(repeat_runs$peak_rss_gib) / mean(repeat_runs$peak_rss_gib))
write.table(repeat_resources,
            file.path(root, "reproducibility_resources.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

perm_runs <- runs[runs$experiment == "permutation", ]
perm_reference <- read_records(perm_runs$run_directory[which.max(perm_runs$nperm)])
permutation <- do.call(rbind, lapply(seq_len(nrow(perm_runs)), function(i) {
  cbind(perm_runs[i, c("nperm", "method_end_to_end_seconds",
                       "process_wall_seconds", "peak_rss_gib",
                       "significant_records_0_05")],
        compare_records(perm_reference, read_records(perm_runs$run_directory[i])))
}))
write.table(permutation, file.path(root, "permutation_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

seed_runs <- runs[runs$experiment == "seed_stability", ]
seed_reference <- read_records(seed_runs$run_directory[1L])
seed_stability <- do.call(rbind, lapply(seq_len(nrow(seed_runs)), function(i) {
  cbind(seed_runs[i, c("seed", "method_end_to_end_seconds",
                       "process_wall_seconds", "peak_rss_gib",
                       "significant_records_0_05")],
        compare_records(seed_reference,
                        read_records(seed_runs$run_directory[i])))
}))
write.table(seed_stability, file.path(root, "seed_stability.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

radius_runs <- runs[runs$experiment == "radius", ]
radius_reference <- read_records(radius_runs$run_directory[
  which(radius_runs$diffusion_radius_um == 35)[1L]])
radius <- do.call(rbind, lapply(seq_len(nrow(radius_runs)), function(i) {
  cbind(radius_runs[i, c("contact_radius_um", "diffusion_radius_um",
                         "contact_edges", "diffusion_edges",
                         "method_end_to_end_seconds", "process_wall_seconds",
                         "peak_rss_gib", "significant_records_0_05")],
        compare_records(radius_reference, read_records(radius_runs$run_directory[i])))
}))
write.table(radius, file.path(root, "radius_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Publication validation complete:", root, "\n")
print(repeat_resources)
print(permutation)
print(seed_stability)
print(radius)
