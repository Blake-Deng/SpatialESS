#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))
out_dir <- file.path(project_dir, "benchmarks", "results", "official_spatial_fixture")

fixture <- readRDS(file.path(out_dir, "fixture.rds"))
official <- readRDS(file.path(out_dir, "official_cellchat_spatial.rds"))
ess <- readRDS(file.path(out_dir, "spatialess.rds"))

distance <- official@images$distance
p_spatial <- 1 / (distance * fixture$scale_distance)
p_spatial[is.na(distance)] <- 0
diag(p_spatial) <- max(p_spatial)
official_support <- which(p_spatial > 0, arr.ind = TRUE)
official_support_keys <- paste(
  rownames(p_spatial)[official_support[, 1]],
  colnames(p_spatial)[official_support[, 2]], sep = "->"
)
ess_support_keys <- paste(
  ess$support$group_pairs$sender_group_name,
  ess$support$group_pairs$receiver_group_name, sep = "->"
)

pairs <- ess$observed$group_pairs
official_lr <- dimnames(official@net$prob)[[3L]]
official_group <- dimnames(official@net$prob)[[1L]]
official_probability <- mapply(function(sender, receiver, interaction) {
  official@net$prob[
    match(sender, official_group),
    match(receiver, official_group),
    match(interaction, official_lr)
  ]
}, pairs$sender_group_name, pairs$receiver_group_name, pairs$interaction_name)
spatial_factor <- mapply(function(sender, receiver) {
  p_spatial[match(sender, rownames(p_spatial)), match(receiver, colnames(p_spatial))]
}, pairs$sender_group_name, pairs$receiver_group_name)
official_molecular <- ifelse(spatial_factor > 0, official_probability / spatial_factor, NA_real_)
comparison <- data.frame(
  interaction_name = pairs$interaction_name,
  sender = pairs$sender_group_name,
  receiver = pairs$receiver_group_name,
  official_probability = official_probability,
  official_spatial_factor = spatial_factor,
  official_molecular_probability = official_molecular,
  spatialess_molecular_probability = pairs$molecular_probability,
  spatialess_mean_weighted_probability = pairs$mean_weighted_probability
)
comparable <- is.finite(comparison$official_molecular_probability) &
  comparison$official_probability > 0
difference <- abs(
  comparison$official_molecular_probability[comparable] -
    comparison$spatialess_molecular_probability[comparable]
)
summary <- data.frame(
  official_supported_group_pairs = length(official_support_keys),
  spatialess_supported_group_pairs = length(ess_support_keys),
  common_supported_group_pairs = length(intersect(official_support_keys, ess_support_keys)),
  support_jaccard = length(intersect(official_support_keys, ess_support_keys)) /
    length(union(official_support_keys, ess_support_keys)),
  molecular_records_compared = sum(comparable),
  max_molecular_absolute_difference = max(difference),
  molecular_pearson = stats::cor(
    comparison$official_molecular_probability[comparable],
    comparison$spatialess_molecular_probability[comparable]
  ),
  weighted_score_spearman = stats::cor(
    comparison$official_probability[comparable],
    comparison$spatialess_mean_weighted_probability[comparable],
    method = "spearman"
  )
)
write.table(comparison, file.path(out_dir, "probability_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_dir, "comparison_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summary)
if (!is.finite(summary$max_molecular_absolute_difference) ||
    summary$max_molecular_absolute_difference > 1e-12) {
  stop("Shared CellChat molecular probabilities failed equivalence.", call. = FALSE)
}

