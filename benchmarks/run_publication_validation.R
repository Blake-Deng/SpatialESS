#!/usr/bin/env Rscript

project_dir <- normalizePath(file.path(dirname(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)), ".."))
.libPaths(c(file.path(project_dir, ".Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(Matrix); library(SpatialESS) })

fixture_root <- Sys.getenv(
  "SECACT_V3_OUT",
  "/data/dzf/SecAct_2026/results/spatialcellchat_v3_headtohead")
out_root <- Sys.getenv(
  "SPATIALESS_PUBLICATION_OUT",
  "/data/dzf/SecAct_2026/results/spatialess_publication_validation")
experiment <- Sys.getenv("SPATIALESS_EXPERIMENT", "reproducibility")
cells <- as.integer(Sys.getenv("SPATIALESS_CELLS", "5000"))
nperm <- as.integer(Sys.getenv("SPATIALESS_NPERM", "20"))
seed <- as.integer(Sys.getenv("SPATIALESS_SEED", "20260804"))
replicate_id <- as.integer(Sys.getenv("SPATIALESS_REPLICATE", "1"))
contact_radius <- as.numeric(Sys.getenv("SPATIALESS_CONTACT_UM", "15"))
diffusion_radius <- as.numeric(Sys.getenv("SPATIALESS_DIFFUSION_UM", "35"))
block_size <- as.numeric(Sys.getenv("SPATIALESS_BLOCK_UM", "100"))
if (!experiment %in% c("reproducibility", "seed_stability",
                       "permutation", "radius") ||
    anyNA(c(cells, nperm, seed, replicate_id)) || cells < 20L || nperm < 1L ||
    replicate_id < 1L ||
    any(!is.finite(c(contact_radius, diffusion_radius, block_size))) ||
    any(c(contact_radius, diffusion_radius, block_size) <= 0)) {
  stop("Invalid publication-validation configuration.", call. = FALSE)
}

fixture_path <- file.path(fixture_root, "secact_cosmx_v3_common_fixture.rds")
fixture <- readRDS(fixture_path)
if (cells > ncol(fixture$expression)) stop("Requested cells exceed fixture.")
idx <- fixture$selection_order[seq_len(cells)]
expression <- fixture$expression[, idx, drop = FALSE]
coordinates <- fixture$coordinates_um[idx, , drop = FALSE]
group <- droplevels(fixture$group[idx])
names(group) <- colnames(expression)

tag <- sprintf(
  "%s_cells%06d_perm%05d_seed%010d_contact%03g_diffusion%03g_block%03g_rep%02d",
  experiment, cells, nperm, seed, contact_radius, diffusion_radius,
  block_size, replicate_id)
run_dir <- file.path(out_root, tag)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

start <- proc.time()[["elapsed"]]
result <- spatialess(
  expression = expression, coordinates = coordinates, group = group,
  lr = fixture$lr_full, complex = fixture$db$complex,
  cofactor = fixture$db$cofactor,
  sample_id = rep("CancerousLiver", cells),
  compartment = rep("all", cells),
  contact_radius = contact_radius, diffusion_radius = diffusion_radius,
  block_size = block_size, inference = "permutation",
  nperm = nperm, seed = seed, Kh = 0.5, n = 1, alpha = 0.05,
  finite_correction = TRUE, tail = "greater_equal",
  graph_weight = "binary", max_edges = 5e8, max_group_pairs = 1e7,
  max_score_entries = 1e8, max_records = 1e8, verbose = FALSE)
elapsed <- proc.time()[["elapsed"]] - start

records <- result$records
key <- paste(records$mechanism, records$sender_group_name,
             records$receiver_group_name, records$interaction_name, sep = "|")
if (anyDuplicated(key)) stop("Output record key is not unique.")
records <- records[order(key), , drop = FALSE]
records$cells <- cells
records$nperm <- nperm
saveRDS(records, file.path(run_dir, "records.rds"), compress = FALSE)
write.table(records, gzfile(file.path(run_dir, "records.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(result$diagnostics,
            file.path(run_dir, "mechanism_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary <- data.frame(
  experiment = experiment, replicate = replicate_id,
  dataset = fixture$parameters$dataset, cells = cells,
  signaling_genes = nrow(expression), groups = nlevels(group),
  lr = nrow(fixture$lr_full), nperm = nperm, seed = seed,
  contact_radius_um = contact_radius,
  diffusion_radius_um = diffusion_radius, block_size_um = block_size,
  contact_edges = if ("contact" %in% names(result$graphs))
    length(result$graphs$contact$neighbors) else 0,
  diffusion_edges = if ("diffusion" %in% names(result$graphs))
    length(result$graphs$diffusion$neighbors) else 0,
  active_records = nrow(records),
  significant_records_0_05 = sum(records$pvalue <= 0.05),
  unconfirmed_significant_records = sum(
    records$pvalue <= 0.05 & !records$permutation_confirmed),
  method_end_to_end_seconds = elapsed,
  stringsAsFactors = FALSE)
write.table(summary, file.path(run_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
writeLines(c(
  sprintf("fixture=%s", normalizePath(fixture_path)),
  sprintf("fixture_md5=%s", unname(tools::md5sum(fixture_path))),
  "threads=one; OMP/BLAS thread variables fixed to one by launcher",
  "null=complete cell profiles permuted within spatial blocks",
  "pvalue=(greater_or_equal_null_count + 1) / (nperm + 1)",
  "all discoveries require fixed-permutation confirmation"
), file.path(run_dir, "provenance.txt"))
print(summary)
