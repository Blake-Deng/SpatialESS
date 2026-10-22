#!/usr/bin/env Rscript
# GSE250346 pulmonary fibrosis Xenium cohort (45 sections, 35 patients):
# per-section spatialess communication scores -> patient random-intercept LMM
# for Control vs Disease, with a slice-level GLM for comparison.
#
# Usage (conda env "spatialess"):
#   Rscript run_spatialess_lmm.R [--out DIR] [--max-sections N] [--max-cells M]
#                                [--nboot K] [--seed S]
# Defaults run the full cohort. --max-sections/--max-cells downsample for
# smoke tests; completed sections are checkpointed under OUT/sections/.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
stopifnot("Run via Rscript run_spatialess_lmm.R" = length(file_arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(dirname(script_dir), "_shared", "spatialess_pipeline.R"))

args <- pipeline_parse_args()
out <- arg_value(args, "out", file.path(script_dir, "spatialess_out"))
max_sections <- if ("max-sections" %in% names(args)) {
  arg_integer(args, "max-sections", 1L)
} else Inf
max_cells <- if ("max-cells" %in% names(args)) {
  arg_integer(args, "max-cells", 1L)
} else Inf
nboot <- arg_integer(args, "nboot", 100L)
seed <- arg_integer(args, "seed", 1L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)

db <- ensure_cellchatdb(file.path(dirname(script_dir), "_shared", "cellchatdb"))
dataset <- load_seurat_dataset(
  file.path(script_dir,
            "GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds"),
  meta_cols = list(
    section = "sample",
    cell_type = c("final_CT", "final_lineage"),
    patient = "patient",
    condition = "disease_status",
    covariates = c("sample_affect", "tma", "run")
  )
)

if ("sections" %in% names(args)) {
  sections <- strsplit(args$sections, ",", fixed = TRUE)[[1L]]
  available <- unique(dataset$metadata$section)
  if (!all(sections %in% available)) {
    stop("Unknown --sections value(s): ",
         paste(setdiff(sections, available), collapse = ", "), call. = FALSE)
  }
} else {
  sections <- sort(unique(dataset$metadata$section))
  if (is.finite(max_sections)) sections <- head(sections, max_sections)
}
log_event("Sections to process:", length(sections))

results <- run_sections(dataset, db, sections, out,
                        max_cells = max_cells, nboot = nboot, seed = seed)

completed <- intersect(sections, names(results))
if (length(completed) < 2L) {
  stop("Fewer than 2 sections completed; cannot fit cohort models.",
       call. = FALSE)
}
if (!setequal(completed, sections)) {
  log_event("WARNING: ", length(sections) - length(completed),
            "section(s) failed; cohort model uses", length(completed))
}

covariates <- intersect(c("sample_affect", "tma", "run"),
                        colnames(dataset$metadata))
section_metadata <- build_section_metadata(dataset$metadata, completed,
                                           covariates = covariates)
section_metadata$condition <- factor(
  section_metadata$condition,
  levels = sort(unique(section_metadata$condition))
)
condition_levels <- levels(section_metadata$condition)
log_event("Condition levels among completed sections:",
          paste(condition_levels, collapse = ", "))

lmm <- glm <- NULL
if (length(condition_levels) < 2L) {
  log_event("Fewer than 2 condition levels among completed sections;",
            "skipping cohort models (section results are checkpointed).")
} else {
  coefficient <- paste0("condition",
                        condition_levels[length(condition_levels)])
  log_event("Testing coefficient:", coefficient,
            " reference:", condition_levels[1L])
  tables <- run_cohort_models(
    results[completed], section_metadata, out,
    design = ~ condition,
    coefficient = coefficient,
    patient_col = "patient_id",
    model = c("lmm", "glm"),
    covariates = covariates,
    prefix = "gse250346"
  )
  lmm <- tables$lmm
  glm <- tables$glm
}

write_run_summary(out, list(
  dataset = "GSE250346_lung_fibrosis",
  sections_requested = length(sections),
  sections_completed = length(completed),
  cells_total = ncol(dataset$counts),
  max_cells_per_section = max_cells,
  nboot = nboot,
  model_lmm = "log1p(score) ~ condition + (1 | patient)",
  model_glm = "log1p(score) ~ condition",
  tested_coefficient = if (is.null(lmm)) "" else
    paste(unique(lmm$coefficient), collapse = ","),
  patients = length(unique(section_metadata$patient_id)),
  features = if (is.null(lmm)) 0L else nrow(lmm),
  lmm_ok = if (is.null(lmm)) 0L else sum(lmm$status == "ok"),
  lmm_significant_q05 = if (is.null(lmm)) 0L else
    sum(lmm$status == "ok" & !is.na(lmm$q_value) & lmm$q_value < 0.05),
  glm_ok = if (is.null(glm)) 0L else sum(glm$status == "ok"),
  glm_significant_q05 = if (is.null(glm)) 0L else
    sum(glm$status == "ok" & !is.na(glm$q_value) & glm$q_value < 0.05),
  elapsed_seconds = proc.time()[["elapsed"]]
))
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
log_event("COMPLETE")
