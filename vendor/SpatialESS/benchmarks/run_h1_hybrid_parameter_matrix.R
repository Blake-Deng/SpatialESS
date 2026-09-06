#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
project_dir <- dirname(dirname(script_path))

env_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(
    Sys.getenv(name, unset = as.character(default))
  ))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(sprintf("%s must be one positive integer.", name), call. = FALSE)
  }
  value
}

env_numeric_vector <- function(name, default, lower = -Inf, upper = Inf) {
  text <- Sys.getenv(name, unset = "")
  value <- if (nzchar(text)) {
    suppressWarnings(as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]]))
  } else {
    as.numeric(default)
  }
  if (!length(value) || any(!is.finite(value)) ||
      any(value < lower | value > upper)) {
    stop(sprintf("%s contains invalid numeric levels.", name), call. = FALSE)
  }
  unique(value)
}

n_replicates <- env_integer("SPATIALESS_MATRIX_REPLICATES", 25L)
nperm <- env_integer("SPATIALESS_MATRIX_NPERM", 499L)
jobs <- env_integer("SPATIALESS_MATRIX_JOBS", 8L)
base_seed <- env_integer("SPATIALESS_MATRIX_SEED", 20260801L)
grid_levels <- as.integer(env_numeric_vector(
  "SPATIALESS_MATRIX_GRID_LEVELS", c(8, 12, 20), lower = 2
))
zero_levels <- env_numeric_vector(
  "SPATIALESS_MATRIX_ZERO_LEVELS", c(0, 0.25, 0.5, 0.75),
  lower = 0, upper = 0.95
)
tie_levels <- env_numeric_vector(
  "SPATIALESS_MATRIX_TIE_LEVELS", c(0, 0.25, 1), lower = 0
)
dependence_levels <- env_numeric_vector(
  "SPATIALESS_MATRIX_DEPENDENCE_LEVELS", c(0, 0.5, 0.85),
  lower = 0, upper = 0.99
)
screening_margin <- env_numeric_vector(
  "SPATIALESS_MATRIX_SCREENING_MARGIN", 0.01,
  lower = 0, upper = 0.949999
)[[1L]]
max_zero_fraction <- env_numeric_vector(
  "SPATIALESS_MATRIX_MAX_ZERO_FRACTION", 0.10,
  lower = 0, upper = 0.95
)[[1L]]

scenarios <- expand.grid(
  grid_side = grid_levels,
  zero_fraction = zero_levels,
  tie_step = tie_levels,
  lr_dependence = dependence_levels,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
scenarios$scenario_id <- sprintf("scenario_%03d", seq_len(nrow(scenarios)))
scenarios$seed <- base_seed + seq_len(nrow(scenarios)) * 1000000L
if (any(as.double(scenarios$seed) +
        as.double(n_replicates) * nperm * 6 > .Machine$integer.max)) {
  stop("Scenario seeds and permutation ranges exceed the R integer limit.")
}

root <- Sys.getenv(
  "SPATIALESS_MATRIX_OUT",
  file.path(
    project_dir, "benchmarks", "results",
    sprintf("h1_hybrid_parameter_matrix_rep%d_perm%d", n_replicates, nperm)
  )
)
dir.create(root, recursive = TRUE, showWarnings = FALSE)
write.table(
  scenarios, file.path(root, "scenario_design.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

benchmark_script <- file.path(
  project_dir, "benchmarks", "run_h1_multilr_fdr_calibration.R"
)
rscript <- file.path(
  "/home/dzf/miniforge3/envs/cellchat-acceleration/bin", "Rscript"
)
path_value <- paste(
  "/home/dzf/miniforge3/envs/cellchat-acceleration/bin",
  "/usr/bin", "/bin", sep = ":"
)

run_scenario <- function(index) {
  scenario <- scenarios[index, , drop = FALSE]
  out_dir <- file.path(root, scenario$scenario_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "run.log")
  time_file <- file.path(out_dir, "time_verbose.txt")
  env_args <- c(
    sprintf("SPATIALESS_MULTILR_REPLICATES=%d", n_replicates),
    sprintf("SPATIALESS_MULTILR_NPERM=%d", nperm),
    sprintf("SPATIALESS_MULTILR_SEED=%d", scenario$seed),
    sprintf("SPATIALESS_MULTILR_GRID_SIDE=%d", scenario$grid_side),
    "SPATIALESS_MULTILR_BLOCK_SIDE=5",
    sprintf(
      "SPATIALESS_MULTILR_ZERO_FRACTION=%.8g",
      scenario$zero_fraction
    ),
    sprintf("SPATIALESS_MULTILR_TIE_STEP=%.8g", scenario$tie_step),
    sprintf(
      "SPATIALESS_MULTILR_DEPENDENCE=%.8g",
      scenario$lr_dependence
    ),
    sprintf(
      "SPATIALESS_MULTILR_SCREENING_MARGIN=%.8g",
      screening_margin
    ),
    sprintf(
      "SPATIALESS_MULTILR_MAX_ZERO_FRACTION=%.8g",
      max_zero_fraction
    ),
    sprintf("SPATIALESS_MULTILR_OUT=%s", out_dir),
    sprintf("PATH=%s", path_value),
    "OPENBLAS_NUM_THREADS=8",
    "OMP_NUM_THREADS=8",
    "MKL_NUM_THREADS=8"
  )
  status <- system2(
    "env",
    c(
      env_args, "/usr/bin/time", "-v", "-o", time_file,
      rscript, benchmark_script
    ),
    stdout = log_file, stderr = log_file
  )
  data.frame(
    scenario_id = scenario$scenario_id,
    status = as.integer(status),
    stringsAsFactors = FALSE
  )
}

message(sprintf(
  "Running %d scenarios with %d concurrent jobs, %d replicates and %d permutations...",
  nrow(scenarios), jobs, n_replicates, nperm
))
started <- Sys.time()
status <- parallel::mclapply(
  seq_len(nrow(scenarios)), run_scenario,
  mc.cores = min(jobs, nrow(scenarios)), mc.preschedule = FALSE
)
status <- do.call(rbind, status)
write.table(
  status, file.path(root, "run_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
if (any(status$status != 0L)) {
  warning(sprintf("%d scenario runs failed.", sum(status$status != 0L)))
}

completed <- status$scenario_id[status$status == 0L]
summaries <- lapply(completed, function(id) {
  path <- file.path(root, id, "summary.tsv")
  if (!file.exists(path)) {
    return(NULL)
  }
  result <- read.delim(path, check.names = FALSE)
  result$scenario_id <- id
  result
})
summaries <- Filter(Negate(is.null), summaries)
if (!length(summaries)) {
  stop("No scenario summary was produced.")
}
summary_matrix <- do.call(rbind, summaries)
summary_matrix <- summary_matrix[
  match(completed, summary_matrix$scenario_id), , drop = FALSE
]
summary_matrix$group_cells <- summary_matrix$cells / summary_matrix$groups
summary_matrix$null_hybrid_not_above_permutation <-
  summary_matrix$null_hybrid_fwer <=
    summary_matrix$null_permutation_fwer + 1e-15
summary_matrix$all_confirmation_gates_pass <-
  summary_matrix$null_all_discoveries_permutation_confirmed &
    summary_matrix$alternative_all_discoveries_permutation_confirmed
summary_matrix$exact_rejection_gate_pass <-
  summary_matrix$all_replicates_exact_rejection_sets

write.table(
  summary_matrix, file.path(root, "matrix_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

gate <- data.frame(
  scenarios_planned = nrow(scenarios),
  scenarios_completed = nrow(summary_matrix),
  scenarios_failed = sum(status$status != 0L),
  all_hybrid_fwer_not_above_permutation = all(
    summary_matrix$null_hybrid_not_above_permutation
  ),
  all_discoveries_permutation_confirmed = all(
    summary_matrix$all_confirmation_gates_pass
  ),
  all_rejection_sets_exact = all(
    summary_matrix$exact_rejection_gate_pass
  ),
  minimum_null_threshold_agreement = min(
    summary_matrix$null_hybrid_threshold_agreement
  ),
  minimum_alternative_threshold_agreement = min(
    summary_matrix$alternative_hybrid_threshold_agreement
  ),
  wall_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  stringsAsFactors = FALSE
)
write.table(
  gate, file.path(root, "gate_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

failures <- summary_matrix[
  !summary_matrix$null_hybrid_not_above_permutation |
    !summary_matrix$all_confirmation_gates_pass |
    !summary_matrix$exact_rejection_gate_pass,
  , drop = FALSE
]
write.table(
  failures, file.path(root, "gate_failures.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(root, "sessionInfo.txt"))
print(gate)
