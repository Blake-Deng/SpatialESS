# Edit paths, then run this file to write config.rds.
# Relative paths in the resulting config are relative to config.rds, not the shell.
# Reference input is created separately; see docs/VISIUM_HD_RCTD.md.

config <- list(
  data_dir = "/server/data/sample/outs",
  reference_file = "/server/data/reference.rds",
  output_dir = "/server/results/sample_visium_hd_rctd_test",
  bin_sizes = c(8, 16, 32),
  doublet_mode = "full",
  max_cores = 8L,
  min_umi = 100,
  min_genes = 10L,
  batch_size = 5000L,
  seed = 1L,
  keep_fit = FALSE,
  save_inputs = TRUE,
  overwrite = FALSE,
  read_args = list(gene_column = "name", in_tissue = TRUE,
                   feature_policy = "filtered",
                   allow_filtered_aggregation = FALSE),
  create_args = list(),
  run_args = list()
)

# Optional explicit settings per resolution. The values below are examples,
# not universally validated QC thresholds or a rule for choosing RCTD mode.
# config$fit_args_by_bin <- list(
#   bin8 = list(doublet_mode = "doublet", min_umi = 100),
#   bin16 = list(doublet_mode = "full", min_umi = 400),
#   bin32 = list(doublet_mode = "full", min_umi = 1600)
# )

# Optional CSV tissue-position replacements for individual source resolutions:
# config$read_args_by_bin <- list(
#   bin8 = list(positions_file = "/server/data/positions_008.csv"),
#   bin16 = list(positions_file = "/server/data/positions_016.csv"),
#   bin32 = list(source_bin_size = 16,
#                positions_file = "/server/data/positions_016.csv")
# )

# Optional downstream communication. lr_file contains a list with lr, complex,
# cofactor tables. This uses dominant labels of retained bins, not a
# proportion-weighted or single-cell reconstruction of expression.
# config$lr_file <- "/server/data/lr_database.rds"
# config$min_weight <- 0.5
# config$normalization <- "log1p_cp10k"
# config$spatialess_args <- list(
#   lr_missing = "filter", nboot = 100L, seed.use = 1L,
#   interaction.range = 250, contact.range = 10
# )

saveRDS(config, "config.rds")
