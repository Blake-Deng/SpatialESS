args <- commandArgs(TRUE)
if (length(args)!=3L) stop('Usage: Rscript fit_multisample.R cohort_input.rds coefficient output.rds')
if (file.exists(args[3])) stop('Output already exists; use a new filename.')
suppressPackageStartupMessages(library(SpatialESS))
input <- readRDS(args[1])
stopifnot(is.list(input$sample_results), is.data.frame(input$sample_metadata))
prepared <- prepare_multisample_communication(results=input$sample_results,
  sample_metadata=input$sample_metadata, unit='sample')
fit <- fit_multisample_communication_lmm(prepared=prepared, design=~condition,
  coefficient=args[2], patient_col='patient_id', min_units=6, min_nonzero=2, REML=FALSE)
saveRDS(fit, args[3])
writeLines(capture.output(sessionInfo()), paste0(args[3], '.sessionInfo.txt'))
