#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
library(Matrix)
library(SpatialESS)
f <- readRDS(file.path(root, "validation", "fixtures", "cervical_1000.rds"))
expected <- readRDS(file.path(root, "validation", "fixtures", "cervical_1000_expected_records.rds"))
coverage <- filter_cellchat_lr(f$db$interaction, rownames(f$expression), f$db$complex, f$db$cofactor)
stopifnot(nrow(coverage$lr) == 1076L, nrow(coverage$excluded) == 2158L)
args <- c(list(expression = f$expression, coordinates = f$coordinates, group = f$group,
             complex = f$db$complex, cofactor = f$db$cofactor,
             min.percent = .1, min.cells.sr = 5L), f$parameters)
explicit <- do.call(spatialess, c(args, list(lr = coverage$lr)))
in_call <- do.call(spatialess, c(args, list(lr = f$db$interaction, lr_missing = "filter")))
stopifnot(identical(explicit$records, in_call$records))
key <- function(x) paste(x$interaction_name, x$sender_group_name, x$receiver_group_name, sep = "|")
r <- explicit$records[match(key(expected), key(explicit$records)), ]
stopifnot(nrow(expected) == nrow(explicit$records),
          setequal(key(expected), key(explicit$records)),
          max(abs(expected$probability - r$probability)) <= 1e-15,
          all(expected$pvalue == r$pvalue),
          all((expected$pvalue <= 0.05) == (r$pvalue <= 0.05)))
cat("PASS: real cervical 1k input; 1076 LR; 100 permutations; explicit/in-call records identical; frozen official V3 probability tolerance and exact p-values passed.\n")
