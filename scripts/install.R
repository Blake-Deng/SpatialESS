#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
external_engine <- Sys.getenv("SPATIALESS_SOURCE", unset = "")
bundled_engine <- file.path(root, "vendor", "SpatialESS")
sibling_engine <- file.path(dirname(root), "SpatialESS")
engine <- if (nzchar(external_engine)) {
  normalizePath(external_engine, mustWork = FALSE)
} else if (dir.exists(bundled_engine)) {
  normalizePath(bundled_engine, mustWork = TRUE)
} else {
  sibling_engine
}
if (!dir.exists(engine)) {
  stop("Could not find the bundled SpatialESS engine. Set SPATIALESS_SOURCE to an external source directory.", call. = FALSE)
}
status_engine <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", shQuote(engine)))
if (status_engine != 0L) stop("SpatialESS installation failed.", call. = FALSE)
status_bridge <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", shQuote(root)))
if (status_bridge != 0L) stop("SpatialESSSPARKLE installation failed.", call. = FALSE)
