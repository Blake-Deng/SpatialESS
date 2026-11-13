#!/usr/bin/env Rscript
# Local --plan/--build-only modes never publish. GitHub publishing is a later job step.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[[1L]] %in% c("main", "sparkle")) {
  stop("Usage: Rscript .github/scripts/release.R main|sparkle [--validate-only|--plan|--build-only]")
}
mode <- args[[1L]]
flags <- args[-1L]
if (length(flags) > 1L || any(!flags %in% c("--validate-only", "--plan", "--build-only"))) {
  stop("Supply at most one supported local mode.")
}
workspace <- normalizePath(Sys.getenv("GITHUB_WORKSPACE", unset = "."), mustWork = TRUE)
tag <- Sys.getenv("RELEASE_TAG")
read_description <- function(path, expected_package) {
  d <- read.dcf(file.path(path, "DESCRIPTION"))
  if (!identical(unname(d[1L, "Package"]), expected_package)) stop("Unexpected package at ", path)
  version <- unname(d[1L, "Version"])
  if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+(\\.[0-9]+)?$", version)) {
    stop("Release version must have three numeric fields, or four for a development version.")
  }
  list(package = expected_package, version = version, source = path)
}
if (mode == "main") {
  packages <- list(read_description(workspace, "SpatialESS"))
  prefix <- "v"
} else {
  packages <- list(read_description(file.path(workspace, "vendor", "SpatialESS"), "SpatialESS"),
                   read_description(workspace, "SpatialESSSPARKLE"))
  prefix <- "sparkle-v"
}
version <- packages[[1L]]$version
if (!all(vapply(packages, function(p) identical(p$version, version), logical(1)))) {
  stop("The bridge and bundled engine versions must be identical.")
}
expected_tag <- paste0(prefix, version)
if (!identical(tag, expected_tag)) stop("Tag must equal ", expected_tag, "; received ", sQuote(tag))
prerelease <- length(strsplit(version, ".", fixed = TRUE)[[1L]]) == 4L
cat("Validated tag:", tag, "\nVersion:", version, "\nPrerelease:", prerelease, "\n")
output <- Sys.getenv("GITHUB_OUTPUT")
if (nzchar(output)) {
  cat("version=", version, "\nprerelease=", tolower(as.character(prerelease)), "\n",
      file = output, append = TRUE, sep = "")
}
if ("--validate-only" %in% flags) quit(status = 0L)
asset_dir <- Sys.getenv("RELEASE_ASSET_DIR", unset = file.path(tempdir(), "spatialess-release-assets"))
if ("--plan" %in% flags) {
  for (p in packages) {
    cat("BUILD:", p$source, "->", file.path(asset_dir, paste0(p$package, "_", p$version, ".tar.gz")), "\n")
    cat("CHECK:", p$package, "--no-manual --no-build-vignettes; error_on=warning\n")
  }
  if (mode == "sparkle") cat("INSTALL: bundled engine into a temporary private R library before checking the bridge.\n")
  cat("No files built, packages installed, checks run, or releases published.\n")
  quit(status = 0L)
}
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)
asset_dir <- normalizePath(asset_dir, mustWork = TRUE)
r <- file.path(R.home("bin"), "R")
build <- function(p) {
  old <- setwd(asset_dir)
  on.exit(setwd(old), add = TRUE)
  status <- system2(r, c("CMD", "build", "--no-build-vignettes", shQuote(p$source)))
  if (status != 0L) stop("Source build failed for ", p$package)
  archive <- file.path(asset_dir, paste0(p$package, "_", p$version, ".tar.gz"))
  if (!file.exists(archive)) stop("Expected archive not produced: ", archive)
  archive
}
archives <- vapply(packages, build, character(1))
if ("--build-only" %in% flags) {
  cat("Build-only local dry run; no checks, installation or publishing performed.\n")
  print(archives)
  quit(status = 0L)
}
if (!requireNamespace("rcmdcheck", quietly = TRUE)) stop("rcmdcheck is required for release checks.")
private_lib <- NULL
if (mode == "sparkle") {
  private_lib <- file.path(Sys.getenv("RUNNER_TEMP", unset = tempdir()), "spatialess-release-Rlib")
  dir.create(private_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(private_lib, .libPaths()))
  status <- system2(r, c("CMD", "INSTALL", shQuote(paste0("--library=", private_lib)), shQuote(archives[[1L]])))
  if (status != 0L) stop("Installing the bundled engine into the private library failed.")
}
for (i in seq_along(packages)) {
  rcmdcheck::rcmdcheck(
    path = archives[[i]], args = c("--no-manual", "--no-build-vignettes"),
    check_dir = file.path(Sys.getenv("RUNNER_TEMP", unset = tempdir()), paste0(packages[[i]]$package, "-release-check")),
    libpath = .libPaths(), error_on = "warning"
  )
}
cat("All release package checks passed. Artifacts ready:\n")
print(archives)
