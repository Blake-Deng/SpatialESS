args <- commandArgs(TRUE)
root <- normalizePath(args[1]); base <- normalizePath(args[2])
.libPaths(c(file.path(root, "lib"), .libPaths()))
suppressPackageStartupMessages({library(SpatialESS); library(Matrix); library(data.table); library(jsonlite)})
stopifnot(as.character(packageVersion("SpatialESS")) == "0.1.4")
set.seed(1)
cfg <- fromJSON(file.path(root, "config/design.json"))
out <- file.path(root, "inputs/cervical")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, "reports/preparation_sessionInfo.txt"))
parent_path <- file.path(base, "inputs/cervical/full.rds")
cat("LOAD archived cervical full input", format(Sys.time()), "\n"); flush.console()
full <- readRDS(parent_path)
meta <- fread(file.path(base, "inputs/cervical/cells.tsv"))
stopifnot(nrow(full$expression) == 18028L, ncol(full$expression) == 715512L,
          nrow(full$lr) == 3081L, full$nboot == 100L, full$seed.use == 1L,
          identical(colnames(full$expression), rownames(full$coordinates)),
          identical(colnames(full$expression), names(full$group)),
          identical(colnames(full$expression), meta$cell_id), !anyNA(full$group),
          identical(unname(full$coordinates), unname(as.matrix(meta[,.(x_centroid,y_centroid)]))))
for (name in names(cfg$parameters)) stopifnot(isTRUE(all.equal(full[[name]], cfg$parameters[[name]])))
validObject(full$expression)
center <- colMeans(full$coordinates)
selection <- order(rowSums(sweep(full$coordinates, 2, center)^2), seq_len(nrow(meta)))
old_barcodes <- fread(file.path(base, "inputs/cervical/central_100000_barcodes.tsv"))$cell_id
stopifnot(identical(meta$cell_id[selection[seq_len(100000L)]], old_barcodes))
parent_md5 <- unname(tools::md5sum(parent_path))
manifest <- list()
constants <- setdiff(names(full), c("expression", "coordinates", "group"))
for (size in cfg$paired_sizes) {
  id <- paste0("central_", size)
  idx <- selection[seq_len(size)]
  path <- file.path(out, paste0(id, ".rds"))
  a <- full
  a$expression <- full$expression[,idx,drop=FALSE]
  a$coordinates <- full$coordinates[idx,,drop=FALSE]
  a$group <- droplevels(full$group[idx]); names(a$group) <- colnames(a$expression)
  stopifnot(ncol(a$expression) == size, !anyDuplicated(colnames(a$expression)),
            identical(names(a$group), rownames(a$coordinates)),
            identical(colnames(a$expression), meta$cell_id[idx]),
            all(vapply(constants, function(k) identical(a[[k]], full[[k]]), logical(1))))
  validObject(a$expression)
  if (file.exists(path)) {
    saved <- readRDS(path)
    stopifnot(identical(saved, a))
    rm(saved)
  } else {
    temporary <- paste0(path, ".partial")
    saveRDS(a, temporary, compress=FALSE)
    stopifnot(file.rename(temporary, path))
  }
  selected_meta <- copy(meta[idx]); selected_meta[,subset_order:=seq_len(.N)]
  barcode_path <- file.path(out, paste0(id, "_barcodes.tsv"))
  fwrite(selected_meta, barcode_path, sep="\t")
  fwrite(selected_meta[,.(cells=.N),by=group], file.path(out,paste0(id,"_groups.tsv")), sep="\t")
  manifest[[id]] <- data.table(tissue="cervical", fixture=id, N=size,
    genes=nrow(a$expression), lr=nrow(a$lr), groups=nlevels(a$group),
    expression_nonzeros=length(a$expression@x), nboot=a$nboot, seed_use=a$seed.use,
    sampling_design="central_nearest_mean_coordinate", sampling_randomness="none",
    center_x=center[1], center_y=center[2],
    max_distance_um=max(sqrt(rowSums(sweep(a$coordinates,2,center)^2))),
    effective_radius_um=a$interaction.range+a$tol,
    scale_distance=a$scale.distance, fixture_md5=unname(tools::md5sum(path)),
    barcodes_md5=unname(tools::md5sum(barcode_path)), parent_full_md5=parent_md5,
    nested_with_previous_100k=TRUE, all_parameters_unchanged=TRUE)
  fwrite(rbindlist(manifest), file.path(out,"fixture_manifest.tsv"), sep="\t")
  cat("FIXTURE VALIDATED", id, "genes", nrow(a$expression), "LR", nrow(a$lr),
      "nonzeros", length(a$expression@x), format(Sys.time()), "\n"); flush.console()
  rm(a, selected_meta); gc()
}
stopifnot(identical(fread(file.path(out,"central_400000_barcodes.tsv"))$cell_id[seq_len(200000L)],
                    fread(file.path(out,"central_200000_barcodes.tsv"))$cell_id))
write_json(list(status="COMPLETE", parent_input=parent_path, parent_md5=parent_md5,
  sampling_rule=cfg$sampling_rule, sampling_seed=NULL, permutation_seed=1L,
  sizes=cfg$paired_sizes, full_cell_count=ncol(full$expression),
  unchanged_fields=constants, nested_100k_200k_400k=TRUE,
  finished=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z")),
  file.path(out,"fixtures_complete.json"), pretty=TRUE, auto_unbox=TRUE)
cat("PREPARATION COMPLETE", format(Sys.time()), "\n")
