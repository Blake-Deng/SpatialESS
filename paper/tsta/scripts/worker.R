args <- commandArgs(TRUE)
root <- normalizePath(args[1]); sample_id <- args[2]
.libPaths(c(file.path(root, "lib"), .libPaths()))
suppressPackageStartupMessages({
  library(SpatialESS); library(Matrix); library(data.table); library(jsonlite); library(Seurat)
})
stopifnot(as.character(packageVersion("SpatialESS")) == "0.1.4",
          normalizePath(find.package("SpatialESS")) == normalizePath(file.path(root, "lib/SpatialESS")))
setDTthreads(1L)
if (requireNamespace("future", quietly=TRUE)) future::plan(future::sequential)
options(future.globals.maxSize=Inf)
cfg <- fromJSON(file.path(root, "config/design.json"))
manifest <- fread(file.path(root, "tables/input_manifest.tsv"), na.strings=NULL)
sample_index <- match(sample_id, manifest$sample_id)
stopifnot(!is.na(sample_index))
item <- manifest[sample_index]
stopifnot(nrow(item) == 1L, item$coordinate_status == "VALID_SAVED_ARRAY_GRID")
out <- file.path(root, "results", sample_id)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(!file.exists(file.path(out, "complete.json")))
started <- proc.time()[["elapsed"]]
stage <- function(label) {
  cat(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sample_id, label, "\n"); flush.console()
  write_json(list(sample_id=sample_id, stage=label, pid=Sys.getpid(),
                  elapsed_seconds=proc.time()[["elapsed"]]-started,
                  updated=format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
             file.path(out, "progress.json"), auto_unbox=TRUE, pretty=TRUE)
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
stage("READ_COUNTS")
con <- gzfile(item$counts_path, "rt")
header <- readLines(con, n=1L); close(con)
spots <- strsplit(header, "\t", fixed=TRUE)[[1]]
if (spots[1] %in% c("", '\"\"')) spots <- spots[-1]
spots <- gsub('^\"|\"$', "", spots)
stopifnot(length(spots) == item$input_spots, !anyDuplicated(spots),
          all(grepl("^[0-9]+x[0-9]+$", spots)))
tab <- fread(cmd=paste("gzip -dc --", shQuote(item$counts_path)),
             header=FALSE, skip=1L, showProgress=FALSE)
stopifnot(ncol(tab) == length(spots)+1L, nrow(tab)>0L)
genes <- as.character(tab[[1L]])
stopifnot(!anyNA(genes), all(nzchar(genes)))
dense <- as.matrix(tab[, -1L]); rm(tab); gc()
storage.mode(dense) <- "double"
stopifnot(!anyNA(dense), all(is.finite(dense)), all(dense >= 0),
          all(dense == floor(dense)))
raw <- as(Matrix(dense, sparse=TRUE), "dgCMatrix"); rm(dense); gc()
colnames(raw) <- spots
duplicate_count <- sum(duplicated(genes))
if (duplicate_count > 0L) {
  unique_genes <- unique(genes)
  map <- match(genes, unique_genes)
  collapse <- sparseMatrix(i=map, j=seq_along(genes), x=1,
                           dims=c(length(unique_genes), length(genes)))
  raw <- as(collapse %*% raw, "dgCMatrix")
  rownames(raw) <- unique_genes
  fwrite(data.table(input_row=seq_along(genes), gene=genes, collapsed_row=map),
         file.path(out, "gene_row_mapping.tsv.gz"), sep="\t")
} else rownames(raw) <- genes
library_size <- Matrix::colSums(raw)
features <- Matrix::colSums(raw > 0)
keep <- library_size > cfg$minimum_library_size_exclusive
if (sum(keep) < 30L) stop("Fewer than 30 spots pass the frozen library-size rule")
spot_meta <- data.table(spot_id=spots, raw_library_size=library_size,
                       detected_genes=features, retained=keep)
raw <- raw[, keep, drop=FALSE]
spots <- colnames(raw)
array_coords <- do.call(rbind, strsplit(spots, "x", fixed=TRUE))
storage.mode(array_coords) <- "double"
stopifnot(all(rowSums(array_coords) %% 2 == 0))
coordinates <- cbind(x_um=array_coords[,1] * 0.5 * sqrt(3) * 100,
                     y_um=array_coords[,2] * 0.5 * 100)
rownames(coordinates) <- spots
stopifnot(!anyNA(coordinates), !anyDuplicated(as.data.frame(coordinates)))
loading_seconds <- proc.time()[["elapsed"]] - started
stage("EXPRESSION_CLUSTERING")
cluster_start <- proc.time()[["elapsed"]]
set.seed(cfg$clustering_seed)
clustering_names <- make.unique(gsub("[_|]", "-", rownames(raw)))
cluster_counts <- raw
rownames(cluster_counts) <- clustering_names
obj <- CreateSeuratObject(counts=cluster_counts, min.cells=0, min.features=0,
                           project=sample_id)
rm(cluster_counts)
obj <- NormalizeData(obj, normalization.method="LogNormalize",
                     scale.factor=cfg$normalization_scale_factor, verbose=FALSE)
obj <- FindVariableFeatures(obj, selection.method="vst", nfeatures=cfg$variable_features,
                             verbose=FALSE)
hvg <- VariableFeatures(obj)
obj <- ScaleData(obj, features=hvg, verbose=FALSE)
npcs <- min(cfg$pca_components, ncol(obj)-1L, length(hvg)-1L)
if (npcs < 2L) stop("Insufficient variable features for frozen clustering method")
obj <- RunPCA(obj, features=hvg, npcs=npcs, seed.use=cfg$clustering_seed,
               verbose=FALSE)
dims <- seq_len(min(cfg$neighbor_pcs, npcs))
obj <- FindNeighbors(obj, dims=dims, k.param=min(cfg$neighbor_k, ncol(obj)-1L), verbose=FALSE)
obj <- FindClusters(obj, resolution=cfg$clustering_resolution, algorithm=1L,
                     random.seed=cfg$clustering_seed, verbose=FALSE)
stopifnot(identical(colnames(obj), spots))
group <- droplevels(factor(paste0("Cluster_", as.character(Idents(obj)))))
names(group) <- spots
fwrite(data.table(original_gene=rownames(raw), clustering_gene=clustering_names,
                  variable_feature=clustering_names %in% hvg),
       file.path(out, "clustering_gene_mapping.tsv.gz"), sep="\t")
rm(obj); gc()
cluster_seconds <- proc.time()[["elapsed"]]-cluster_start
expression <- raw
expression@x <- log1p(expression@x * rep(
  cfg$normalization_scale_factor / Matrix::colSums(raw), diff(raw@p)))
stopifnot(!anyNA(expression@x), all(is.finite(expression@x)),
          identical(colnames(expression), names(group)),
          identical(colnames(expression), rownames(coordinates)))
spot_meta[, `:=`(array_row=NA_real_, array_col=NA_real_, x_um=NA_real_, y_um=NA_real_,
                 expression_cluster=NA_character_)]
spot_meta[keep, `:=`(array_row=array_coords[,1], array_col=array_coords[,2],
                     x_um=coordinates[,1], y_um=coordinates[,2], expression_cluster=as.character(group))]
fwrite(spot_meta, file.path(out, "spot_metadata.tsv.gz"), sep="\t")
fwrite(data.table(group=levels(group), spots=as.integer(table(group))),
       file.path(out, "group_sizes.tsv"), sep="\t")
stage("LR_COVERAGE")
db <- readRDS(file.path(root, "inputs/database.rds"))
coverage <- filter_cellchat_lr(db$interaction, rownames(expression), db$complex, db$cofactor)
fwrite(coverage$audit, file.path(out, "lr_coverage.tsv.gz"), sep="\t")
fwrite(coverage$missing_genes, file.path(out, "lr_missing_genes.tsv.gz"), sep="\t")
lr <- coverage$lr
stopifnot(nrow(lr)>0L)
fwrite(data.table(database_row=rownames(lr), lr), file.path(out, "lr.tsv.gz"), sep="\t")
fwrite(data.table(gene=rownames(expression)), file.path(out, "measured_genes.tsv.gz"), sep="\t")
parameters <- cfg$parameters
parameters$min.cells.sr <- as.integer(parameters$min.cells.sr)
parameters$nboot <- as.integer(cfg$nboot)
parameters$seed.use <- as.integer(cfg$seed_use)
saveRDS(list(sample=as.list(item), parameters=parameters,
             normalization=list(method="log1p(count/library_size*10000)", scale_factor=cfg$normalization_scale_factor),
             grouping_scope=cfg$grouping_scope, seed=cfg$clustering_seed,
             clustering_resolution=cfg$clustering_resolution, actual_pcs=dims,
             gene_symbols=rownames(expression), group_levels=levels(group),
             source_counts_md5=unname(tools::md5sum(item$counts_path)),
             rng_kind=RNGkind(), duplicate_gene_rows_summed=duplicate_count),
        file.path(out, "input_contract.rds"))
rm(raw); gc()
stage("SPATIALESS_INFERENCE_100_PERMUTATIONS")
inference_start <- proc.time()[["elapsed"]]
result <- do.call(spatialess, c(list(expression=expression, coordinates=coordinates,
    group=group, lr=lr, complex=db$complex, cofactor=db$cofactor), parameters))
inference_seconds <- proc.time()[["elapsed"]]-inference_start
stage("VALIDATE_AND_SAVE")
checks <- validate_spatial_csr_graph(result$graph)
graph <- result$graph; N <- ncol(expression); E <- length(graph$neighbors)
stopifnot(isTRUE(attr(checks, "valid")), length(graph$offsets)==N+1,
          tail(graph$offsets, 1)==E, length(graph$distances)==E,
          all(is.finite(graph$distances)), all(graph$distances>0))
records <- as.data.table(result$records)
if (nrow(records)) {
  keys <- c("sender_group_name", "receiver_group_name", "interaction_name")
  stopifnot(!anyDuplicated(records[, ..keys]), all(is.finite(records$probability)),
            all(records$probability >= 0), all(is.finite(records$pvalue)),
            all(records$pvalue>=0 & records$pvalue<=1))
  setorderv(records, keys)
}
fwrite(records, file.path(out, "records.tsv.gz"), sep="\t")
saveRDS(result, file.path(out, "SpatialESS_result.rds"), compress=FALSE)
writeLines(capture.output(checks), file.path(out, "graph_validation.txt"))
byte_names <- c("offsets", "neighbors", "distances", "weights", "degree", "weighted_degree")
array_bytes <- vapply(graph[byte_names], function(x) length(x)*if(is.integer(x)) 4 else 8, numeric(1))
summary <- data.table(sample_id=sample_id, Study=item$Study, Cancer=item$Cancer,
  Sample_Name=item$Sample_Name, engine="SpatialESS", package_version="0.1.4",
  unit="Visium_spot", input_spots=item$input_spots, N=N, genes=nrow(expression),
  groups=nlevels(group), lr=nrow(result$lr), database_lr=nrow(db$interaction),
  nboot=cfg$nboot, seed_use=cfg$seed_use, clustering_seed=cfg$clustering_seed,
  duplicate_gene_rows_summed=duplicate_count,
  loading_seconds=loading_seconds, clustering_seconds=cluster_seconds,
  inference_seconds=inference_seconds, worker_elapsed_seconds=proc.time()[["elapsed"]]-started,
  active_records=if(nrow(records))sum(records$probability>0)else 0L,
  raw_p_below_005=if(nrow(records))sum(records$pvalue<.05)else 0L,
  E=E, mean_out_degree=E/N, graph_density=E/as.double(N)^2,
  csr_core_array_bytes=sum(array_bytes[1:3]), csr_all_array_bytes=sum(array_bytes),
  interaction_range_um=cfg$parameters$interaction.range,
  contact_dependent=FALSE, coordinate_source="saved_array_ID_SecAct_calWeights_transform",
  group_source="per_sample_expression_graph_clusters_not_cell_types")
fwrite(summary, file.path(out, "summary.tsv"), sep="\t")
write_json(list(status="COMPLETE", sample_id=sample_id, N=N, E=E,
                groups=nlevels(group), lr=nrow(result$lr),
                inference_seconds=inference_seconds,
                finished=format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
           file.path(out, "complete.json"), pretty=TRUE, auto_unbox=TRUE)
stage("COMPLETE")
