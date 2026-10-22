args <- commandArgs(TRUE)
root <- normalizePath(args[1]); tissue <- args[2]
.libPaths(c(file.path(root, 'lib'), .libPaths()))
suppressPackageStartupMessages({library(SpatialESS); library(Matrix); library(data.table); library(jsonlite)})
stopifnot(as.character(packageVersion('SpatialESS')) == '0.1.4')
cfg <- fromJSON(file.path(root, 'config/design.json'), simplifyVector=FALSE)
out <- file.path(root, 'inputs', tissue)
if (file.exists(file.path(out, 'fixtures_complete.json'))) quit(status=0)
m <- fromJSON(file.path(out, 'matrix_complete.json'))
meta <- fread(file.path(out, 'cells.tsv')); genes <- fread(file.path(out, 'genes.tsv'))
read_binary <- function(name, what, n, size) {
  con <- file(file.path(out, name), 'rb'); on.exit(close(con))
  ans <- readBin(con, what, n=n, size=size, endian='little')
  if (length(ans) != n) stop('Truncated binary array ', name)
  ans
}
cat('LOAD sparse binary', tissue, m$nnz, 'nonzeros\n')
x <- read_binary('x.bin', 'double', m$nnz, 8)
i <- read_binary('i.bin', 'integer', m$nnz, 4)
p <- read_binary('p.bin', 'integer', m$eligible_cells+1, 4)
expression <- new('dgCMatrix', x=x, i=i, p=p,
                  Dim=as.integer(c(m$genes, m$eligible_cells)),
                  Dimnames=list(genes$gene, meta$cell_id))
validObject(expression); rm(x,i,p); gc()
stopifnot(max(expression@x) == m$normalized_expression_max, nrow(expression)==18028L)
db_env <- new.env(); data('CellChatDB.human', package='SpatialCellChat', envir=db_env)
db <- db_env$CellChatDB.human
if (!file.exists(file.path(root,'inputs/database.rds'))) saveRDS(db,file.path(root,'inputs/database.rds'))
coverage <- filter_cellchat_lr(db$interaction, rownames(expression), db$complex, db$cofactor)
saveRDS(coverage, file.path(out,'lr_coverage.rds'))
fwrite(coverage$audit,file.path(out,'lr_coverage.tsv'),sep='\t')
fwrite(coverage$missing_genes,file.path(out,'lr_missing_genes.tsv'),sep='\t')
lr <- coverage$lr
# Freeze the database ordering used by both engines; do not rank by results.
annotation <- factor(lr$annotation,levels=c('Secreted Signaling','ECM-Receptor',
                                          'Non-protein Signaling','Cell-Cell Contact'))
lr <- lr[order(annotation),,drop=FALSE]
lr$annotation <- as.character(lr$annotation)
fwrite(data.table(database_row=rownames(lr),lr),file.path(out,'lr.tsv'),sep='\t')
coordinates <- as.matrix(meta[,.(x_centroid,y_centroid)])
rownames(coordinates) <- meta$cell_id
center <- colMeans(coordinates)
distance <- rowSums(sweep(coordinates,2,center)^2)
selection <- order(distance, seq_len(nrow(meta)))
params <- cfg$parameters; params$scale.distance <- m$scale_distance
params$nboot <- as.integer(cfg$nboot); params$seed.use <- as.integer(cfg$seed_use)
params$max_edges <- cfg$max_edges; params$max_active_edges_per_lr <- cfg$max_active_edges_per_lr
manifest <- list()
save_fixture <- function(id, idx, design, seed=NA_integer_, anchor=center) {
  path <- file.path(out,paste0(id,'.rds'))
  barcode_file <- file.path(out,paste0(id,'_barcodes.tsv'))
  sub_meta <- copy(meta[idx]); sub_meta[,subset_order:=seq_len(.N)]
  if (!file.exists(path)) {
    a <- c(list(expression=expression[,idx,drop=FALSE],
                coordinates=coordinates[idx,,drop=FALSE],
                group=droplevels(factor(meta$group[idx])),lr=lr,
                complex=db$complex,cofactor=db$cofactor), params)
    names(a$group) <- colnames(a$expression)
    saveRDS(a,path,compress=FALSE); rm(a); gc()
    fwrite(sub_meta,barcode_file,sep='\t')
  }
  group_support <- sub_meta[,.(cells=.N),by=group]
  fwrite(group_support,file.path(out,paste0(id,'_groups.tsv')),sep='\t')
  record <- data.frame(tissue=tissue,fixture=id,cells=length(idx),genes=nrow(expression),
    groups=nrow(group_support),lr=nrow(lr),sampling_design=design,sampling_seed=seed,
    center_x=anchor[1],center_y=anchor[2],max_distance_from_anchor=max(sqrt(rowSums(
      sweep(coordinates[idx,,drop=FALSE],2,anchor)^2))),
    scale_distance=m$scale_distance,path=path,barcodes=barcode_file,
    fixture_md5=unname(tools::md5sum(path)),barcodes_md5=unname(tools::md5sum(barcode_file)))
  manifest[[length(manifest)+1L]] <<- record
  cat('FIXTURE',id,'cells',length(idx),'groups',nrow(group_support),'LR',nrow(lr),'\n')
}
for (size in unlist(cfg$central_sizes)) save_fixture(paste0('central_',size),
                      selection[seq_len(size)],'central_nearest_mean_coordinate')
save_fixture('full',seq_len(nrow(meta)),'full_eligible_tissue')
types <- sort(unique(meta$group)); counts <- table(factor(meta$group,levels=types))
by_type <- split(seq_len(nrow(meta)),factor(meta$group,levels=types))
quotas_all <- list()
for (size in unlist(cfg$sampling_sizes)) {
  raw_quota <- size * as.numeric(counts) / sum(counts)
  quota <- floor(raw_quota); remainder <- size-sum(quota)
  if (remainder>0) {
    extra <- order(-(raw_quota-quota),seq_along(types))[seq_len(remainder)]
    quota[extra] <- quota[extra]+1L
  }
  stopifnot(sum(quota)==size,all(quota<=counts))
  quotas_all[[as.character(size)]] <- data.frame(group=types,raw_quota=raw_quota,quota=quota)
  for (rep in seq_len(cfg$sampling_repetitions)) {
    seed <- cfg$sampling_seed+rep-1L
    set.seed(seed)
    idx <- unlist(lapply(seq_along(types),function(j) {
      pool<-by_type[[j]]; pool[sample.int(length(pool),quota[j])]
    }),use.names=FALSE)
    idx <- idx[sample.int(length(idx))]
    save_fixture(paste0('stratified_',size,'_rep',rep),idx,
                 'proportional_official_group_dispersed',seed)
    # Additional deterministic anchors test tissue-region dependence.
    target <- c(unname(quantile(coordinates[,1],probs=c(.25,.5,.75)[rep])),
                unname(median(coordinates[,2])))
    anchor_idx <- which.min(rowSums(sweep(coordinates,2,target)^2))
    anchor <- coordinates[anchor_idx,]
    idx <- order(rowSums(sweep(coordinates,2,anchor)^2),seq_len(nrow(meta)))[seq_len(size)]
    save_fixture(paste0('block_',size,'_rep',rep),idx,
                 'nearest_observed_quantile_anchor',anchor=anchor)
  }
}
fwrite(rbindlist(manifest),file.path(out,'fixture_manifest.tsv'),sep='\t')
saveRDS(quotas_all,file.path(out,'sampling_quotas.rds'))
for (name in names(quotas_all)) fwrite(quotas_all[[name]],
                     file.path(out,paste0('sampling_quotas_',name,'.tsv')),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'preparation_sessionInfo.txt'))
write_json(list(status='COMPLETE',fixtures=length(manifest),cells=nrow(meta),genes=nrow(expression),
                lr=nrow(lr),database_total_lr=nrow(db$interaction),
                central_rule='Nearest cells to arithmetic mean of all eligible centroids; squared Euclidean distance, original order breaks ties; no randomness.',
                stratified_rule='Largest-remainder proportional quotas from all eligible official groups; sorted group order breaks ties; independent R sample.int within groups and final shuffle.',
                block_rule='Nearest observed centroid to (x quartile, median y) is anchor; nearest N centroids by squared distance; original order breaks ties.',
                expression_scope='All 18028 measured Gene Expression rows; normalization uses all measured expression rows; full rows are shared by both engines.',
                lr_rule=cfg$lr_selection),file.path(out,'fixtures_complete.json'),pretty=TRUE,auto_unbox=TRUE)
cat('FIXTURES COMPLETE',tissue,'\n')
