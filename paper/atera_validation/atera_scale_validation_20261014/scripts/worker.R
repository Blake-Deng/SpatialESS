args <- commandArgs(TRUE)
root <- normalizePath(args[1]); tissue <- args[2]; fixture <- args[3]
engine <- args[4]; phase <- args[5]; nboot <- as.integer(args[6]); tag <- args[7]
.libPaths(c(file.path(root, 'lib'), .libPaths()))
suppressPackageStartupMessages({library(SpatialESS); library(Matrix); library(data.table); library(jsonlite)})
stopifnot(as.character(packageVersion('SpatialESS')) == '0.1.4')
options(future.globals.maxSize=Inf)
if (requireNamespace('future', quietly=TRUE)) future::plan(future::sequential)
out <- file.path(root, 'results', tag); dir.create(out, recursive=TRUE, showWarnings=FALSE)
cfg <- fromJSON(file.path(root, 'config/design.json'))
load_start <- proc.time()[['elapsed']]
a <- readRDS(file.path(root, 'inputs', tissue, paste0(fixture, '.rds')))
if (phase == 'smoke') {
  idx <- seq_len(min(250L, ncol(a$expression)))
  a$expression <- a$expression[,idx,drop=FALSE]
  a$coordinates <- a$coordinates[idx,,drop=FALSE]
  a$group <- droplevels(a$group[idx]); names(a$group) <- colnames(a$expression)
  # Deterministic category-spanning pilot, independent of inferred results.
  picked <- unlist(lapply(split(seq_len(nrow(a$lr)), a$lr$annotation), function(i) {
    i[unique(round(seq(1,length(i),length.out=min(8L,length(i)))))]
  }), use.names=FALSE)
  a$lr <- a$lr[sort(picked),,drop=FALSE]
}
a$nboot <- nboot
stopifnot(identical(colnames(a$expression), rownames(a$coordinates)),
          identical(colnames(a$expression), names(a$group)), !anyNA(a$group),
          nrow(a$expression) == 18028L, nrow(a$lr)>0)
loading_seconds <- proc.time()[['elapsed']] - load_start
cat('START', format(Sys.time()), tissue, fixture, engine, phase,
    'N', ncol(a$expression), 'genes', nrow(a$expression), 'LR', nrow(a$lr),
    'permutations', nboot, '\n'); flush.console()
writeLines(capture.output(sessionInfo()), file.path(out,'sessionInfo.txt'))
saveRDS(list(tissue=tissue, fixture=fixture, engine=engine, phase=phase,
             parameters=a[setdiff(names(a),c('expression','coordinates','group','complex','cofactor'))],
             expression_max=max(a$expression@x), genes=rownames(a$expression),
             barcodes=colnames(a$expression), group_levels=levels(a$group),
             group_counts=table(a$group), rng_kind=RNGkind(), package_version='0.1.4'),
        file.path(out,'contract.rds'))
start <- proc.time()[['elapsed']]
stages <- list()
graph_metrics <- NULL
if (engine == 'SpatialESS') {
  result <- do.call(spatialess, a)
  inference_seconds <- proc.time()[['elapsed']] - start
  records <- as.data.table(result$records)
  graph <- result$graph; E <- length(graph$neighbors); N <- ncol(a$expression)
  checks <- validate_spatial_csr_graph(graph)
  stopifnot(isTRUE(attr(checks,'valid')), length(graph$offsets)==N+1,
            tail(graph$offsets,1)==E, length(graph$distances)==E,
            all(is.finite(graph$distances)), all(graph$distances>0),
            E<=as.double(N)^2)
  array_names <- c('offsets','neighbors','distances','weights','degree','weighted_degree')
  bytes <- vapply(graph[array_names], function(x) length(x)*if (is.integer(x)) 4 else 8, numeric(1))
  object_bytes <- vapply(graph[array_names],function(x) as.double(object.size(x)),numeric(1))
  graph_metrics <- data.table(tissue=tissue,fixture=fixture,phase=phase,N=N,E=E,
    N_squared=as.double(N)^2,mean_out_degree=E/N,density=E/as.double(N)^2,
    offset_bytes=bytes['offsets'],neighbor_bytes=bytes['neighbors'],
    distance_bytes=bytes['distances'],csr_core_array_bytes=sum(bytes[1:3]),
    additional_required_array_bytes=sum(bytes[4:6]),all_array_bytes=sum(bytes),
    all_array_object_bytes=sum(object_bytes),dense_float64_bytes=as.double(N)^2*8,
    effective_radius_um=a$interaction.range+a$tol,scale_distance=a$scale.distance)
  fwrite(graph_metrics,file.path(out,'graph_metrics.tsv'),sep='\t')
  writeLines(capture.output(checks),file.path(out,'graph_validation.txt'))
  if (phase == 'smoke') {
    dense <- as.matrix(dist(a$coordinates*a$ratio)); diag(dense) <- Inf
    expected <- sum(dense<=a$interaction.range+a$tol)
    stopifnot(expected == E)
    cat('GRAPH ORACLE PASS',E,'directed non-self edges\n')
  }
  if (fixture == 'full') {
    matrix_info <- fromJSON(file.path(root,'inputs',tissue,'matrix_complete.json'))
    stopifnot(E == matrix_info$full_graph_edges_ckdtree_diagnostic)
    cat('INDEPENDENT CKDTREE EDGE COUNT PASS',E,'\n')
  }
  fwrite(data.table(interaction_name=result$lr$interaction_name,
                   active_edges=result$diagnostics$active_edges_per_lr),
         file.path(out,'lr_active_edges.tsv'),sep='\t')
  saveRDS(result$diagnostics,file.path(out,'diagnostics.rds'))
  # Save a reusable full graph once, rather than retaining every inference object.
  if (fixture == 'full' && !file.exists(file.path(root,'inputs',tissue,'full_graph.rds'))) {
    saveRDS(graph,file.path(root,'inputs',tissue,'full_graph.rds'),compress=FALSE)
  }
} else if (engine == 'SpatialCellChat') {
  suppressPackageStartupMessages(library(SpatialCellChat))
  db <- readRDS(file.path(root,'inputs/database.rds'))
  prep <- proc.time()[['elapsed']]
  chat <- createSpatialCellChat(object=a$expression,
    meta=data.frame(labels=a$group,row.names=colnames(a$expression)),group.by='labels',
    datatype='spatial',coordinates=a$coordinates,
    spatial.factors=list(ratio=a$ratio,tol=a$tol),do.sparse=TRUE)
  chat@DB <- db
  # Fixed-LR benchmark: both engines normalize over the same full measured panel.
  # This deliberately bypasses overexpression-based LR selection.
  chat@data.signaling <- a$expression
  chat@LR$LRsig <- a$lr
  stages$reference_preparation_seconds <- proc.time()[['elapsed']] - prep
  individual <- proc.time()[['elapsed']]
  chat <- computeCommunProb(chat,LR.use=a$lr,raw.use=TRUE,Kh=a$Kh,n=a$n,
    distance.use=TRUE,interaction.range=a$interaction.range,scale.distance=a$scale.distance,
    use.AGAN=a$use.AGAN,contact.dependent=a$contact.dependent,contact.range=a$contact.range,
    contact.dependent.forced=a$contact.dependent.forced)
  stages$reference_individual_seconds <- proc.time()[['elapsed']] - individual
  averaging <- proc.time()[['elapsed']]
  chat <- computeAvgCommunProb(chat,avg.type='avg',min.percent=a$min.percent,
    min.cells.sr=a$min.cells.sr,do.permutation=nboot>0,nboot=max(1L,nboot),
    seed.use=a$seed.use,colocalization.use=FALSE)
  stages$reference_aggregation_permutation_seconds <- proc.time()[['elapsed']] - averaging
  inference_seconds <- proc.time()[['elapsed']] - start
  prob <- chat@net$prob; pval <- chat@net$pval
  stopifnot(!anyNA(prob),all(is.finite(prob)),all(prob>=0))
  idx <- which(prob>0,arr.ind=TRUE)
  records <- data.table(sender_group_name=dimnames(prob)[[1]][idx[,1]],
    receiver_group_name=dimnames(prob)[[2]][idx[,2]],
    interaction_name=dimnames(prob)[[3]][idx[,3]],probability=prob[idx],pvalue=pval[idx],
    active_edges=NA_real_)
  lr_idx <- match(records$interaction_name,a$lr$interaction_name)
  records[,ligand:=a$lr$ligand[lr_idx]]; records[,receptor:=a$lr$receptor[lr_idx]]
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  write_json(list(package_version=as.character(packageVersion('SpatialCellChat')),
                  full_panel_signaling_override=TRUE,workflow='fixed_LR_computational_benchmark',
                  stages=stages),file.path(out,'reference_provenance.json'),pretty=TRUE,auto_unbox=TRUE)
} else stop('Unknown engine')
required <- c('sender_group_name','receiver_group_name','interaction_name','probability','pvalue')
key_cols <- required[1:3]
if (nrow(records)==0L && !all(required %in% names(records))) {
  records <- data.table(sender_group_name=character(),receiver_group_name=character(),
    interaction_name=character(),probability=double(),pvalue=double(),active_edges=double())
}
stopifnot(all(required %in% names(records)), !anyNA(records$probability),
          all(is.finite(records$probability)),all(records$probability>=0),
          !anyNA(records$pvalue),all(records$pvalue>=0),all(records$pvalue<=1),
          !anyDuplicated(records[,..key_cols]))
setorderv(records,c('sender_group_name','receiver_group_name','interaction_name'))
saveRDS(records,file.path(out,'records.rds'))
fwrite(records,file.path(out,'records.tsv'),sep='\t')
summary <- data.table(tag=tag,tissue=tissue,fixture=fixture,engine=engine,phase=phase,
  N=ncol(a$expression),genes=nrow(a$expression),lr=nrow(a$lr),groups=nlevels(a$group),
  nboot=nboot,seed_use=a$seed.use,loading_seconds=loading_seconds,
  inference_seconds=inference_seconds,active_records=sum(records$probability>0),
  retained_records=nrow(records),raw_p_below_005=if(nboot>0)sum(records$pvalue<.05)else NA_integer_,
  graph_edges=if(is.null(graph_metrics))NA_real_ else graph_metrics$E,
  input_md5=unname(tools::md5sum(file.path(root,'inputs',tissue,paste0(fixture,'.rds')))),
  comparison_scope='fixed LR; full measured panel; no overexpression selection')
fwrite(summary,file.path(out,'summary.tsv'),sep='\t')
write_json(list(status='COMPLETE',tag=tag,inference_seconds=inference_seconds,
                finished=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')),
           file.path(out,'complete.json'),pretty=TRUE,auto_unbox=TRUE)
cat('COMPLETE',tag,'inference_seconds',inference_seconds,'records',nrow(records),'\n')
