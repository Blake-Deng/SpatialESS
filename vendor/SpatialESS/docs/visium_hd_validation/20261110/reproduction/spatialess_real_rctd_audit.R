#!/usr/bin/env Rscript
# Same-input numerical audit. This tests wrapper/upstream compatibility, not
# biological accuracy or equivalence to independently produced Space Ranger bins.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1L && args[[1L]] %in% c('--help', '-h')) {
  cat(paste0(
    'Usage: Rscript spatialess_real_rctd_audit.R config.rds [OUTPUT_DIR]\n\n',
    'Config RDS fields:\n',
    '  hd_file, reference_file, output_dir: paths (relative to config directory).\n',
    '  bin_size: expected input resolution, 8/16/32; never resamples this audit input.\n',
    '  modes: c("full","doublet") by default; max_cores must be 1.\n',
    '  seed=1L, min_umi=100, min_genes=10L, batch_size=64L.\n',
    '  create_args=list(), run_args=list(): exactly shared upstream settings.\n',
    '  roi32=list(row_start=?,col_start=?,n_rows=5L,n_cols=5L): a rectangle\n',
    '    of complete 32um macro-bins; origin defaults near tissue median.\n',
    '  max_input_bins=500L; require_raw_roi=TRUE; overwrite=FALSE.\n',
    '  tolerance=1e-8: an explicit audit tolerance, not an official RCTD standard.\n',
    '  library_paths: optional character vector prepended to .libPaths().\n\n',
    'reference_file contains list(counts=gene_by_reference_cell, cell_types=named_labels).\n',
    'hd_file is a SpatialESSVisiumHD object with full all_positions. The strict\n',
    'ROI default requires unfiltered raw barcode counts, in_tissue=FALSE provenance.\n',
    'Absent count columns may be zero-read bins, but the positions grid must be complete.\n',
    'Prepared/filtered inputs can be audited only with explicit require_raw_roi=FALSE;\n',
    'that compares implementations on those inputs and does not establish raw ROI completeness.\n',
    'No inference is executed by --help. No network connection is made by this script.\n'
  ))
  quit(status = 0L)
}
if (!length(args) %in% c(1L, 2L)) stop('Expected CONFIG_RDS and optional OUTPUT_DIR; use --help.', call. = FALSE)
config_file <- normalizePath(args[[1L]], mustWork = TRUE)
config_dir <- dirname(config_file)
config <- readRDS(config_file)
allowed <- c('hd_file','reference_file','output_dir','bin_size','modes','max_cores',
  'seed','min_umi','min_genes','batch_size','create_args','run_args','roi32',
  'max_input_bins','require_raw_roi','overwrite','tolerance','library_paths')
if (!is.list(config) || is.null(names(config)) || anyNA(names(config)) ||
    any(!nzchar(names(config))) || anyDuplicated(names(config)) ||
    length(setdiff(names(config), allowed))) stop('Config must be a uniquely named list containing only documented fields.', call. = FALSE)
resolve <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
    stop('Expected a single nonempty path.', call. = FALSE)
  path <- path.expand(path)
  if (!grepl('^(/|[A-Za-z]:[/\\\\])', path)) path <- file.path(config_dir, path)
  path
}
if (length(args) == 2L) config$output_dir <- args[[2L]]
for (key in c('hd_file','reference_file','output_dir')) config[[key]] <- resolve(config[[key]])
if (!is.null(config$library_paths)) {
  if (!is.character(config$library_paths) || anyNA(config$library_paths)) stop('library_paths must be character paths.')
  .libPaths(c(vapply(config$library_paths, resolve, character(1L)), .libPaths()))
}
suppressPackageStartupMessages({ library(Matrix); library(SpatialESS); library(spacexr) })
required <- c('Reference','SpatialRNA','create.RCTD','run.RCTD')
if (!all(required %in% getNamespaceExports('spacexr'))) stop('The original spacexr API is required; the Bioconductor API is different.')
if (!'run_visium_hd_rctd' %in% getNamespaceExports('SpatialESS')) stop('Install the patched SpatialESS first.')
number <- function(x, label, lower, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower ||
      (integer && (x != floor(x) || x > .Machine$integer.max))) stop(label, ' is invalid.', call. = FALSE)
  if (integer) as.integer(x) else x
}
flag <- function(x, label) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(label, ' must be TRUE/FALSE.', call. = FALSE)
  x
}
named_list <- function(x, label) {
  if (is.null(x)) return(list())
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x))))) stop(label, ' must be a uniquely named list.')
  x
}
value <- function(key, default) if (is.null(config[[key]])) default else config[[key]]
settings <- list(
  modes = value('modes', c('full','doublet')),
  max_cores = number(value('max_cores', 1L), 'max_cores', 1, TRUE),
  seed = number(value('seed', 1L), 'seed', 0, TRUE),
  min_umi = number(value('min_umi', 100), 'min_umi', 0),
  min_genes = number(value('min_genes', 10L), 'min_genes', 0, TRUE),
  batch_size = number(value('batch_size', 64L), 'batch_size', 1, TRUE),
  max_input_bins = number(value('max_input_bins', 500L), 'max_input_bins', 1, TRUE),
  require_raw_roi = flag(value('require_raw_roi', TRUE), 'require_raw_roi'),
  overwrite = flag(value('overwrite', FALSE), 'overwrite'),
  tolerance = number(value('tolerance', 1e-8), 'tolerance', 0),
  create_args = named_list(config$create_args, 'create_args'),
  run_args = named_list(config$run_args, 'run_args')
)
if (settings$max_cores != 1L) stop('Use max_cores=1 for this reproducibility audit; performance tests are separate.')
if (!is.character(settings$modes) || !length(settings$modes) || anyNA(settings$modes) ||
    anyDuplicated(settings$modes) || !all(settings$modes %in% c('full','doublet'))) stop('modes must contain full and/or doublet, without duplicates.')
if (any(names(settings$create_args) %in% c('spatialRNA','reference','max_cores'))) stop('create_args cannot override shared data or cores.')
if (!all(names(settings$run_args) %in% c('MAX_MULTI_TYPES','CONFIDENCE_THRESHOLD','DOUBLET_THRESHOLD','MIN_CHANGE_REG'))) stop('Unsupported run_args field.')
if (!'UMI_min' %in% names(settings$create_args)) settings$create_args$UMI_min <- settings$min_umi
outdir <- config$output_dir
if (dir.exists(outdir) && length(list.files(outdir, all.files = TRUE, no.. = TRUE)) && !settings$overwrite)
  stop('Output directory is not empty; use a fresh path or explicitly set overwrite=TRUE.', call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(outdir)) stop('Cannot create output directory.')
saveRDS(list(config = config, effective_settings = settings, library_paths = .libPaths()), file.path(outdir,'resolved_audit_config.rds'))
writeLines(capture.output(sessionInfo()), file.path(outdir,'sessionInfo.txt'))
writeLines(c('This is same-input implementation parity, not biological validation.',
  'Tolerance is a project audit setting and is not an official RCTD acceptance standard.',
  'No native Space Ranger target comparison is performed by this script.'), file.path(outdir,'AUDIT_SCOPE.txt'))

as_counts <- function(x, label) {
  if (!is.matrix(x) && !inherits(x,'Matrix')) stop(label, ' must be a matrix.')
  x <- methods::as(methods::as(Matrix::Matrix(x,sparse=TRUE),'generalMatrix'),'CsparseMatrix')
  x <- methods::as(x,'dMatrix')
  if (any(dim(x)==0L) || is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) || any(!nzchar(rownames(x))) ||
      any(!nzchar(colnames(x))) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) stop(label, ' requires unique gene and barcode names.')
  if (any(!is.finite(x@x)) || any(x@x<0) || any(abs(x@x-round(x@x))>1e-8)) stop(label,' requires raw integer counts.')
  Matrix::drop0(x)
}
hd <- readRDS(config$hd_file)
reference_input <- readRDS(config$reference_file)
if (!inherits(hd,'SpatialESSVisiumHD') || !all(c('counts','coordinates','metadata','all_positions','bin_size') %in% names(hd))) stop('hd_file must be the current prepared HD object with full all_positions.')
hd$counts <- as_counts(hd$counts,'HD counts')
if (!is.list(reference_input) || !all(c('counts','cell_types') %in% names(reference_input))) stop('reference_file requires list(counts,cell_types).')
reference_counts <- as_counts(reference_input$counts,'Reference counts')
cell_types <- reference_input$cell_types
if (is.null(names(cell_types)) || anyNA(names(cell_types)) || anyDuplicated(names(cell_types)) ||
    !setequal(names(cell_types),colnames(reference_counts))) stop('Reference labels must match all reference barcodes exactly.')
cell_types <- as.character(cell_types[colnames(reference_counts)])
if (anyNA(cell_types) || any(!nzchar(trimws(cell_types))) || any(grepl('/',cell_types,fixed=TRUE))) stop('Reference labels must be nonempty, nonmissing, and contain no slash.')
cell_types <- factor(cell_types); names(cell_types) <- colnames(reference_counts)
if (nlevels(cell_types)<2L) stop('At least two reference cell types are required.')
size <- number(hd$bin_size,'HD bin_size',2,TRUE)
if (!size %in% c(8L,16L,32L)) stop('This audit ROI uses input bin_size 8,16,or32 only.')
if (!is.null(config$bin_size) && as.numeric(config$bin_size)!=size) stop('Configured bin_size differs from the input; this script does not resample.')
ids <- colnames(hd$counts)
if (!identical(rownames(hd$coordinates),ids) || !identical(rownames(hd$metadata),ids)) stop('Input metadata/coordinate order differs from counts.')
positions <- as.data.frame(hd$all_positions,stringsAsFactors=FALSE)
if (!all(c('barcode','array_row','array_col','in_tissue') %in% names(positions)) || anyDuplicated(positions$barcode) || anyNA(positions$barcode)) stop('Full positions need unique barcode, array_row, array_col and in_tissue.')
for (key in c('array_row','array_col')) if (!is.numeric(positions[[key]]) || any(!is.finite(positions[[key]])) || any(positions[[key]]<0 | positions[[key]]!=floor(positions[[key]]))) stop('Positions require nonnegative integer array indices.')
if (anyNA(positions$in_tissue) || any(!positions$in_tissue %in% c(0L,1L))) stop('Full tissue flags must be zero/one.')
suffix <- function(x) {
  known <- grepl('-[0-9]+$',x); out <- rep('1',length(x))
  out[known] <- sub('^.*-([0-9]+)$','\\1',x[known]); out
}
if (length(unique(suffix(positions$barcode))) != 1L) stop('Audit one sample/barcode suffix at a time.')
if (anyDuplicated(paste(positions$array_row,positions$array_col,suffix(positions$barcode),sep=':'))) stop('Full positions contain duplicate source grid locations.')
raw_source_confirmed <- identical(hd$provenance$matrix_kind,'raw_feature_bc_matrix') && identical(hd$provenance$tissue_filtered,FALSE)
if (settings$require_raw_roi && !raw_source_confirmed) stop('Strict ROI audit needs raw_feature_bc_matrix with in_tissue=FALSE provenance. Set require_raw_roi=FALSE only for explicitly auditing an already prepared/filtered input.')

# Select a physical rectangle of complete 32um grid blocks, retaining every
# observed count column in the rectangle rather than taking head()/random bins.
roi <- named_list(config$roi32,'roi32')
if (length(setdiff(names(roi),c('row_start','col_start','n_rows','n_cols')))) stop('Unknown roi32 field.')
factor32 <- 32/size
positions$macro_row <- floor(positions$array_row/factor32)
positions$macro_col <- floor(positions$array_col/factor32)
row_bounds <- range(positions$macro_row); col_bounds <- range(positions$macro_col)
n_rows <- number(if(is.null(roi$n_rows)) min(5L,diff(row_bounds)+1L) else roi$n_rows,'roi32$n_rows',1,TRUE)
n_cols <- number(if(is.null(roi$n_cols)) min(5L,diff(col_bounds)+1L) else roi$n_cols,'roi32$n_cols',1,TRUE)
if (n_rows>diff(row_bounds)+1 || n_cols>diff(col_bounds)+1) stop('ROI dimensions exceed available positions.')
inside <- positions$in_tissue==1L
if (!any(inside)) stop('No tissue positions are available for choosing the ROI.')
choose_start <- function(values,bounds,width) max(bounds[1L],min(bounds[2L]-width+1L,floor(stats::median(values))-floor(width/2)))
row_start <- if(is.null(roi$row_start)) choose_start(positions$macro_row[inside],row_bounds,n_rows) else number(roi$row_start,'roi32$row_start',0,TRUE)
col_start <- if(is.null(roi$col_start)) choose_start(positions$macro_col[inside],col_bounds,n_cols) else number(roi$col_start,'roi32$col_start',0,TRUE)
in_roi <- positions$macro_row>=row_start & positions$macro_row<row_start+n_rows & positions$macro_col>=col_start & positions$macro_col<col_start+n_cols
selected_positions <- positions[in_roi,,drop=FALSE]
key <- paste(selected_positions$macro_row,selected_positions$macro_col,sep=':')
coverage <- table(key)
expected_positions <- factor32^2
if (length(coverage)!=n_rows*n_cols || any(coverage!=expected_positions)) stop('The chosen ROI lacks a complete source positions grid for one or more 32um blocks; select another aligned rectangle.')
count_index <- match(ids,positions$barcode)
if (anyNA(count_index)) stop('Count barcodes are missing from full positions.')
keep <- in_roi[count_index]
if (!any(keep)) stop('No observed count bins in the chosen ROI.')
if (sum(keep)>settings$max_input_bins) stop('ROI exceeds max_input_bins; reduce roi32 dimensions or explicitly increase the audit cap.')
hd$counts <- hd$counts[,keep,drop=FALSE]
hd$coordinates <- hd$coordinates[keep,,drop=FALSE]
hd$metadata <- hd$metadata[keep,,drop=FALSE]
hd$all_positions <- selected_positions[,setdiff(names(selected_positions),c('macro_row','macro_col')),drop=FALSE]
roi_info <- list(bin_size=size,row_start32=row_start,col_start32=col_start,n_rows32=n_rows,n_cols32=n_cols,
  aligned_um_bounds=c(row_min=row_start*32,row_max=(row_start+n_rows)*32,col_min=col_start*32,col_max=(col_start+n_cols)*32),
  position_rows=nrow(selected_positions),observed_count_bins=ncol(hd$counts),expected_positions_per_32um=expected_positions,
  raw_source_confirmed=raw_source_confirmed,feature_policy=hd$provenance$feature_policy,
  scope='Complete positions grid does not itself prove input counts are unfiltered or biologically accurate.')
saveRDS(roi_info,file.path(outdir,'roi_geometry.rds'))
saveRDS(hd,file.path(outdir,'audit_input_roi.rds'))
roi_table <- hd$metadata
roi_table$nUMI_all_input_genes <- as.numeric(Matrix::colSums(hd$counts))
utils::write.table(roi_table,file.path(outdir,'roi_bins.tsv'),sep='\t',quote=TRUE,row.names=FALSE)

# Compute the exact shared input, preserving full-library denominators before
# matching features. These are the same operations as the documented wrapper.
counts <- hd$counts
nUMI <- Matrix::colSums(counts)
nUMI_ref <- Matrix::colSums(reference_counts)
n_genes <- Matrix::colSums(counts>0)
eligible <- nUMI>=settings$min_umi & n_genes>=settings$min_genes & nUMI>0
if (!any(eligible)) stop('No ROI bins pass the configured QC; thresholds are user choices, not universal defaults.')
if (any(nUMI_ref<=0)) stop('Reference contains zero-library cells.')
common <- intersect(rownames(counts),rownames(reference_counts))
if (length(common)<10L) stop('Fewer than 10 shared genes; check identifiers and panel compatibility.')
saveRDS(list(spatial_all_gene_umi=nUMI,reference_all_gene_umi=nUMI_ref,
  shared_genes=common,eligible_barcodes=colnames(counts)[eligible],cell_types=cell_types),file.path(outdir,'shared_input_manifest.rds'))

capture_run <- function(expr,logfile) {
  con <- file(logfile,'wt'); sink(con); sink(con,type='message')
  on.exit({sink(type='message');sink();close(con)},add=TRUE)
  force(expr)
}
normalize_rows <- function(x) {
  x <- as.matrix(x)
  if (any(!is.finite(x)) || any(x<0)) stop('Original fit produced invalid weights.')
  totals <- rowSums(x)
  x[totals>0,] <- x[totals>0,,drop=FALSE]/totals[totals>0]
  x[totals<=0,] <- 0
  x
}
empty_summary <- function(mode) data.frame(mode=mode,input_bins=ncol(counts),eligible_bins=sum(eligible),
  original_fitted_bins=NA_integer_,wrapper_fitted_bins=NA_integer_,shared_genes=length(common),
  max_abs_weight_difference=NA_real_,max_abs_annotation_max_weight_difference=NA_real_,weight_tolerance=settings$tolerance,
  annotation_status_disagreements=NA_integer_,dominant_type_disagreements=NA_integer_,
  native_classification_equal=NA,native_pair_max_difference=NA_real_,
  wrapper_seconds=NA_real_,original_seconds=NA_real_,status='error',error='',stringsAsFactors=FALSE)
rows <- list()
for (mode in settings$modes) {
  mode_dir <- file.path(outdir,mode); dir.create(mode_dir,showWarnings=FALSE)
  row <- empty_summary(mode)
  problem <- tryCatch({
    t0 <- proc.time()[['elapsed']]
    fit <- capture_run(SpatialESS::run_visium_hd_rctd(hd,reference_counts,cell_types,
      doublet_mode=mode,max_cores=1L,min_umi=settings$min_umi,min_genes=settings$min_genes,
      create_args=settings$create_args,run_args=settings$run_args,seed=settings$seed,
      batch_size=settings$batch_size,keep_fit=FALSE),file.path(mode_dir,'wrapper.log'))
    row$wrapper_seconds <- proc.time()[['elapsed']]-t0
    saveRDS(fit,file.path(mode_dir,'wrapper_result.rds'))
    t0 <- proc.time()[['elapsed']]
    original <- capture_run(local({
      set.seed(settings$seed)
      ref <- spacexr::Reference(reference_counts[common,,drop=FALSE],cell_types,nUMI_ref)
      reference_audit <- list(retained_cells=ncol(ref@counts),retained_cells_per_type=table(ref@cell_types),
                              retained_cell_barcodes=colnames(ref@counts))
      puck <- spacexr::SpatialRNA(coords=as.data.frame(hd$coordinates[eligible,,drop=FALSE]),
                                  counts=counts[common,eligible,drop=FALSE],nUMI=nUMI[eligible])
      model <- do.call(spacexr::create.RCTD,c(list(spatialRNA=puck,reference=ref,max_cores=1L),settings$create_args))
      model@config$RCTDmode <- mode
      for (field in names(settings$run_args)) model@config[[field]] <- settings$run_args[[field]]
      model <- spacexr::run.RCTD(model,doublet_mode=mode)
      list(results=model@results,config=model@config,fitted_barcodes=colnames(model@spatialRNA@counts),
           reference_audit=reference_audit)
    }),file.path(mode_dir,'original.log'))
    row$original_seconds <- proc.time()[['elapsed']]-t0
    saveRDS(original,file.path(mode_dir,'original_results.rds'))
    fitted <- original$fitted_barcodes
    types <- colnames(original$results$weights)
    if (!length(fitted) || anyDuplicated(fitted) || !all(fitted %in% colnames(counts)) ||
        is.null(types) || anyDuplicated(types) || !setequal(types,colnames(fit$weights))) stop('Original/wrapper fitted barcode or type naming is invalid.')
    if (!is.null(fit$diagnostics$reference) &&
        fit$diagnostics$reference$retained_cells!=original$reference_audit$retained_cells)
      stop('Reference constructor retained different numbers of cells.')
    row$original_fitted_bins <- length(fitted)
    row$wrapper_fitted_bins <- fit$diagnostics$fitted_bins
    if (row$original_fitted_bins!=row$wrapper_fitted_bins) stop('Fitted bin counts differ.')
    if (!identical(rownames(fit$weights),colnames(counts))) stop('Wrapper output does not preserve every ROI barcode in input order.')
    if (mode=='full') {
      expected_fitted <- normalize_rows(original$results$weights[fitted,types,drop=FALSE])
      native_equal <- TRUE
    } else {
      native <- as.data.frame(original$results$results_df[fitted,,drop=FALSE])
      pair <- as.matrix(original$results$weights_doublet[fitted,c('first_type','second_type'),drop=FALSE])
      expected_fitted <- matrix(0,length(fitted),length(types),dimnames=list(fitted,types))
      # Independent implementation of the documented derived doublet policy.
      for (i in seq_along(fitted)) {
        label <- as.character(native$spot_class[[i]])
        if (label=='reject') next
        first <- as.character(native$first_type[[i]])
        if (label=='singlet') {
          expected_fitted[i,first] <- 1
        } else if (label %in% c('doublet_certain','doublet_uncertain')) {
          candidates <- c(first,as.character(native$second_type[[i]]))
          for (j in 1:2) if (is.finite(pair[i,j]) && pair[i,j]>0) {
            if (!candidates[j] %in% types) stop('Original doublet candidate does not match profile types.')
            expected_fitted[i,candidates[j]] <- expected_fitted[i,candidates[j]]+pair[i,j]
          }
        } else stop('Unexpected native doublet spot_class.')
      }
      expected_fitted <- normalize_rows(expected_fitted)
      classification <- fit$classification
      if (!is.list(classification) || !all(c('results_df','weights_doublet') %in% names(classification)) ||
          !setequal(rownames(classification$results_df),fitted) ||
          !setequal(rownames(classification$weights_doublet),fitted)) stop('Wrapper did not retain the native doublet outputs.')
      native_equal <- isTRUE(all.equal(classification$results_df[fitted,,drop=FALSE],native,
                                       tolerance=settings$tolerance,check.attributes=TRUE))
      delta_pair <- as.matrix(classification$weights_doublet[fitted,c('first_type','second_type'),drop=FALSE])-pair
      pair_missing_equal <- identical(is.na(as.matrix(classification$weights_doublet[fitted,c('first_type','second_type'),drop=FALSE])),is.na(pair))
      row$native_pair_max_difference <- if (any(is.finite(delta_pair))) max(abs(delta_pair[is.finite(delta_pair)])) else 0
      native_equal <- native_equal && pair_missing_equal && row$native_pair_max_difference<=settings$tolerance
      utils::write.table(data.frame(barcode=fitted,native,check.names=FALSE),file.path(mode_dir,'original_classification.tsv'),
                         sep='\t',quote=TRUE,row.names=FALSE)
    }
    row$native_classification_equal <- native_equal
    expected <- matrix(0,ncol(counts),length(types),dimnames=list(colnames(counts),types))
    expected[fitted,] <- expected_fitted
    observed <- as.matrix(fit$weights[colnames(counts),types,drop=FALSE])
    delta <- abs(observed-expected)
    row$max_abs_weight_difference <- max(delta)
    expected_status <- ifelse(eligible,'not_fit','filtered')
    expected_label <- rep(NA_character_,ncol(counts))
    for (i in seq_along(fitted)) {
      k <- match(fitted[i],colnames(counts)); w <- expected_fitted[i,]
      status <- 'ok'
      if (mode=='doublet') {
        label <- as.character(original$results$results_df[fitted[i],'spot_class'])
        if (label=='reject') status <- 'reject'
        if (label=='doublet_uncertain') status <- 'ambiguous'
      }
      maximum <- max(w)
      if (maximum<=0 && status=='ok') status <- 'not_fit'
      tied <- sum(abs(w-maximum)<=1e-10)>1L
      if (maximum>0 && tied && status=='ok') status <- 'ambiguous'
      if (maximum>0 && !tied && status=='ok') expected_label[k] <- types[which.max(w)]
      expected_status[k] <- status
    }
    annotations <- fit$annotations[match(colnames(counts),fit$annotations$barcode),,drop=FALSE]
    row$max_abs_annotation_max_weight_difference <- max(abs(annotations$max_weight-apply(expected,1L,max)))
    status_equal <- as.character(annotations$status)==expected_status
    label_equal <- (is.na(annotations$dominant_type) & is.na(expected_label)) |
      (!is.na(annotations$dominant_type) & !is.na(expected_label) & annotations$dominant_type==expected_label)
    row$annotation_status_disagreements <- sum(!status_equal)
    row$dominant_type_disagreements <- sum(!label_equal)
    compare <- data.frame(barcode=colnames(counts),eligible=eligible,
      max_abs_weight_difference=apply(delta,1L,max),wrapper_status=annotations$status,
      expected_status=expected_status,wrapper_dominant_type=annotations$dominant_type,
      expected_dominant_type=expected_label,status_equal=status_equal,dominant_type_equal=label_equal,
      stringsAsFactors=FALSE,row.names=NULL)
    utils::write.table(compare,file.path(mode_dir,'per_bin_comparison.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
    if (!is.finite(row$max_abs_weight_difference) || row$max_abs_weight_difference>settings$tolerance ||
        !is.finite(row$max_abs_annotation_max_weight_difference) || row$max_abs_annotation_max_weight_difference>settings$tolerance ||
        row$annotation_status_disagreements!=0L || row$dominant_type_disagreements!=0L || !native_equal)
      stop('Same-input parity did not meet the configured numerical/annotation checks; inspect per-bin and native outputs.')
    row$status <- 'pass'
    NULL
  },error=function(e)e)
  if (!is.null(problem)) {
    row$error <- conditionMessage(problem)
    writeLines(conditionMessage(problem),file.path(mode_dir,'ERROR.txt'))
  }
  rows[[mode]] <- row
  utils::write.table(do.call(rbind,rows),file.path(outdir,'audit_summary.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
  print(row)
  rm(list=intersect(c('fit','original','expected_fitted','expected','observed','delta'),ls()))
  invisible(gc())
}
result <- do.call(rbind,rows)
writeLines(capture.output(sessionInfo()),file.path(outdir,'sessionInfo.txt'))
if (any(result$status!='pass')) stop('Audit incomplete/failed; inspect audit_summary.tsv and per-mode ERROR/log files.',call.=FALSE)
cat('Same-input implementation audit passed for the configured ROI and modes.\n')
cat('This does not validate biological accuracy, cell-number fractions, or native Space Ranger bin equivalence.\n')
