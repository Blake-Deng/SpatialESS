#!/usr/bin/env Rscript
# Local/server filesystem preparation only. No SSH or network requests.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==1L && args[1L] %in% c('--help','-h')) {
  cat('Usage: Rscript --vanilla spatialess_prepare_real_hd_roi.R RUN_DIR [SETTINGS_RDS]\n',
      'Defaults: 5x5 complete 32um blocks; >=80% tissue, observed source bins and QC-depth bins.\n',
      'Source8 depth QC: >=100 UMI and >=10 detected retained Gene Expression feature IDs.\n',
      'Settings may explicitly override roi_rows,roi_cols,min_tissue_fraction,min_observed_fraction,\n',
      'min_depth_fraction,min_umi,min_genes,seed,batch_size,overwrite. These are audit settings, not optimal thresholds.\n',
      'Reads RUN_DIR/inputs/official_outs and ovarian_reference.rds; writes only scoped ROI/config/report files.\n',sep='')
  quit(status=0L)
}
if (!length(args)%in%c(1L,2L)) stop('Expected RUN_DIR and optional settings RDS; use --help.')
run_dir <- normalizePath(args[1L],mustWork=TRUE)
suppressPackageStartupMessages({library(Matrix);library(SpatialESS);library(hdf5r);library(arrow)})
settings <- list(roi_rows=5L,roi_cols=5L,min_tissue_fraction=.8,min_observed_fraction=.8,
  min_depth_fraction=.8,min_umi=100,min_genes=10L,seed=1L,batch_size=64L,overwrite=FALSE)
if (length(args)==2L) {
  supplied <- readRDS(args[2L])
  if (!is.list(supplied) || is.null(names(supplied)) || anyDuplicated(names(supplied)) ||
      anyNA(names(supplied)) || length(setdiff(names(supplied),names(settings)))) stop('Invalid settings fields.')
  settings <- utils::modifyList(settings,supplied)
}
for (key in c('roi_rows','roi_cols','min_genes','seed','batch_size')) {
  x<-settings[[key]]
  if(!is.numeric(x)||length(x)!=1L||!is.finite(x)||x!=floor(x)||x<if(key=='seed')0 else 1) stop('Invalid ',key)
}
for (key in c('min_tissue_fraction','min_observed_fraction','min_depth_fraction')) {
  x<-settings[[key]];if(!is.numeric(x)||length(x)!=1L||!is.finite(x)||x<0||x>1)stop('Invalid ',key)
}
if (!is.numeric(settings$min_umi)||length(settings$min_umi)!=1L||!is.finite(settings$min_umi)||settings$min_umi<0) stop('Invalid min_umi.')
if (!is.logical(settings$overwrite)||length(settings$overwrite)!=1L||is.na(settings$overwrite))stop('Invalid overwrite.')
inputs <- file.path(run_dir,'inputs');reports<-file.path(run_dir,'reports')
dir.create(inputs,recursive=TRUE,showWarnings=FALSE);dir.create(reports,recursive=TRUE,showWarnings=FALSE)
owned <- c(file.path(inputs,paste0('roi',c(8,16,32),'.rds')),
  file.path(inputs,paste0('audit_bin',c(8,16,32),'_config.rds')),
  file.path(inputs,'real_roi_cli_config.rds'),file.path(reports,'roi_preparation.rds'))
if(any(file.exists(owned))&&!settings$overwrite)stop('ROI/config output exists; use a fresh RUN_DIR or explicit overwrite=TRUE.')
official <- file.path(inputs,'official_outs')
reference_file <- file.path(inputs,'ovarian_reference.rds')
if(!dir.exists(official)||!file.exists(reference_file))stop('Need official_outs directory/link and ovarian_reference.rds.')
find_square <- function(size) {
  name<-sprintf('square_%03dum',size)
  candidates<-unique(c(file.path(official,'binned_outputs',name),file.path(official,name),file.path(official,'outs','binned_outputs',name)))
  candidates<-candidates[dir.exists(candidates)]
  if(length(candidates)!=1L)stop('Expected one existing official ',name,' directory.')
  normalizePath(candidates)
}
paths <- lapply(c(8,16),function(size){
  square<-find_square(size)
  result<-list(size=size,square=square,raw=file.path(square,'raw_feature_bc_matrix.h5'),
    filtered=file.path(square,'filtered_feature_bc_matrix.h5'),positions=file.path(square,'spatial','tissue_positions.parquet'))
  if(!all(file.exists(unlist(result[c('raw','filtered','positions')]))))stop('Native ',size,' source files are incomplete.')
  result
});names(paths)<-c('8','16')
read_h5<-getFromNamespace('.visium_hd_read_h5','SpatialESS')
read_features<-getFromNamespace('.visium_hd_read_features','SpatialESS')
read_feature_table<-function(path)read_features(list(format='h5',path=path))
feature_tables<-lapply(paths,function(p)read_feature_table(p$filtered))
kept_ids<-lapply(feature_tables,function(f)f$id[f$type=='Gene Expression'])
if(!setequal(kept_ids[[1]],kept_ids[[2]])) {
  only8<-setdiff(kept_ids[[1]],kept_ids[[2]]);only16<-setdiff(kept_ids[[2]],kept_ids[[1]])
  differences<-rbind(data.frame(scope=rep('only_native8',length(only8)),feature_id=only8),
    data.frame(scope=rep('only_native16',length(only16)),feature_id=only16))
  utils::write.table(differences,file.path(reports,'feature_universe_difference.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
  stop('Native8/native16 filtered Gene Expression ID universes differ; no forced intersection performed.')
}
feature_manifest<-feature_tables[[1]][feature_tables[[1]]$id%in%kept_ids[[1]],,drop=FALSE]
utils::write.table(feature_manifest,file.path(reports,'roi_retained_feature_ids.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
# Fast metadata fingerprints are explicitly NOT whole-file checksums.
header_fingerprint<-function(path){
  info<-file.info(path);con<-file(path,'rb');on.exit(close(con))
  bytes<-readBin(con,'raw',n=min(65536,info$size));tmp<-tempfile('hd-header-')
  on.exit(unlink(tmp),add=TRUE);writeBin(bytes,tmp)
  data.frame(path=normalizePath(path),bytes=info$size,mtime=as.character(info$mtime),
    header_bytes=length(bytes),header_md5=unname(tools::md5sum(tmp)),
    checksum_scope='FIRST_64_KIB_ONLY_NOT_FULL_FILE_MD5',stringsAsFactors=FALSE)
}
source_files<-c(unlist(lapply(paths,function(p)unlist(p[c('raw','filtered','positions')]))),reference_file)
utils::write.table(do.call(rbind,lapply(source_files,header_fingerprint)),file.path(reports,'source_header_fingerprints.tsv'),
  sep='\t',quote=TRUE,row.names=FALSE)
writeLines('Large-file/archive full MD5 is not recomputed here; use the separate download/archive integrity report.',
  file.path(reports,'source_checksum_scope.txt'))
load_full<-function(p,keep){
  message('Reading raw native',p$size,' by filtered feature ID universe')
  counts<-read_h5(p$raw,gene_column='id',keep_feature_ids=keep)
  actual_ids<-attr(counts,'visium_hd_feature_ids',exact=TRUE);attr(counts,'visium_hd_feature_ids')<-NULL
  if(!setequal(actual_ids,keep)||!setequal(rownames(counts),keep))stop('Raw matrix does not retain the expected feature ID universe.')
  positions<-as.data.frame(arrow::read_parquet(p$positions))
  hd<-SpatialESS::prepare_visium_hd(counts,positions,bin_size=p$size,in_tissue=FALSE)
  hd$provenance$matrix_kind<-'raw_feature_bc_matrix';hd$provenance$tissue_filtered<-FALSE
  hd$provenance$feature_policy<-'filtered';hd$provenance$gene_column<-'id'
  hd$provenance$source_directory<-p$square;hd$provenance$matrix_source<-p$raw
  hd$provenance$source_feature_id_source<-p$filtered;hd$provenance$feature_ids_kept<-actual_ids
  hd$provenance$feature_id_count<-length(actual_ids)
  hd$provenance$feature_filter_order<-'feature ID selection before name mapping; ID rows retained in audit'
  hd$provenance$incomplete_count_source<-FALSE;hd$provenance$positions_file<-p$positions
  hd$provenance$data_origin<-'Official downloaded native Space Ranger outputs; archive integrity recorded separately'
  hd
}
# Integral-image window sums keep ROI search small despite full-slide positions.
window_sum<-function(m,nr,nc){
  a<-apply(m,2,cumsum);a<-t(apply(a,1,cumsum))
  p<-matrix(0,nrow(m)+1L,ncol(m)+1L);p[-1L,-1L]<-a
  rr<-seq_len(nrow(m)-nr+1L);cc<-seq_len(ncol(m)-nc+1L)
  p[rr+nr,cc+nc,drop=FALSE]-p[rr,cc+nc,drop=FALSE]-p[rr+nr,cc,drop=FALSE]+p[rr,cc,drop=FALSE]
}
full8<-load_full(paths[['8']],kept_ids[['8']])
pos<-full8$all_positions;factor32<-4L
mr<-floor(pos$array_row/factor32);mc<-floor(pos$array_col/factor32)
nr<-max(mr)+1L;nc<-max(mc)+1L;nb<-nr*nc
if(settings$roi_rows>nr||settings$roi_cols>nc)stop('Requested ROI exceeds available capture grid.')
index<-mr+1L+mc*nr
count_index<-match(colnames(full8$counts),pos$barcode)
if(anyNA(count_index))stop('Raw count columns missing complete positions.')
umi<-Matrix::colSums(full8$counts);ngenes<-diff(full8$counts@p)
eligible<-umi>=settings$min_umi & ngenes>=settings$min_genes & umi>0
as_grid<-function(v)matrix(v,nrow=nr,ncol=nc)
position_grid<-as_grid(tabulate(index,nbins=nb))
tissue_grid<-as_grid(tabulate(index[pos$in_tissue==1L],nbins=nb))
observed_grid<-as_grid(tabulate(index[count_index],nbins=nb))
depth_grid<-as_grid(tabulate(index[count_index[eligible]],nbins=nb))
position_windows<-window_sum(position_grid,settings$roi_rows,settings$roi_cols)
tissue_windows<-window_sum(tissue_grid,settings$roi_rows,settings$roi_cols)
observed_windows<-window_sum(observed_grid,settings$roi_rows,settings$roi_cols)
depth_windows<-window_sum(depth_grid,settings$roi_rows,settings$roi_cols)
expected<-16L*settings$roi_rows*settings$roi_cols
valid<-position_windows==expected & tissue_windows/expected>=settings$min_tissue_fraction &
  observed_windows/expected>=settings$min_observed_fraction & depth_windows/expected>=settings$min_depth_fraction
candidate<-which(matrix(TRUE,nrow(valid),ncol(valid)),arr.ind=TRUE)
selection<-data.frame(row_start32=candidate[,1L]-1L,col_start32=candidate[,2L]-1L,
  position_count=as.numeric(position_windows),tissue_fraction=as.numeric(tissue_windows)/expected,
  observed_fraction=as.numeric(observed_windows)/expected,depth_fraction=as.numeric(depth_windows)/expected,
  eligible=as.logical(valid))
tissue<-pos$in_tissue==1L
center<-c(row=stats::median(pos$array_row[tissue]/4),col=stats::median(pos$array_col[tissue]/4))
selection$distance_to_tissue_median_um<-32*sqrt((selection$row_start32+settings$roi_rows/2-center['row'])^2+
  (selection$col_start32+settings$roi_cols/2-center['col'])^2)
ranking<-order(!selection$eligible,selection$distance_to_tissue_median_um,-selection$depth_fraction,
  -selection$tissue_fraction,selection$row_start32,selection$col_start32)
utils::write.table(head(selection[ranking,],100L),file.path(reports,'roi_selection_candidates_top100.tsv'),
  sep='\t',quote=TRUE,row.names=FALSE)
if(!any(selection$eligible))stop('No ROI meets the explicitly configured tissue/coverage/depth criteria; inspect candidate report before choosing new settings.')
chosen<-selection[ranking[1L],,drop=FALSE]
bounds<-list(row_start32=chosen$row_start32,col_start32=chosen$col_start32,
  n_rows32=settings$roi_rows,n_cols32=settings$roi_cols,
  row_min_um=chosen$row_start32*32,row_max_um=(chosen$row_start32+settings$roi_rows)*32,
  col_min_um=chosen$col_start32*32,col_max_um=(chosen$col_start32+settings$roi_cols)*32)
crop_hd<-function(hd,bounds){
  f<-32/hd$bin_size
  inside<-function(p) p$array_row>=bounds$row_start32*f & p$array_row<(bounds$row_start32+bounds$n_rows32)*f &
    p$array_col>=bounds$col_start32*f & p$array_col<(bounds$col_start32+bounds$n_cols32)*f
  k<-inside(hd$metadata);kp<-inside(hd$all_positions)
  if(!any(k)||sum(kp)!=bounds$n_rows32*bounds$n_cols32*f^2)stop('ROI does not have complete positions at native resolution.')
  hd$counts<-hd$counts[,k,drop=FALSE];hd$coordinates<-hd$coordinates[k,,drop=FALSE]
  hd$metadata<-hd$metadata[k,,drop=FALSE];hd$all_positions<-hd$all_positions[kp,,drop=FALSE]
  hd$provenance$roi_bounds<-bounds;hd$provenance$retained_bins<-ncol(hd$counts)
  hd$provenance$roi_scope<-'Complete 32um-aligned physical rectangle; all native positions retained; raw barcode counts with filtered feature IDs'
  hd
}
hd8<-crop_hd(full8,bounds)
saveRDS(hd8,file.path(inputs,'roi8.rds'))
rm(full8,pos,mr,mc,index,count_index,umi,ngenes,eligible,position_grid,tissue_grid,observed_grid,depth_grid,
   position_windows,tissue_windows,observed_windows,depth_windows,selection,candidate,valid);invisible(gc())
full16<-load_full(paths[['16']],kept_ids[['16']])
hd16<-crop_hd(full16,bounds);rm(full16);invisible(gc())
saveRDS(hd16,file.path(inputs,'roi16.rds'))
check_equal<-function(native,rebuilt,label,native_evidence){
  row<-data.frame(check=label,evidence_scope=native_evidence,
    genes_equal=setequal(rownames(native$counts),rownames(rebuilt$counts)),
    count_barcodes_equal=setequal(colnames(native$counts),colnames(rebuilt$counts)),
    counts_exact=FALSE,max_abs_count_difference=NA_real_,
    position_barcodes_equal=setequal(native$all_positions$barcode,rebuilt$all_positions$barcode),
    tissue_mask_equal=FALSE,status='failed',stringsAsFactors=FALSE)
  if(row$genes_equal && row$count_barcodes_equal){
    d<-native$counts-rebuilt$counts[rownames(native$counts),colnames(native$counts),drop=FALSE]
    row$max_abs_count_difference<-if(length(d@x))max(abs(d@x))else 0
    row$counts_exact<-row$max_abs_count_difference==0
  }
  if(row$position_barcodes_equal)row$tissue_mask_equal<-identical(as.integer(native$all_positions$in_tissue),
    as.integer(rebuilt$all_positions$in_tissue[match(native$all_positions$barcode,rebuilt$all_positions$barcode)]))
  if(all(row$genes_equal,row$count_barcodes_equal,row$counts_exact,row$position_barcodes_equal,row$tissue_mask_equal))row$status<-'pass'
  row
}
rebuilt16<-SpatialESS::aggregate_visium_hd(hd8,16)
check16<-check_equal(hd16,rebuilt16,'official_native16_vs_raw8_aggregation_ROI',
  'REAL_OFFICIAL_NATIVE16_COMPARISON_WITHIN_SELECTED_ROI_ONLY')
utils::write.table(check16,file.path(reports,'real_roi_native16_equivalence.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
if(check16$status!='pass')stop('Real ROI native16 equivalence failed; inspect report. No native16 gene/column forcing was performed.')
hd32<-SpatialESS::aggregate_visium_hd(hd16,32)
from8_32<-SpatialESS::aggregate_visium_hd(hd8,32)
check32<-check_equal(hd32,from8_32,'raw8_to32_vs_native16_to32_ROI',
  'REAL_DATA_CROSS_ROUTE_ASSOCIATIVITY_ONLY_NOT_OFFICIAL_NATIVE32')
utils::write.table(check32,file.path(reports,'real_roi_32_crossroute.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
if(check32$status!='pass')stop('ROI 32 cross-route equality failed.')
hd32$provenance$native_output_equivalence<-'No official native32 supplied; cross-route 8->32 versus16->32 checked only'
saveRDS(hd32,file.path(inputs,'roi32.rds'))
writeLines('No official square_032um was supplied. The32 comparison uses real data through two aggregation routes; it is NOT native Space Ranger32 equivalence.',
  file.path(reports,'NATIVE32_NOT_AVAILABLE.txt'))
# Derived MEX is for running the production CLI on this exact local ROI; it is
# never labelled an independently generated Space Ranger native reference.
write_mex<-function(hd,feature_table,square){
  dir.create(file.path(square,'spatial'),recursive=TRUE,showWarnings=FALSE)
  info<-feature_table[match(rownames(hd$counts),feature_table$id),c('id','name','type'),drop=FALSE]
  if(anyNA(info$id))stop('MEX feature metadata lookup failed.')
  for(kind in c('raw_feature_bc_matrix','filtered_feature_bc_matrix')){
    d<-file.path(square,kind);dir.create(d,recursive=TRUE,showWarnings=FALSE)
    k<-if(kind=='filtered_feature_bc_matrix')hd$metadata$in_tissue==1L else rep(TRUE,ncol(hd$counts))
    if(!any(k))stop('No tissue count bins in ROI for CLI export.')
    Matrix::writeMM(hd$counts[,k,drop=FALSE],file.path(d,'matrix.mtx'))
    utils::write.table(info,file.path(d,'features.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
    writeLines(colnames(hd$counts)[k],file.path(d,'barcodes.tsv'))
  }
  utils::write.csv(hd$all_positions,file.path(square,'spatial','tissue_positions.csv'),row.names=FALSE)
}
roi_root<-file.path(inputs,'real_roi','binned_outputs')
write_mex(hd8,feature_tables[['8']],file.path(roi_root,'square_008um'))
write_mex(hd16,feature_tables[['16']],file.path(roi_root,'square_016um'))
writeLines('Derived ROI export from official native8/native16 raw sources with filtered feature IDs; not independently regenerated Space Ranger output.',
  file.path(inputs,'real_roi','DERIVED_ROI_SCOPE.txt'))
reference<-readRDS(reference_file)
if(!is.list(reference)||!all(c('counts','cell_types')%in%names(reference)))stop('Invalid prepared reference.')
if(!is.null(names(reference$cell_types)) && !setequal(names(reference$cell_types),colnames(reference$counts)))stop('Reference labels/barcodes differ.')
objects<-list('8'=hd8,'16'=hd16,'32'=hd32)
summary_rows<-list();fit_by_bin<-list()
for(key in names(objects)){
  hd<-objects[[key]];u<-Matrix::colSums(hd$counts);g<-diff(hd$counts@p)
  e<-u>=settings$min_umi & g>=settings$min_genes & u>0
  if(!any(e))stop('No eligible bins at resolution',key)
  sigma_min<-if(mean(u[e]>300)>=.8)300 else 100
  if(!any(u[e]>sigma_min))stop('No bins exceed the explicit sigma UMI threshold at bin',key)
  create<-list(UMI_min_sigma=sigma_min)
  audit_config<-list(hd_file=file.path(inputs,paste0('roi',key,'.rds')),reference_file=reference_file,
    output_dir=file.path(run_dir,paste0('audit_bin',key)),bin_size=as.numeric(key),modes=c('full','doublet'),max_cores=1L,
    min_umi=settings$min_umi,min_genes=settings$min_genes,seed=settings$seed,batch_size=settings$batch_size,
    create_args=create,roi32=list(row_start=bounds$row_start32,col_start=bounds$col_start32,
      n_rows=settings$roi_rows,n_cols=settings$roi_cols),max_input_bins=max(500L,ncol(hd$counts)),
    require_raw_roi=TRUE,overwrite=FALSE,tolerance=1e-8)
  saveRDS(audit_config,file.path(inputs,paste0('audit_bin',key,'_config.rds')))
  fit_by_bin[[paste0('bin',key)]]<-list(create_args=create)
  summary_rows[[key]]<-data.frame(bin_size_um=as.numeric(key),count_bins=ncol(hd$counts),position_bins=nrow(hd$all_positions),
    tissue_count_bins=sum(hd$metadata$in_tissue==1L),retained_gene_ids=nrow(hd$counts),total_umi=sum(hd$counts),
    qc_eligible=sum(e),min_umi=settings$min_umi,min_genes=settings$min_genes,UMI_min_sigma=sigma_min,
    sigma_eligible=sum(u[e]>sigma_min),shared_reference_gene_ids=length(intersect(rownames(hd$counts),rownames(reference$counts))),
    reference_cells=ncol(reference$counts),reference_types=length(unique(as.character(reference$cell_types))),
    stringsAsFactors=FALSE)
}
main_config<-list(data_dir=file.path(inputs,'real_roi'),reference_file=reference_file,
  output_dir=file.path(run_dir,'real_roi_cli_results'),bin_sizes=c(8,16,32),doublet_mode='full',max_cores=1L,
  min_umi=settings$min_umi,min_genes=settings$min_genes,seed=settings$seed,batch_size=settings$batch_size,
  fit_args_by_bin=fit_by_bin,read_args=list(gene_column='id',feature_policy='filtered'),save_inputs=TRUE,overwrite=FALSE)
saveRDS(main_config,file.path(inputs,'real_roi_cli_config.rds'))
summary<-do.call(rbind,summary_rows)
utils::write.table(summary,file.path(reports,'real_roi_input_summary.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
utils::write.table(chosen,file.path(reports,'selected_roi_quality.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
preparation<-list(settings=settings,bounds=bounds,selected_quality=chosen,source_paths=paths,
  feature_ids_kept=kept_ids,reference_file=reference_file,input_summary=summary,
  native16_roi_check=check16,crossroute32_check=check32,
  sigma_threshold_rule='300 if >=80% eligible bins exceed300UMI; otherwise100. Explicit validation setting, not an optimized biological threshold.',
  native32_available=FALSE,scope='Real native16 evidence within selected ROI plus real-data crossroute32, not full-slice or genuine native32 validation.')
saveRDS(preparation,file.path(reports,'roi_preparation.rds'))
writeLines(capture.output(sessionInfo()),file.path(reports,'roi_preparation_sessionInfo.txt'))
print(chosen);print(summary);print(check16);print(check32)
cat('Prepared ROI inputs, strict same-input audit configs, and main CLI config under ',run_dir,'\n',sep='')
