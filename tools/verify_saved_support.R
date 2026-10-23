args <- commandArgs(trailingOnly=TRUE);root<-args[1];out<-args[2]
x<-readRDS(file.path(root,'paper/cohort_sensitivity/data/prepared_slice_scores.rds'))
s<-x$score;m<-x$sample_metadata
stopifnot(all(dim(s)==c(45L,50040L)),all(s>=0),!anyNA(s),identical(rownames(s),as.character(m$sample_id)))
d<-read.delim(gzfile(file.path(root,'paper/submission_materials/Supplementary_Data_4_Lung_Support.tsv.gz')),check.names=FALSE)
d<-d[match(colnames(s),d$feature_id),];stopifnot(identical(colnames(s),as.character(d$feature_id)))
audit<-data.frame(feature_id=colnames(s),check.names=FALSE)
for(cond in c('Control','Disease')) {
 idx<-which(m$condition==cond);pm<-split(idx,m$patient_id[idx])
 anyscore<-sapply(pm,function(ii) colSums(s[ii,,drop=FALSE]!=0)>0)
 allscore<-sapply(pm,function(ii) colSums(s[ii,,drop=FALSE]!=0)==length(ii))
 prefix<-if(cond=='Control') 'control' else 'pf'
 audit[[paste0(prefix,'_nonzero_slices')]]<-colSums(s[idx,,drop=FALSE]!=0)
 audit[[paste0(prefix,'_nonzero_patients')]]<-rowSums(anyscore)
 audit[[paste0(prefix,'_all_available_sections_nonzero_patients')]]<-rowSums(allscore)
 audit[[paste0(prefix,'_single_section_patients')]]<-sum(lengths(pm)==1L)
 stopifnot(identical(as.numeric(audit[[paste0(prefix,'_nonzero_slices')]]),as.numeric(d[[paste0(prefix,'_nonzero_slices')]])),identical(as.numeric(audit[[paste0(prefix,'_nonzero_patients')]]),as.numeric(d[[paste0(prefix,'_nonzero_patients')]])))
}
write.table(audit,gzfile(file.path(out,'DESCRIPTIVE_SUPPORT_AUDIT.tsv.gz'),'wt'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PASS 50040 features: nonzero slices and at-least-one-section patient counts match Data4; all-available-section counts include single-section patients. No model or BH recalculation.\n')
