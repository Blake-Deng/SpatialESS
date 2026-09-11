#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: summarize_bridge_regression.R VALIDATION_ROOT")
root <- normalizePath(args[1]); out <- file.path(root,"work/SpatialESS-SPARKLE/validation/bridge_interface")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
key <- function(x) paste(x$interaction_name,x$sender_group_name,x$receiver_group_name,sep="|")
metrics <- function(a,b) {
  ka<-key(a);kb<-key(b);stopifnot(!anyDuplicated(ka),!anyDuplicated(kb))
  common<-intersect(ka,kb); i<-match(common,ka);j<-match(common,kb)
  sa<-ka[a$pvalue<0.05];sb<-kb[b$pvalue<0.05]
  data.frame(active_a=nrow(a),active_b=nrow(b),
    active_jaccard=if(length(union(ka,kb)))length(common)/length(union(ka,kb)) else 1,
    max_abs_probability_difference=if(length(common))max(abs(a$probability[i]-b$probability[j])) else NA_real_,
    pvalue_exact_fraction=if(length(common))mean(a$pvalue[i]==b$pvalue[j]) else NA_real_,
    significant_a=length(sa),significant_b=length(sb),
    significance_disagreements=length(setdiff(sa,sb))+length(setdiff(sb,sa)),significance_rule="pvalue < 0.05")
}
same <- function(a,b,file) identical(readRDS(file.path(a,file)),readRDS(file.path(b,file)))
rows<-list();res<-list();hist<-list()
source <- Sys.getenv("OVARIAN_SOURCE","/home/dzf/cellchat_acceleration/SpatialESS_RNA_leakage_20260902")
for(d in list.dirs(file.path(root,"runs"),recursive=FALSE)) {
  old<-file.path(d,"old_bridge");new<-file.path(d,"new_bridge")
  for(mode in c("old_bridge","new_bridge","direct_main")) {
    p<-file.path(d,mode);s<-read.delim(file.path(p,"summary.tsv"));t<-readLines(file.path(p,"time.txt"))
    s$peak_rss_gib<-as.numeric(sub(".*: ","",grep("Maximum resident set size",t,value=TRUE)))/1024^2
    time<-strsplit(sub(".*\\): ","",grep("Elapsed ",t,value=TRUE)),":")[[1]]
    s$process_wall_seconds<-sum(rev(as.numeric(time))*60^(seq_along(time)-1))
    s$exit_code<-as.integer(sub(".*: ","",grep("Exit status",t,value=TRUE)))
    res[[length(res)+1L]]<-s
    if(mode=="old_bridge") next
    reference<-if(mode=="new_bridge")old else new
    m<-cbind(s[,c("condition","cells","mode")],metrics(readRDS(file.path(reference,"records.rds")),readRDS(file.path(p,"records.rds"))))
    m$records_identical<-same(reference,p,"records.rds")
    m$input_identical<-same(reference,p,"input.rds")
    m$contract_identical<-same(reference,p,"contract.rds")
    m$rng_identical<-same(reference,p,"rng.rds")
    m$pass<-m$records_identical & m$input_identical & m$contract_identical & m$rng_identical &
      m$active_jaccard==1 & m$max_abs_probability_difference==0 & m$pvalue_exact_fraction==1 & m$significance_disagreements==0
    rows[[length(rows)+1L]]<-m
    if(mode=="new_bridge" && s$cells==16247L) {
      paths<-list(historical_spatialess=if(s$condition=="sparkle") "results/spatialess_rctd_sparkle_full_20260904" else "results/spatialess_rctd_raw")
      if(s$condition=="sparkle") paths$archived_official_V3<-"results/official_rctd_sparkle_full_20260904"
      for(label in names(paths)) {
        q<-file.path(source,paths[[label]],"records.rds")
        h<-cbind(data.frame(condition=s$condition,cells=s$cells,comparator=label,source_records=q,archived_run=TRUE),
                  metrics(readRDS(q),readRDS(file.path(p,"records.rds"))))
        h$pass<-h$active_jaccard==1 & h$max_abs_probability_difference<=1e-15 & h$pvalue_exact_fraction==1 & h$significance_disagreements==0
        hist[[length(hist)+1L]]<-h
      }
    }
  }
}
comparisons<-do.call(rbind,rows);resources<-do.call(rbind,res);history<-do.call(rbind,hist)
for(name in c("comparisons","resources","history")) write.table(get(name),file.path(out,paste0(name,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
stopifnot(nrow(resources)==18L,all(resources$exit_code==0),nrow(comparisons)==12L,all(comparisons$pass))
cat("PASS:",nrow(resources),"runs;",nrow(comparisons),"strict bridge/version comparisons.\n")
print(history)
if(!all(history$pass)) warning("Historical preparation differences need explicit review; see history.tsv.")
