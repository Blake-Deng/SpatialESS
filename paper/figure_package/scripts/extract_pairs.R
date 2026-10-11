# Export actual saved binary probabilities at 17-digit precision, without inference.
args <- commandArgs(TRUE)
root <- normalizePath(args[1])
project <- if(length(args)>1L) normalizePath(args[2]) else dirname(dirname(root))
release <- file.path(project,'releases/SpatialESS_v0.1.4_20261001/repositories/SpatialESS')
manifest <- read.delim(file.path(root,'provenance/build_inputs.tsv'),stringsAsFactors=FALSE)
cases <- c('cosmx_50000',
  paste0('cervical_',rep(c('center','stratified'),each=2),'_',rep(c(1000,5000),2),'_default'))
extract <- function(path) {
  d <- readRDS(path)
  if(is.list(d)&&!is.data.frame(d))d<-d$records
  stopifnot(is.data.frame(d)); d
}
key <- function(d)do.call(paste,c(d[c('sender_group_name','receiver_group_name','interaction_name')],sep='\r'))
checks <- list()
for(case in cases) {
  row <- manifest[manifest$case==case,,drop=FALSE];stopifnot(nrow(row)==1L)
  official <- extract(row$official)
  ess <- extract(file.path(release,'results/records',paste0(case,'.rds')))
  ko<-key(official);ke<-key(ess)
  stopifnot(!anyDuplicated(ko),!anyDuplicated(ke),setequal(ko,ke))
  ix<-match(ke,ko)
  d<-ess[c('sender_group_name','receiver_group_name','interaction_name')]
  d$probability_official<-official$probability[ix];d$probability_ess<-ess$probability
  d$pvalue_official<-official$pvalue[ix];d$pvalue_ess<-ess$pvalue;d$case<-case
  checks[[case]]<-data.frame(case=case,n=nrow(d),unmatched=0L,
    zeros=sum(d$probability_official==0|d$probability_ess==0),
    nonfinite=sum(!is.finite(d$probability_official)|!is.finite(d$probability_ess)),
    max_abs_difference=max(abs(d$probability_official-d$probability_ess)),
    pearson=cor(d$probability_official,d$probability_ess),
    exact_p_fraction=mean(d$pvalue_official==d$pvalue_ess),
    official_file=row$official,official_md5=unname(tools::md5sum(row$official)))
  for(n in names(d))if(is.double(d[[n]]))d[[n]]<-sprintf('%.17g',d[[n]])
  write.table(d,file.path(root,'source_data/current',paste0(case,'_paired.tsv')),
    sep='\t',row.names=FALSE,quote=TRUE)
}
write.table(do.call(rbind,checks),file.path(root,'reports/pairing_checks.tsv'),
  sep='\t',row.names=FALSE,quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(root,'logs/extraction_sessionInfo.txt'))
cat('Exported',length(cases),'binary-result pair tables. No inference rerun.\n')
