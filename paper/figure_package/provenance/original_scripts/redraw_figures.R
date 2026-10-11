#!/usr/bin/env Rscript
# Frozen numerical tables are read-only; all changes concern figure presentation.
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(cowplot); library(ggrepel)
  library(scales); library(ComplexHeatmap); library(circlize); library(grid)
})
args <- commandArgs(TRUE)
root <- normalizePath(if(length(args)) args[1] else file.path(getwd()))
out <- file.path(root,'figures'); datadir <- file.path(root,'source_data')
dir.create(file.path(out,'panels'),recursive=TRUE,showWarnings=FALSE)
rd <- function(n) {
  con <- gzfile(file.path(datadir,paste0(n,'.tsv.gz')),'rt'); on.exit(close(con))
  as.data.table(read.delim(con,check.names=FALSE,stringsAsFactors=FALSE))
}
red <- '#C53F50'; blue <- '#2166AC'; teal <- blue; gold <- '#7B858F'
neutral <- '#515C64'; ink <- '#111820'; pale <- '#E5EBEE'; muted <- '#85949D'
method_cols <- c('SpatialESS'=blue,'SpatialCellChat V3'=red)
font <- 'Arial'
stopifnot(any(systemfonts::system_fonts()$family==font))
do.call(grDevices::pdfFonts,setNames(list(grDevices::pdfFonts('Helvetica')[[1]]),font))
theme_set(theme_classic(base_size=8,base_family=font)+theme(
  text=element_text(colour=ink),axis.line=element_line(linewidth=.34,colour=ink),
  axis.ticks=element_line(linewidth=.25,colour=neutral),axis.ticks.length=unit(.9,'mm'),
  axis.text=element_text(size=7,colour=ink),axis.title=element_text(size=7.6),
  legend.text=element_text(size=7,colour=ink),legend.title=element_text(size=7,colour=ink),
  legend.key.size=unit(2.8,'mm'),legend.spacing.x=unit(.8,'mm'),
  plot.margin=margin(2,4,2,3),plot.title=element_text(size=8,face='bold'),
  strip.background=element_blank(),strip.text=element_text(size=7,face='bold'),
  panel.spacing=unit(2,'mm')))
num <- function(x,digits=2) trimws(ifelse(is.na(x),'NA',ifelse(x==0,'0',
  ifelse(abs(x)<.001,formatC(x,format='e',digits=digits),formatC(x,format='fg',digits=digits+1)))))
count <- function(x) format(round(x),big.mark=',',scientific=FALSE,trim=TRUE)
cell_lab <- function(x) vapply(x,function(z) if(z<1000)as.character(z) else
  paste0(format(z/1000,trim=TRUE,scientific=FALSE),'k'),character(1))
effect_scale <- function(lim) scale_fill_gradient2(low=blue,mid='white',high=red,
  midpoint=0,limits=c(-lim,lim),oob=squish,na.value=pale)
panel <- function(p,letter,title,note=NULL) {
  ggdraw()+draw_plot(p,0,if(is.null(note))0 else .065,1,if(is.null(note)).91 else .835)+
    draw_label(letter,.004,.97,hjust=0,vjust=1,fontfamily=font,fontface='bold',size=10)+
    draw_label(title,.065,.982,hjust=0,vjust=1,fontfamily=font,fontface='bold',size=8)+
    {if(!is.null(note))draw_label(note,.02,.012,hjust=0,vjust=0,fontfamily=font,size=6.5,colour=neutral)}
}
artwork <- list()
export <- function(p,name,w=180,h=180) {
  artwork[[name]] <<- list(plot=p,width=w,height=h)
  f<-file.path(out,name)
  cairo_pdf(paste0(f,'.pdf'),width=w/25.4,height=h/25.4,family=font,onefile=FALSE)
  print(p);dev.off()
  svglite::svglite(paste0(f,'.svg'),width=w/25.4,height=h/25.4,system_fonts=list(sans=font))
  print(p);dev.off()
  for(dpi in c(600,180)) {
    ragg::agg_png(paste0(f,if(dpi==180)'_preview'else'', '.png'),width=w,height=h,units='mm',res=dpi,background='white')
    print(p);dev.off()
  }
  message('Exported ',name,' at ',w,' x ',h,' mm, Arial, PDF/SVG/600 dpi PNG')
}
matrix_plot <- function(d,row_levels,col_levels,value='label',fill='tint',base=7) {
  d<-copy(d);d[,row:=factor(row,levels=rev(row_levels))];d[,col:=factor(col,levels=col_levels)]
  ggplot(d,aes(col,row))+geom_tile(aes(fill=.data[[fill]]),width=.96,height=.88)+
    geom_text(aes(label=.data[[value]]),family=font,size=base/2.845)+scale_fill_identity()+
    scale_x_discrete(position='top',expand=expansion(add=.1))+
    scale_y_discrete(expand=expansion(add=.1))+labs(x=NULL,y=NULL)+
    theme_minimal(base_family=font,base_size=8)+theme(panel.grid=element_blank(),
      axis.text.x=element_text(size=6.8,colour=ink),axis.text.y=element_text(size=7,colour=ink),
      plot.margin=margin(1,3,1,2))
}
density_agreement <- function(d,xcol,ycol,xlab,ylab,annotation,bins=45) {
  d<-copy(d)[is.finite(get(xcol))&is.finite(get(ycol))&get(xcol)>0&get(ycol)>0]
  d[,`:=`(x=log10(get(xcol)),y=log10(get(ycol)))];rg<-range(d$x,d$y);pad<-diff(rg)*.07
  ggplot(d,aes(x,y))+geom_hex(bins=bins,linewidth=0)+
    scale_fill_gradient(low='#80B4DE',high='#164F86',trans='log10',name='Records')+
    geom_abline(intercept=0,slope=1,linewidth=.3,linetype='dashed',colour=ink)+
    annotate('text',x=rg[1]+.025*diff(rg),y=rg[2]+pad*.3,label=annotation,
      hjust=0,vjust=1,family=font,size=2.4,lineheight=1.1)+
    coord_fixed(xlim=rg+c(-pad,pad),ylim=rg+c(-pad,pad))+
    labs(x=xlab,y=ylab)+guides(fill=guide_colourbar(display='rectangles'))+
    theme(legend.position='right',legend.key.height=unit(4,'mm'),legend.key.width=unit(2,'mm'))
}
prov<-list()
provenance<-function(fig,panel,source,columns,rule,n) {
  prov[[length(prov)+1L]]<<-data.table(figure=fig,panel=panel,source_file=source,
    columns=columns,filtering_and_statistical_rule=rule,n_observations=n)
}

# The supplied illustrated PDF is inserted by compose_original_workflow.py.
provenance('Figure1','A','source_snapshot/SpatialESS_workflow_comparison.pdf','workflow concepts only',
  'Original vector workflow; vertical spacing compacted; colors strengthened; no numerical observations',NA_integer_)

sc<-rd('cosmx_scaling');ac<-rd('cosmx_fidelity');cp<-rd('cosmx_probability_records');nc<-max(ac$cells)
f1b<-panel(density_agreement(cp,'probability_v3','probability_ess','V3, log10(probability)',
  'SpatialESS\nlog10(probability)',paste0('r = ',num(ac[cells==nc,probability_pearson]),
  '\nMax err. ',num(ac[cells==nc,max_abs_probability_difference]),'\nn = ',count(nrow(cp)))),
  'B',paste0('Probability fidelity (',cell_lab(nc),')'))
metrics<-c('Pearson\nr','Max |prob.\ndifference|','Exact p\nfraction','Active\nJaccard','Call\ndisagreements')
vals<-c('probability_pearson','max_abs_probability_difference','exact_pvalue_fraction','active_jaccard','significance_threshold_disagreements')
fmat<-rbindlist(lapply(seq_along(vals),function(i)data.table(row=cell_lab(ac$cells),col=metrics[i],
  label=num(ac[[vals[i]]]),tint=if(i==2)'#DE98A2'else'#8BB3D9')))
f1c<-panel(matrix_plot(fmat,cell_lab(ac$cells),metrics,base=7.2),'C','Agreement across paired scales',
  'Matched cells, LR, parameters and seed')
provenance('Figure1','B','source_data/cosmx_probability_records.tsv.gz; source_data/cosmx_fidelity.tsv.gz',
  'probability_v3;probability_ess;cells;probability_pearson;max_abs_probability_difference',
  'Largest paired scale; positive probabilities, log10 hexbin; summary error from binary-result audit',nrow(cp))
provenance('Figure1','C','source_data/cosmx_fidelity.tsv.gz',paste(c('cells',vals),collapse=';'),
  'All six completed paired scales; calls use unadjusted permutation p <= 0.05',nrow(ac))
longsc<-rbindlist(lapply(c('spatialess','official'),function(m)data.table(cells=sc$cells,
  method=if(m=='spatialess')'SpatialESS'else'SpatialCellChat V3',status=sc[[paste0(m,'_status')]],
  seconds=sc[[paste0(m,'_wall_seconds')]],gb=sc[[paste0(m,'_peak_rss_gib')]]*2^30/1e9)))
scalep<-function(value,ylab) {
  d<-longsc[status=='completed'];s50<-sc[cells==nc]
  ann<-if(value=='seconds')paste0(cell_lab(nc),': ',count(s50$official_wall_seconds/s50$spatialess_wall_seconds),'-fold faster')else
    paste0(cell_lab(nc),': ',formatC(100*(1-s50$spatialess_peak_rss_gib/s50$official_peak_rss_gib),format='f',digits=1),'% less peak RSS')
  p<-ggplot(d,aes(cells,.data[[value]],colour=method,group=method))+
    geom_line(linewidth=.85)+geom_point(size=2.05,stroke=.35)+scale_color_manual(values=method_cols)+
    scale_x_log10(breaks=c(1000,10000,100000),labels=cell_lab,limits=range(sc$cells)*c(.83,1.2))+
    scale_y_log10(labels=label_number(big.mark=','))+labs(x='Cells (log scale)',y=ylab,colour=NULL)+
    theme(legend.position='top',legend.margin=margin(-3,0,-3,0),legend.key.size=unit(2.4,'mm'))+
    annotate('text',x=min(sc$cells)*1.08,y=max(d[[value]])*.6,label=ann,hjust=0,colour=blue,family=font,size=2.6)
  fail<-sc[official_status!='completed',.(cells,reason=ifelse(official_status=='failed_signal_11','sig. 11','alloc.'))]
  band<-ggplot(fail,aes(cells,1))+geom_point(shape=4,size=1.5,stroke=.6,colour=neutral)+
    geom_text(aes(label=reason),y=.37,size=2.1,family=font,colour=neutral)+
    annotate('text',x=min(sc$cells),y=.8,label='V3 failed',hjust=0,size=2.3,family=font,colour=neutral)+
    scale_x_log10(limits=range(sc$cells)*c(.83,1.2))+ylim(0,1.6)+theme_void()+theme(plot.margin=margin(0,4,0,3))
  plot_grid(band,p,ncol=1,rel_heights=c(.14,1),align='v',axis='lr')
}
f1d<-panel(scalep('seconds','Wall-clock time (s)'),'D','Runtime scaling')
f1e<-panel(scalep('gb','Peak RSS (GB)'),'E','Measured process memory')
for(z in c('D','E'))provenance('Figure1',z,'source_data/cosmx_scaling.tsv.gz',
  if(z=='D')'cells;*_status;*_wall_seconds;official_failure_reason'else'cells;*_status;*_peak_rss_gib;official_failure_reason',
  'Completed observations only in curves; failure marks have no numeric resource ordinate; GiB converted to decimal GB',
  nrow(longsc[status=='completed']))

# The former S3 evidence is now part of the main computational figure.
spfile<-'spatial_sparsity_validation_20260922/data/spatial_sparsity_summary.tsv'
sp<-fread(file.path(root,spfile),integer64='double');stopifnot(nrow(sp)==14L,all(sp$validation_status=='PASS'))
stopifnot(all(sp$N_squared==as.double(sp$N)^2),all(abs(sp$mean_out_degree-sp$E/sp$N)<1e-8),
  all(abs(sp$graph_density-sp$E/sp$N_squared)<1e-12),
  all(sp$csr_total_bytes==sp$offset_bytes+sp$neighbor_bytes+sp$distance_bytes),
  all(sp$dense_float64_bytes==8*sp$N_squared))
series<-c('CosMx / central','Xenium / central','Xenium / stratified','Visium HD / ovarian')
sp[,series:=factor(fcase(grepl('^cosmx',fixture),series[1],grepl('cervical_center',fixture),series[2],
  grepl('cervical_stratified',fixture),series[3],default=series[4]),levels=series)]
platform_cols<-setNames(c(blue,red,red,gold),series);platform_shapes<-setNames(c(16,17,15,18),series)
setorder(sp,N,fixture)
sp_base<-function(d,y,ylab,logy=TRUE) {
  p<-ggplot(d,aes(N,.data[[y]],colour=series,shape=series))+
    geom_line(data=d[grepl('^cosmx',fixture)],aes(group=series),linewidth=.65,show.legend=FALSE)+
    geom_point(size=2,stroke=.45)+scale_colour_manual(values=platform_cols)+scale_shape_manual(values=platform_shapes)+
    scale_x_log10(breaks=c(1000,10000,100000),labels=cell_lab,limits=range(sp$N)*c(.85,1.15))+
    labs(x='Cells (log scale)',y=ylab,colour=NULL,shape=NULL)+theme(legend.position='none',
      axis.title=element_text(size=7),axis.text=element_text(size=6.8),plot.margin=margin(2,3,1,2))
  p+if(logy)scale_y_log10(labels=label_log())else scale_y_continuous(limits=c(0,NA),breaks=c(0,10,20,30))
}
pairsp<-unique(sp[,.(N,N_squared)])
pairs<-sp_base(sp,'E','Pairs / directed edges')+
  geom_line(data=pairsp,aes(N,N_squared),inherit.aes=FALSE,linetype='dashed',linewidth=.35,colour=neutral)+
  annotate('text',x=max(sp$N),y=max(sp$N_squared)*1.3,label='N^2',parse=TRUE,hjust=1,size=2.6,colour=neutral)+
  scale_y_log10(breaks=c(1e2,1e5,1e8,1e11),labels=label_log())
degree<-sp_base(sp,'mean_out_degree','Mean neighbors, E/N',FALSE)
density<-sp_base(sp,'graph_density',expression('Graph density, '*E/N^2))+
  scale_y_log10(breaks=c(1e-4,1e-3,1e-2,1e-1),labels=label_log())
legend<-get_legend(sp_base(sp,'E','E')+theme(legend.position='top',legend.text=element_text(size=6.8),
  legend.margin=margin(0,0,0,0),legend.key.width=unit(3,'mm'))+guides(colour=guide_legend(nrow=1),shape=guide_legend(nrow=1)))
f1f<-ggdraw()+draw_plot(plot_grid(pairs,degree,density,nrow=1),0,.02,1,.73)+draw_plot(legend,.02,.74,.98,.15)+
  draw_label('F',.004,.99,hjust=0,vjust=1,size=10,fontfamily=font,fontface='bold')+
  draw_label('Local graph sparsity across platforms',.065,.98,hjust=0,vjust=1,size=8,fontfamily=font,fontface='bold')+
  draw_label('Effective radius: 35 \u00b5m',.99,.98,hjust=1,vjust=1,size=7,fontfamily=font,colour=neutral)
for(y in c('E','mean_out_degree','graph_density'))provenance('Figure1',paste0('F / ',y),spfile,
  paste('fixture;N;N_squared;validation_status;',y),
  'All 14 PASS fixtures; only nested CosMx points connected; E counts directed non-self edges',nrow(sp))
reps<-sp[match(c('cosmx_50000',sp[grepl('^cosmx',fixture)][which.max(N),fixture],
  'cervical_center_5000','cervical_stratified_5000','ovarian_16247'),fixture)]
reps[,label:=paste0(c('CosMx','CosMx','Xenium central','Xenium stratified','Visium HD'),'\n',count(N))]
reps[,label:=factor(label,levels=label)]
storage<-ggplot(reps,aes(label,colour=series))+
  geom_segment(aes(y=csr_total_bytes/1e6,yend=dense_float64_bytes/1e6,xend=label),linewidth=.75)+
  geom_point(aes(y=csr_total_bytes/1e6,fill=series),shape=21,size=2.6,stroke=.7)+
  geom_point(aes(y=dense_float64_bytes/1e6),shape=21,fill='white',size=2.6,stroke=.7)+
  scale_colour_manual(values=platform_cols,guide='none')+scale_fill_manual(values=platform_cols,guide='none')+
  scale_y_log10(labels=label_log(),breaks=10^seq(-2,6,2))+
  scale_x_discrete(expand=expansion(add=.6))+labs(x=NULL,y='Storage (MB)')+
  theme(axis.text.x=element_text(size=7),plot.margin=margin(2,5,0,3))
f1g<-panel(storage,'G','Theoretical dense storage versus measured CSR arrays',
  'Open: dense float64 (8N\u00b2 bytes). Filled: offsets + neighbors + distances. Array payload, not peak RSS.')
provenance('Figure1','G',spfile,'fixture;N;csr_total_bytes;dense_float64_bytes;offset_bytes;neighbor_bytes;distance_bytes',
  'Five prespecified representative fixtures; decimal MB; core array payload excludes R overhead and extra arrays',nrow(reps))
strip<-copy(longsc);strip[,x:=match(cells,sc$cells)]
strip[,method:=factor(method,levels=c('SpatialESS','SpatialCellChat V3'))]
strip[,state:=ifelse(status!='completed','Failed',ifelse(cells>max(ac$cells),'ESS-only','Paired'))]
fstrip<-ggplot(strip,aes(x,method))+
  annotate('rect',xmin=sum(sc$cells<=max(ac$cells))+.5,xmax=nrow(sc)+.5,ymin=.5,ymax=2.5,fill=pale,colour=NA)+
  geom_point(aes(shape=state,colour=method),size=1.8,stroke=.7)+scale_colour_manual(values=method_cols,guide='none')+
  scale_shape_manual(values=c('Failed'=4,'ESS-only'=1,'Paired'=16))+
  scale_x_continuous(breaks=seq_len(nrow(sc)),labels=cell_lab(sc$cells),position='top',expand=expansion(add=.4))+
  labs(x=NULL,y=NULL,shape=NULL)+theme_minimal(base_family=font,base_size=7)+theme(panel.grid=element_blank(),
    axis.text=element_text(colour=ink,size=6.8),legend.position='right',legend.text=element_text(size=6.5),
    legend.key.height=unit(2.5,'mm'),legend.margin=margin(0,0,0,0),plot.margin=margin(0,3,0,7))
f1h<-ggdraw()+draw_plot(fstrip,0,0,1,.73)+
  draw_label('H',.004,.99,hjust=0,vjust=1,size=10,fontfamily=font,fontface='bold')+
  draw_label('Completed comparisons and unpaired scaling',.065,.98,hjust=0,vjust=1,size=8,fontfamily=font,fontface='bold')
provenance('Figure1','H','source_data/cosmx_scaling.tsv.gz; source_data/cosmx_fidelity.tsv.gz',
  'cells;spatialess_status;official_status','Filled paired; open ESS-only; crosses failed; unpaired scales provide no fidelity evidence',nrow(strip))

# Dimensions in mm keep row heights explicit and make the layout reproducible.
fig1<-ggdraw()+draw_plot(f1b,0,138/180,.40,42/180)+draw_plot(f1c,.405,138/180,.595,42/180)+
  draw_plot(f1d,0,97/180,.50,41/180)+draw_plot(f1e,.50,97/180,.50,41/180)+
  draw_plot(f1f,0,55/180,1,42/180)+draw_plot(f1g,0,18/180,1,37/180)+draw_plot(f1h,0,0,1,18/180)
export(fig1,'Figure1_data',180,180)

# Retain the audited primary-cohort selections and feature identities exactly.
code<-paste(readLines(file.path(root,'source_snapshot/figure2_executed.R')),collapse='\n')
replace_once<-function(old,new) {
  if(!grepl(old,code,fixed=TRUE))stop('Missing archived code fragment: ',old)
  code<<-sub(old,new,code,fixed=TRUE)
}
replace_once("file.path(root,'audit/FIGURE2D_IDENTITIES.tsv')","file.path(root,'source_snapshot/previous_audit/FIGURE2D_IDENTITIES.tsv')")
replace_once("label:=paste(ligand,gsub('_','+',receptor,fixed=TRUE),sep='-')",
  "label:=ifelse(grepl('_',receptor,fixed=TRUE),paste0(ligand,'-(',gsub('_','+',receptor,fixed=TRUE),')'),paste(ligand,receptor,sep='-'))")
replace_once("forest[,label:=paste0(display_id,'  ',ligand,'-',gsub('_','+',receptor,fixed=TRUE))",
  "forest[,label:=paste0(display_id,'  ',ifelse(grepl('_',receptor,fixed=TRUE),paste0(ligand,'-(',gsub('_','+',receptor,fixed=TRUE),')'),paste(ligand,receptor,sep='-')))")
replace_once("lr_label:=paste0(display_id,'  ',ligand,'-',gsub('_','+',receptor,fixed=TRUE))",
  "lr_label:=paste0(display_id,'  ',ifelse(grepl('_',receptor,fixed=TRUE),paste0(ligand,'-(',gsub('_','+',receptor,fixed=TRUE),')'),paste(ligand,receptor,sep='-')))")
replace_once("key<-paste(strwrap(key,width=34),collapse='\\n')",
  "key<-paste0('ID: ',paste(strwrap(gsub(':',': ',z$feature_id,fixed=TRUE),width=37),collapse='\\n'))")
replace_once("size=5.15,fontfamily=font","size=6,fontfamily=font")
replace_once("y=.75,vjust=1","y=.74,vjust=1")
replace_once("draw_plot(q,0,0,1,.60)","draw_plot(q,0,0,1,.58)")
replace_once("size=6.15,fontfamily=font","size=6.5,fontfamily=font")
replace_once("Support = nonzero slices / 45 (0-1); ICC spans 0-1. Both models use the same cohort.",
  "Nonzero-slice fraction = nonzero / 45; bars span 0-1. LMM and OLS: same cohort.")
code<-gsub("size=2.3,family=font","size=2.1,family=font",code,fixed=TRUE)
code<-gsub("size=.42,alpha=.55","size=.5,alpha=.8",code,fixed=TRUE)
code<-gsub("'#BBC2C6'","'#A4ADB3'",code,fixed=TRUE)
code<-gsub("size=.7,alpha=.3","size=.85,alpha=.5",code,fixed=TRUE)
code<-gsub("colour=neutral)","colour=ink)",code,fixed=TRUE)
code<-gsub("linewidth=1.6,colour=ink","linewidth=1.8,colour=teal",code,fixed=TRUE)
code<-gsub("max_size=3.4","max_size=3.7",code,fixed=TRUE)
code<-gsub("axis.text.x=element_text(angle=0,hjust=.5,vjust=1,size=6.2)",
  "axis.text.x=element_text(angle=0,hjust=.5,vjust=1,size=6.2,colour=ink)",code,fixed=TRUE)
code<-gsub("axis.text.y=element_text(size=6.5)","axis.text.y=element_text(size=6.5,colour=ink)",code,fixed=TRUE)
code<-gsub("guide_colourbar(","guide_colourbar(display='rectangles',",code,fixed=TRUE)
replace_once("draw_plot(f2a,0,.78,.56,.22)+draw_plot(f2b,.56,.78,.44,.22)+\n  draw_plot(f2c,0,.405,1,.375)+draw_plot(f2d,0,.15,1,.255)+draw_plot(f2e,0,0,1,.15)",
  "draw_plot(f2a,0,189/238,.55,49/238)+draw_plot(f2b,.555,189/238,.445,49/238)+\n  draw_plot(f2c,0,85/238,1,104/238)+draw_plot(f2d,0,31/238,1,54/238)+draw_plot(f2e,0,0,1,31/238)")
replace_once("export(fig2,'Figure2',180,236)","export(fig2,'Figure2',180,238)")
writeLines(code,file.path(root,'scripts/figure2_rendered_source.R'))
eval(parse(text=code))
stopifnot(nrow(hm)==28,nrow(forest)==12,nrow(ex)==3)
for(z in list(c('A','lmm;volcano_labels','status;estimate;q_value;feature_id;ligand;receptor',nrow(valid)),
  c('B','forest','feature_id;estimate;lower;upper;singular',nrow(forest)),
  c('C','effect_heatmap;selected_interactions','feature_id;ligand;receptor;standardized_effect;glm_standardized_effect;q_value;glm_q;support_fraction;icc;singular',nrow(hm)),
  c('D','patient_examples;slice_scores;patient_scores','feature_id;sender_group_name;receiver_group_name;condition;log1p_score;estimate;lower;upper',nrow(ss)),
  c('E','pathway_pair_representatives','pathway_name;sender_group_name;receiver_group_name;standardized_effect;q_value',nrow(bubbles)))) {
  provenance('Figure2',z[1],paste0('source_data/',strsplit(z[2],';')[[1]],'.tsv.gz',collapse=';'),z[3],
    'Archived selections unchanged; BH across valid status==ok fits; singular fits retained and marked; effect colors centered at zero',as.integer(z[4]))
}

# Supplementary figures retain their data and panels; use the same visual grammar.
original<-readLines(file.path(root,'source_snapshot/redraw_figures_20260920.R'))
start<-grep('^# Supplement 1:',original)[1];end<-grep("^artwork\\[\\['Figure1A'",original)[1]-1L
code<-paste(original[start:end],collapse='\n')
code<-gsub('#EAF0F3','#8BB3D9',code,fixed=TRUE)
code<-gsub('#F8ECE8','#DE98A2',code,fixed=TRUE)
code<-gsub("'#4B9891','#A18545','#765E82','#CED2D5'","'#5C96C8','#D98A96','#515C64','#C1C8CD'",code,fixed=TRUE)
code<-gsub("colour=red","colour=blue",code,fixed=TRUE)
code<-gsub("size=.42,alpha=.8","size=.49,alpha=1",code,fixed=TRUE)
code<-gsub("size=.7,colour=blue,alpha=.65","size=.9,colour=blue,alpha=.9",code,fixed=TRUE)
code<-gsub("'#EAF1F5'","'#C4DDEA'",code,fixed=TRUE)
code<-gsub("'#F2F3F4'","'#D8E0E3'",code,fixed=TRUE)
replace_once("'E: old / new engine | B: bridge / direct'","'E: old / new engine; B: bridge / direct'")
replace_once("'Same corrected input, labels, LR and inference settings'","'Shared corrected input, cell order, LR and inference settings'")
replace_once("yend=.28,linewidth=.4,colour=neutral","yend=.28,linewidth=.4,colour=red")
replace_once("label='SpatialCellChat V3',hjust=0,family=font,size=2.6,colour=neutral",
  "label='SpatialCellChat V3',hjust=0,family=font,size=2.6,colour=red")
replace_once("'Independent cervical Xenium tissue'","'Independent cervical Xenium'")
replace_once("' stratified cells; public clusters, not curated cell types'","' cells; cell-type-stratified sampling'")
replace_once("'Paired fidelity across sampling designs'","'Fidelity across sampling designs'")
replace_once("'Default filters; matched records; no spatial-null interpretation'","'Default filters; matched records'")
replace_once("cr[,case_label:=gsub('_',' ',case)]","cr[,case_label:=paste(ifelse(grepl('center',case),'Central','Stratified'),cell_lab(as.numeric(sub('.*_','',case))))]")
replace_once("draw_plot(s1a,0,.78,1,.22)+draw_plot(s1b,0,.44,.51,.34)+draw_plot(s1c,.51,.44,.49,.34)+\n  draw_plot(s1d,0,.08,.47,.36)+draw_plot(s1e,.48,.26,.52,.18)+draw_plot(s1f,.48,0,.52,.26)",
  "draw_plot(s1a,0,.84,1,.16)+draw_plot(s1b,0,.43,.54,.41)+draw_plot(s1c,.55,.43,.45,.41)+\n  draw_plot(s1d,0,0,.47,.43)+draw_plot(s1e,.48,.25,.52,.18)+draw_plot(s1f,.48,0,.52,.25)")
replace_once("export(figs1,'FigureS1',180,184)","export(figs1,'FigureS1',180,170)")
replace_once("draw_plot(s2a,0,.52,.44,.48)+draw_plot(s2b,.45,.52,.55,.48)+\n  draw_plot(s2c,0,0,.44,.52)+draw_plot(s2d,.45,0,.55,.52)",
  "draw_plot(s2a,0,.54,.43,.46)+draw_plot(s2b,.44,.54,.56,.46)+\n  draw_plot(s2c,0,0,.43,.54)+draw_plot(s2d,.44,0,.56,.54)")
replace_once("export(figs2,'FigureS2',180,174)","export(figs2,'FigureS2',180,170)")
writeLines(code,file.path(root,'scripts/supplement_rendered_source.R'))
eval(parse(text=code))
for(z in list(c('A','bridge_regression','inputs_identical;full_records_identical',NA),
 c('B','bridge_regression','case;comparison;max_abs_probability_difference;pvalue_exact_fraction;active_jaccard;significance_disagreements',nrow(br)),
 c('C','sparkle_probability_records;sparkle_fidelity','official_probability;spatialess_probability;max_abs_probability_difference',nrow(spr)),
 c('D','sparkle_probability_records;sparkle_fidelity','official_pvalue;spatialess_pvalue;pvalue_exact_fraction',nrow(spr)),
 c('E','sparkle_probability_records;sparkle_fidelity','official_pvalue;spatialess_pvalue;active_jaccard;significance_disagreements',nrow(spr)),
 c('F','raw_corrected_transitions;raw_corrected_metrics','category;records;metric;value',nrow(trans)))) {
 provenance('FigureS1',z[1],paste0('source_data/',strsplit(z[2],';')[[1]],'.tsv.gz',collapse=';'),z[3],
  'Original corrected-input comparison and saved-input regressions; unadjusted p < 0.05 for calls; RAW/corrected is input sensitivity',as.integer(z[4]))
}
for(z in list(c('A','cervical_spatial_metadata','x_centroid;y_centroid;group',nrow(cm)),
 c('B','cervical_fidelity','case;filter;cells;max_abs_probability_difference;pvalue_exact_fraction;active_jaccard;significance_disagreements',nrow(cv)),
 c('C','cervical_probability_records','case;probability_official;probability_ess',nrow(cr)),
 c('D','cervical_fidelity','case;filter;ess_active;official_active;ess_significant_lt_005;official_significant_lt_005',nrow(cv)))) {
 provenance('FigureS2',z[1],paste0('source_data/',z[2],'.tsv.gz'),z[3],
  'Four frozen sampling fixtures; V3 dimensions-only repair; default and permissive filters; calls use unadjusted p < 0.05',as.integer(z[4]))
}
fwrite(rbindlist(prov),file.path(root,'audit/FIGURE_PROVENANCE.tsv'),sep='\t',na='NA')
print(rbindlist(prov))
fwrite(data.table(role=c('SpatialESS','SpatialCellChat V3','Positive effect / disease','Negative effect / control','CosMx','Xenium','Visium HD','Text'),
  hex=c(blue,red,red,blue,blue,red,gold,ink)),file.path(root,'audit/PALETTE.tsv'),sep='\t')
saveRDS(artwork,file.path(out,'editable_artwork.rds'))
writeLines(capture.output(sessionInfo()),file.path(root,'logs/plotting_sessionInfo.txt'))
message('All four figures generated from unchanged tables; former S3 is in Figure1 F-G.')
