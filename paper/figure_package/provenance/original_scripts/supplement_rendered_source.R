# Supplement 1: corrected-input fidelity, bridge regression, input sensitivity.
br <- rd('bridge_regression'); spr <- rd('sparkle_probability_records'); spm <- rd('sparkle_fidelity')
trans <- rd('raw_corrected_transitions'); rawmetrics <- rd('raw_corrected_metrics')
swork <- ggplot()+xlim(0,1)+ylim(0,1)+theme_void()+
  annotate('rect',xmin=.015,xmax=.24,ymin=.27,ymax=.75,fill='#C4DDEA',colour=NA)+
  annotate('text',x=.127,y=.51,label='SPARKLE\ncorrected H5AD',family=font,size=2.6,lineheight=1.35)+
  annotate('segment',x=.25,xend=.34,y=.51,yend=.51,linewidth=.4,arrow=arrow(length=unit(1.3,'mm')))+
  annotate('rect',xmin=.35,xmax=.60,ymin=.27,ymax=.75,fill='#D8E0E3',colour=NA)+
  annotate('text',x=.475,y=.51,label='Sparse bridge\nR matrix + metadata',family=font,size=2.6,lineheight=1.35)+
  annotate('segment',x=.61,xend=.69,y=.51,yend=.71,linewidth=.4,colour=blue,arrow=arrow(length=unit(1.3,'mm')))+
  annotate('segment',x=.61,xend=.69,y=.51,yend=.28,linewidth=.4,colour=red,arrow=arrow(length=unit(1.3,'mm')))+
  annotate('text',x=.71,y=.71,label='SpatialESS',hjust=0,family=font,size=2.8,fontface='bold',colour=blue)+
  annotate('text',x=.71,y=.28,label='SpatialCellChat V3',hjust=0,family=font,size=2.6,colour=red)+
  annotate('text',x=.015,y=.045,label='Shared corrected input, cell order, LR and inference settings',hjust=0,family=font,size=2.4,colour=neutral)
s1a <- panel(swork,'A','Corrected-input comparison')
br[,row:=paste(gsub('sparkle','SPK',gsub('_',' ',case)),ifelse(comparison=='old_vs_new','E','B'),sep=' / ')]
cols <- c('Max |prob.\ndifference|','Exact p\nfraction','Active\nJaccard','Significance\ndisagreements')
bmat <- rbindlist(lapply(seq_along(cols),function(i) {
  v <- c('max_abs_probability_difference','pvalue_exact_fraction','active_jaccard','significance_disagreements')[i]
  data.table(row=br$row,col=cols[i],label=num(br[[v]]),tint='#8BB3D9')
}))
s1b <- panel(matrix_plot(bmat,br$row,cols,base=6.2),'B','Bridge and historical regression',
  'E: old / new engine; B: bridge / direct')
s1c <- panel(density_agreement(spr,'official_probability','spatialess_probability',
  'V3, log10(probability)','SpatialESS, log10(probability)',
  paste0('n = ',count(nrow(spr)),'\nMax error ',num(spm$max_abs_probability_difference))),
  'C','Corrected-input probability fidelity')
pcounts <- spr[,.(records=.N),by=.(official_pvalue,spatialess_pvalue)]
pp <- ggplot(pcounts,aes(official_pvalue,spatialess_pvalue))+
  geom_abline(slope=1,intercept=0,linewidth=.3,linetype='dashed',colour=neutral)+
  geom_point(aes(size=records),colour=blue,alpha=.75)+scale_size_area(max_size=4.5,name='Records')+
  coord_fixed(xlim=c(-.025,1.025),ylim=c(-.025,1.025))+
  annotate('text',x=.03,y=.98,label=paste0('Exact p fraction = ',num(spm$pvalue_exact_fraction)),
    hjust=0,vjust=1,family=font,size=2.6)+
  scale_x_continuous(breaks=c(0,.5,1))+scale_y_continuous(breaks=c(0,.5,1))+
  labs(x='SpatialCellChat V3 p-value',y='SpatialESS p-value')+
  theme(legend.position='right')
s1d <- panel(pp,'D','Permutation p-value fidelity')
sig_off <- spr$official_pvalue<.05; sig_ess <- spr$spatialess_pvalue<.05
sj <- sum(sig_off & sig_ess)/sum(sig_off | sig_ess)
calls <- data.table(row=c('Active records','p-values','Significance calls'),
  col='Agreement',label=c(paste0('Jaccard = ',num(spm$active_jaccard)),
    paste0('Exact fraction = ',num(spm$pvalue_exact_fraction)),paste0('Jaccard = ',num(sj))),tint='#8BB3D9')
s1e <- panel(matrix_plot(calls,calls$row,'Agreement',base=7),'E','Interaction-call agreement',
  paste0('Significance disagreements: ',num(spm$significance_disagreements),' | same corrected input'))
tp <- ggplot(trans,aes(records,factor(category,levels=rev(category))))+
  geom_col(aes(fill=category),width=.58)+geom_text(aes(label=count(records)),hjust=-.1,size=2.8,family=font)+
  scale_fill_manual(values=setNames(c(neutral,blue,red),trans$category),guide='none')+
  scale_x_continuous(expand=expansion(mult=c(0,.23)),labels=label_comma())+
  labs(x='Significant interactions',y=NULL)+theme(axis.text.y=element_text(size=7))
s1f <- panel(tp,'F','RAW / SPARKLE input sensitivity',
  paste0('Significant Jaccard = ',formatC(rawmetrics[metric=='significant_jaccard',value],format='f',digits=3),' | input sensitivity'))
figs1 <- ggdraw()+draw_plot(s1a,0,.84,1,.16)+draw_plot(s1b,0,.43,.54,.41)+draw_plot(s1c,.55,.43,.45,.41)+
  draw_plot(s1d,0,0,.47,.43)+draw_plot(s1e,.48,.25,.52,.18)+draw_plot(s1f,.48,0,.52,.25)
export(figs1,'FigureS1',180,170)

# Supplement 2: independent public cervical tissue and exact paired tests.
cm <- rd('cervical_spatial_metadata'); cv <- rd('cervical_fidelity'); cr <- rd('cervical_probability_records')
tops <- cm[, .N,by=group][order(-N,group),head(group,5)]
cm[,cluster_display:=ifelse(group %in% tops,group,'Other clusters')]
mapcols <- setNames(c(blue,red,'#5C96C8','#D98A96','#515C64','#C1C8CD'),c(tops,'Other clusters'))
mp <- ggplot(cm,aes(x_centroid,y_centroid))+
  geom_point(aes(colour=factor(cluster_display,levels=c(tops,'Other clusters'))),size=.49,alpha=1,stroke=0)+
  scale_colour_manual(values=mapcols,name='Public clusters')+coord_fixed()+
  labs(x='x (um)',y='y (um)')+theme(legend.position='bottom',legend.text=element_text(size=6),
    legend.title=element_text(size=6.5),legend.margin=margin(0,0,0,0))+
  guides(colour=guide_legend(ncol=3,override.aes=list(size=1.7,alpha=1)))
s2a <- panel(mp,'A','Independent cervical Xenium',paste0(count(nrow(cm)),' cells; cell-type-stratified sampling'))
cv[,row:=paste(ifelse(grepl('center',case),'Central','Stratified'),cell_lab(cells),
  ifelse(filter=='default','0.1 / 5','0 / 1'))]
cvals <- c('max_abs_probability_difference','pvalue_exact_fraction','active_jaccard','significance_disagreements')
cmt <- rbindlist(lapply(seq_along(cvals),function(i) data.table(row=cv$row,col=cols[i],
  label=num(cv[[cvals[i]]]),tint=if(i==1)'#DE98A2' else '#8BB3D9')))
s2b <- panel(matrix_plot(cmt,cv$row,cols,base=6.4),'B','Fidelity across sampling designs',
  'V3: dims-only repair | row suffix: min.percent / min.cells.sr')
cr <- cr[probability_official>0 & probability_ess>0]
cr[,case_label:=paste(ifelse(grepl('center',case),'Central','Stratified'),cell_lab(as.numeric(sub('.*_','',case))))]
scatter <- ggplot(cr,aes(log10(probability_official),log10(probability_ess)))+
  geom_abline(slope=1,intercept=0,linewidth=.25,linetype='dashed',colour=neutral)+
  geom_point(size=.9,colour=blue,alpha=.9,stroke=0)+facet_wrap(~case_label,nrow=2)+
  labs(x='V3, log10(probability)',y='SpatialESS, log10(probability)')+
  theme(axis.text=element_text(size=6.2),panel.spacing=unit(3,'mm'))
s2c <- panel(scatter,'C','Independent probability agreement','Default filters; matched records')
ct <- rbindlist(lapply(seq_len(nrow(cv)),function(i) {
  z<-cv[i]; rbindlist(lapply(c('Active','Significant'),function(st) data.table(
    sampling=ifelse(grepl('center',z$case),'Central','Stratified'),cells=cell_lab(z$cells),
    filter=ifelse(z$filter=='default','0.1 / 5','0 / 1'),measure=st,
    method=c('SpatialESS','SpatialCellChat V3'),records=if(st=='Active')c(z$ess_active,z$official_active) else
      c(z$ess_significant_lt_005,z$official_significant_lt_005))))
}))
ct[,case_label:=paste(sampling,cells)]
ct[,measure:=factor(measure,levels=c('Active','Significant'))]
counts <- ggplot(ct,aes(filter,records,fill=method))+
  geom_col(position=position_dodge(width=.75),width=.65)+
  facet_grid(measure~case_label,scales='free_y')+scale_fill_manual(values=method_cols)+
  scale_y_continuous(labels=label_comma(),expand=expansion(mult=c(0,.1)))+
  labs(x='min.percent / min.cells.sr',y='Records',fill=NULL)+theme(legend.position='bottom',
    axis.text.x=element_text(size=5.8,angle=45,hjust=1),strip.text=element_text(size=6.1),
    legend.margin=margin(-3,0,-3,0),panel.spacing=unit(1.5,'mm'))
s2d <- panel(counts,'D','Active and significant records','Significant: permutation p < 0.05, not BH-adjusted')
figs2 <- ggdraw()+draw_plot(s2a,0,.54,.43,.46)+draw_plot(s2b,.44,.54,.56,.46)+
  draw_plot(s2c,0,0,.43,.54)+draw_plot(s2d,.44,0,.56,.54)
export(figs2,'FigureS2',180,170)
