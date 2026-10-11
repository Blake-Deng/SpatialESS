# Figure 2: verified primary cohort, without implying four independent cohorts.
lmm <- rd('lmm'); meta <- rd('sample_metadata'); sel <- rd('selected_interactions')
hm <- rd('effect_heatmap'); forest <- rd('forest'); ex <- rd('patient_examples')
ss <- rd('slice_scores'); ps <- rd('patient_scores'); bubbles <- rd('pathway_pair_representatives')
condition_cols <- c(Control=blue,Disease=red)
valid <- lmm[status=='ok']; valid[,direction:=ifelse(q_value>=.05,'Not significant',ifelse(estimate>0,'PF higher','Control higher'))]
valid[,logq:=-log10(q_value)]
stopifnot(all(is.finite(valid$logq)))
lab <- head(rd('volcano_labels'),8); lab[,label:=ifelse(grepl('_',receptor,fixed=TRUE),paste0(ligand,'-(',gsub('_','+',receptor,fixed=TRUE),')'),paste(ligand,receptor,sep='-'))]
setorder(lab,estimate,feature_id)
# Fixed label slots retain the archived selection and render identically as
# editable PowerPoint text; ggrepel's deferred grobs were lost on that device.
lab[, label_x:=rep(c(-.00035,.00105,.00240),each=3)[seq_len(.N)]]
lab[, label_y:=rep(c(26,21,16),3)[seq_len(.N)]]
vp <- ggplot(valid,aes(estimate,logq))+
  geom_point(aes(colour=direction),size=.5,alpha=.8,stroke=0)+
  geom_hline(yintercept=-log10(.05),linewidth=.25,linetype='dashed',colour=ink)+
  geom_vline(xintercept=0,linewidth=.25,colour=ink)+
  geom_segment(data=lab,aes(x=estimate,y=-log10(q_value),xend=label_x,yend=label_y-.7),
    inherit.aes=FALSE,linewidth=.2,colour=ink)+
  geom_text(data=lab,aes(label_x,label_y,label=label),inherit.aes=FALSE,
    size=2.05,family=font,colour=ink)+
  scale_color_manual(values=c('PF higher'=red,'Control higher'=blue,'Not significant'='#A4ADB3'))+
  scale_x_continuous(labels=label_number(accuracy=.001),expand=expansion(mult=c(.07,.12)))+
  labs(x='PF - Control effect (log1p score)',y='-log10(BH FDR)',colour=NULL)+
  theme(legend.position='top',legend.key.width=unit(2.5,'mm'),legend.margin=margin(-4,0,-4,0))+
  guides(colour=guide_legend(override.aes=list(size=1.6,alpha=1)))
f2a <- panel(vp,'A','Primary LMM associations',paste0(count(nrow(valid)),' valid fits; ',count(sum(valid$q_value<.05)),' BH < 0.05'))

abbrev <- function(x) {
  x <- gsub('Proliferating ','Prolif. ',x,fixed=TRUE)
  x <- gsub('Macrophages - IFN-activated','IFN-activated Mph',x,fixed=TRUE)
  x <- gsub('Alveolar Macrophages','Alveolar Mph',x,fixed=TRUE)
  x <- gsub('Activated Fibrotic FBs','Act. fibrotic FBs',x,fixed=TRUE)
  x <- gsub('Interstitial Macrophages','Interst. Mph',x,fixed=TRUE)
  x
}
hm <- hm[match(sel$feature_id,feature_id)]
mat <- as.matrix(hm[,.(standardized_effect,glm_standardized_effect)])
colnames(mat) <- c('Patient\nLMM','Slice\nOLS')
rownames(mat) <- hm$feature_id
qmat <- cbind(hm$q_value,hm$glm_q)
rowlabels <- paste0(hm$display_id,'  ',hm$ligand,'-',hm$receptor,' | ',
  abbrev(hm$sender_group_name),' -> ',abbrev(hm$receiver_group_name))
lim <- max(abs(mat),na.rm=TRUE)
hmfun <- colorRamp2(c(-lim,0,lim),c(blue,'white',red))
hm[, row_label:=rowlabels]
stopifnot(all(hm$sender_group_name==hm$receiver_group_name))
hm[, lr_label:=paste0(display_id,'  ',ifelse(grepl('_',receptor,fixed=TRUE),paste0(ligand,'-(',gsub('_','+',receptor,fixed=TRUE),')'),paste(ligand,receptor,sep='-')))]
hm[, group_label:=abbrev(sender_group_name)]
hm[, yi:=-seq_len(.N)-cumsum(c(TRUE,diff(as.integer(factor(pathway_name,levels=unique(pathway_name)))))!=0)*.42]
hl <- rbind(data.table(yi=hm$yi,x=121,effect=mat[,1],q=qmat[,1]),
            data.table(yi=hm$yi,x=133,effect=mat[,2],q=qmat[,2]))
path_text <- hm[,.(yi=mean(yi)),by=pathway_name]
hp <- ggplot()+geom_tile(data=hl,aes(x,yi,fill=effect),width=9.5,height=.87)+
  effect_scale(lim)+
  geom_point(data=hl[is.finite(q)&q<.05],aes(x,yi),size=.45,colour=ink)+
  geom_text(data=hm,aes(x=21,y=yi,label=lr_label),hjust=0,size=2.1,family=font)+
  geom_text(data=hm,aes(x=75,y=yi,label=group_label),hjust=0,size=2.1,family=font)+
  geom_text(data=path_text,aes(x=0,y=yi,label=pathway_name),hjust=0,size=2.1,family=font,colour=ink)+
  geom_segment(data=hm,aes(x=143,xend=153,y=yi,yend=yi),linewidth=1.6,colour=pale)+
  geom_segment(data=hm,aes(x=143,xend=143+10*support_fraction,y=yi,yend=yi),linewidth=1.8,colour=teal)+
  geom_segment(data=hm,aes(x=160,xend=169,y=yi,yend=yi),linewidth=1.6,colour=pale)+
  geom_segment(data=hm,aes(x=160,xend=160+9*icc,y=yi,yend=yi),linewidth=1.6,colour=blue)+
  geom_point(data=hm[singular==TRUE],aes(x=178,y=yi),shape=15,size=1.45,colour=muted)+
  annotate('text',x=c(0,21,75),y=1.1,label=c('Pathway','Ligand-receptor','Cell group'),hjust=0,size=2.1,family=font)+
  annotate('text',x=75,y=-.15,label='(sender = receiver)',hjust=0,size=2.05,family=font,colour=ink)+
  annotate('text',x=c(121,133,148,164.5,178),y=.72,label=c('Patient\nLMM','Slice\nOLS','Nonzero-slice\nfraction','Patient\nICC','Singular'),
    size=2.25,family=font,lineheight=1.05)+
  coord_cartesian(xlim=c(0,184),ylim=c(min(hm$yi)-.6,2.9),expand=FALSE)+
  labs(fill='Effect / SD')+theme_void(base_family=font)+theme(legend.position='bottom',
    legend.title=element_text(size=7),legend.text=element_text(size=6.5),plot.margin=margin(0,3,0,3),
    legend.margin=margin(-1,0,0,0))+
  guides(fill=guide_colourbar(display='rectangles',direction='horizontal',barwidth=unit(28,'mm'),barheight=unit(2.5,'mm')))
f2c <- panel(hp,'C','Communication effects and model comparison',
  'Dot: BH q < 0.05; square: singular fit. Nonzero-slice fraction = nonzero / 45; bars span 0-1. LMM and OLS: same cohort.')

forest[,label:=paste0(display_id,'  ',ifelse(grepl('_',receptor,fixed=TRUE),paste0(ligand,'-(',gsub('_','+',receptor,fixed=TRUE),')'),paste(ligand,receptor,sep='-')))]
forest[,label:=factor(label,levels=rev(label))]
fp <- ggplot(forest,aes(estimate,label,colour=estimate>0))+
  geom_vline(xintercept=0,colour=neutral,linewidth=.3)+
  geom_segment(aes(x=lower,xend=upper,yend=label),linewidth=.5)+
  geom_point(aes(shape=singular),size=1.9,stroke=.6)+
  scale_colour_manual(values=c('TRUE'=red,'FALSE'=blue),guide='none')+
  scale_shape_manual(values=c('TRUE'=1,'FALSE'=16),labels=c('FALSE'='Non-singular','TRUE'='Singular'),name=NULL)+
  scale_x_continuous(labels=label_number(accuracy=.001))+
  labs(x='PF - Control effect (95% Wald CI)',y=NULL)+
  theme(legend.position='bottom',axis.text.y=element_text(size=6.5,colour=ink),legend.margin=margin(-3,0,-3,0))
f2b <- panel(fp,'B','Primary-cohort effect sizes')

bubbles[,pair_short:=paste(abbrev(sender_group_name),abbrev(receiver_group_name),sep='\n> ')]
bubbles[,logq:=pmin(-log10(q_value),10)]
bp <- ggplot(bubbles,aes(pair_short,pathway_name))+
  geom_point(aes(size=logq,fill=standardized_effect),shape=21,colour='white',stroke=.25)+
  effect_scale(max(abs(bubbles$standardized_effect)))+
  scale_size_area(max_size=3.7,breaks=c(2,5,10),limits=c(0,10))+
  labs(x=NULL,y=NULL,fill='Effect / SD',size='-log10(FDR)')+
  theme_minimal(base_family=font,base_size=7)+theme(panel.grid=element_blank(),
    axis.text.x=element_text(angle=0,hjust=.5,vjust=1,size=6.2,colour=ink),axis.text.y=element_text(size=6.5,colour=ink),
    legend.position='right',legend.key.size=unit(2,'mm'),legend.text=element_text(size=6),
    legend.title=element_text(size=6),legend.margin=margin(-3,0,-3,0),plot.margin=margin(3,3,1,3))+
  guides(fill=guide_colourbar(display='rectangles',direction='vertical',barwidth=unit(2,'mm'),barheight=unit(12,'mm')),
    size=guide_legend(override.aes=list(fill=neutral)))
f2e <- panel(bp,'E','Pathway / cell-pair representatives','Strongest qualifying LR per dot; FDR is feature-level')

identities<-fread(file.path(root,'source_snapshot/previous_audit/FIGURE2D_IDENTITIES.tsv'))
stopifnot(setequal(ex$feature_id,identities$feature_id))
ex<-merge(ex,identities[,.(feature_id,display_id)],by='feature_id',sort=FALSE)
setorder(ex,ligand,receptor,feature_id)
ss[,condition:=factor(condition,levels=c('Control','Disease'))]
ps[,condition:=factor(condition,levels=c('Control','Disease'))]
dplots<-lapply(seq_len(nrow(ex)),function(i){
 z<-ex[i];id<-z$feature_id
 s<-ss[feature_id==id];p<-ps[feature_id==id]
 stopifnot(nrow(s)==z$total_units,nrow(p)==z$patient_count)
 group_line<-paste0(z$sender_group_name,' -> ',z$receiver_group_name,' | ',z$display_id)
 group_line<-paste(strwrap(group_line,width=30),collapse='\n')
 key<-paste0('ID: ',z$feature_id)
 key<-paste0('ID: ',paste(strwrap(gsub(':',': ',z$feature_id,fixed=TRUE),width=37),collapse='\n'))
 q<-ggplot(s,aes(condition,log1p_score,colour=condition))+
   geom_point(position=position_jitter(width=.15,seed=13),size=.85,alpha=.5,stroke=0)+
   geom_point(data=p,position=position_jitter(width=.10,seed=19),size=1.45,alpha=.9,shape=21,fill='white',stroke=.45)+
   annotate('text',x=1,y=Inf,label=paste0('LMM effect ',num(z$estimate),'\n95% CI ',num(z$lower),' to ',num(z$upper)),hjust=0,vjust=1.1,size=2.05,family=font,colour=ink)+
   scale_colour_manual(values=condition_cols,guide='none')+scale_x_discrete(labels=c(Control='Control',Disease='PF'))+
   scale_y_continuous(expand=expansion(mult=c(.04,.40)))+
   labs(x=NULL,y=if(i==1)'log1p(score)'else NULL)+
   theme(axis.text=element_text(size=6.5),axis.title=element_text(size=7))
 ggdraw()+draw_plot(q,0,0,1,.58)+
   draw_label(paste(z$ligand,z$receptor,sep='-'),x=.5,y=.985,vjust=1,size=7.2,fontfamily=font,fontface='bold')+
   draw_label(group_line,x=.5,y=.91,vjust=1,size=6.5,fontfamily=font,lineheight=1.05)+
   draw_label(key,x=.5,y=.74,vjust=1,size=6,fontfamily=font,lineheight=1.05,colour=ink)
})
ep<-plot_grid(plotlist=dplots,nrow=1)
f2d<-panel(ep,'D','Slices and patient summaries','Light points: slices; outlined points: patient means. Intervals show the LMM condition effect.')
fig2 <- ggdraw()+draw_plot(f2a,0,189/238,.55,49/238)+draw_plot(f2b,.555,189/238,.445,49/238)+
  draw_plot(f2c,0,85/238,1,104/238)+draw_plot(f2d,0,31/238,1,54/238)+draw_plot(f2e,0,0,1,31/238)
export(fig2,'Figure2',180,238)
for (nm in c('a','b','c','d','e')) {
  p<-get(paste0('f2',nm)); ggsave(file.path(out,'panels',paste0('Figure2',toupper(nm),'.pdf')),p,
    width=if(nm%in%c('c','d','e'))180 else 90,height=if(nm=='c')92 else 65,units='mm',device=cairo_pdf)
}

