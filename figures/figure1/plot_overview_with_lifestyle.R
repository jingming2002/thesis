#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(grid)})

args <- commandArgs(FALSE)
script <- sub('^--file=','',args[grepl('^--file=',args)])
if (sys.nframe() == 0 && length(script)) setwd(dirname(normalizePath(script[1])))
# Display-only revision: keep source data, colours and computations unchanged.
FONT_SCALE <- 1  # Increase to 1.1 for a further 10% enlargement.
fs <- function(pt) pt * FONT_SCALE
mm_text <- function(pt) fs(pt) / (72.27 / 25.4)
d <- read.delim('confidence_plot_data.tsv')
h <- read.delim('host_stacked_data.tsv')
cl <- read.delim('class_display.tsv')
fa <- read.delim('family_display.tsv')
stopifnot(sum(cl$vOTUs)==5384,sum(fa$vOTUs)==5384)
stopifnot(all(tapply(d$count,interaction(d$category,d$stage),sum)==tapply(d$total,interaction(d$category,d$stage),unique)))
d$category <- factor(d$category,c('AMG','APG','AReG'))
d$stage <- factor(d$stage,c('Candidate set','Filtered'))
d$confidence <- factor(d$confidence,c('high','medium','low'))
th <- theme_classic(base_size=fs(13.5),base_family='sans') +
 theme(axis.text=element_text(colour='#303030'),axis.ticks.x=element_blank(),
 axis.title.x=element_blank(),plot.title=element_text(size=fs(15),face='bold',margin=margin(b=9)),
 plot.title.position='plot',plot.margin=margin(9,12,9,10),legend.position='none',
 strip.background=element_blank(),strip.text=element_text(face='bold',size=fs(13.5)),
 panel.spacing=unit(10,'pt'))
tot <- unique(d[c('category','stage','total')])
conf <- ggplot(d,aes(stage,fill=confidence))+facet_grid(.~category)+th+
 scale_fill_manual(values=c(high='#64BC8C',medium='#DB9D85',low='#C5C9CD'),
 labels=c('High','Medium','Low'),name='CheckAMG confidence')+
 scale_x_discrete(labels=c('Candidate set'='Candidate\nset','Filtered'='Filtered'))
a <- conf+geom_col(aes(y=count),width=.64,position=position_stack(reverse=TRUE))+
 geom_text(data=tot,aes(stage,total,label=format(total,big.mark=',',trim=TRUE)),
 inherit.aes=FALSE,vjust=-.45,size=mm_text(12))+
 scale_y_continuous(labels=scales::label_comma(),breaks=seq(0,20000,5000),limits=c(0,21000),expand=c(0,0))+
 labs(title='A  AVG candidate counts',y='Candidate records')
b <- conf+geom_col(aes(y=percent),width=.64,position=position_stack(reverse=TRUE))+
 scale_y_continuous(breaks=seq(0,100,25),labels=function(x)paste0(x,'%'),expand=c(0,0))+
 labs(title='B  Confidence composition',y='Proportion')
extract_legend <- function(p,ncol) {
 g<-ggplotGrob(p+theme(legend.position='bottom',legend.title=element_text(size=fs(13),margin=margin(b=5)),
 legend.text=element_text(size=fs(12.5),margin=margin(l=3,r=12)),legend.key.size=unit(3.5,'mm'),legend.margin=margin(0,0,0,0))+
 guides(fill=guide_legend(ncol=ncol,byrow=TRUE,title.position='top',title.hjust=.5)))
 idx <- which(g$layout$name=='guide-box-bottom')
 if(!length(idx)) idx <- which(g$layout$name=='guide-box')
 g$grobs[[idx[1]]]
}
leg_conf <- extract_legend(b,3)
cols <- setNames(c('#9B85C5','#45A5AE','#E2AC4A','#D77B82','#639CC7','#7FA576','#B8B8B8')[seq_len(nrow(cl))],cl$taxon)
cols['No prediction'] <- '#DEDEDE';cols['Other classes'] <- '#A8B4BE'
fcols <- setNames(grDevices::hcl.colors(nrow(fa),'Set 3'),fa$taxon)
fcols['Family not reported'] <- '#B9BFC5';fcols['No prediction'] <- '#E5E5E5'
donut <- function(z,pal,title) {
 z$taxon<-factor(z$taxon,levels=z$taxon)
 lab<-paste0(gsub('_like','-like',z$taxon,fixed=TRUE),'\n',format(z$vOTUs,big.mark=',',trim=TRUE),
             ' (',sprintf('%.1f',100*z$vOTUs/5384),'%)')
 ggplot(z,aes(x=2,y=vOTUs,fill=taxon))+geom_col(width=.85,colour='white',linewidth=.22)+
 coord_polar(theta='y')+xlim(.5,2.5)+
 annotate('text',x=.5,y=0,label='5,384\nvOTUs',size=mm_text(13),lineheight=1.2)+
 scale_fill_manual(values=pal,labels=setNames(lab,z$taxon),name=NULL)+labs(title=title)+
 theme_void(base_size=fs(13.5),base_family='sans')+
 theme(plot.title=element_text(face='bold',size=fs(15),margin=margin(b=8)),
 plot.title.position='plot',plot.margin=margin(9,12,9,10),legend.position='right',legend.text=element_text(size=fs(11.8)),
 legend.key.size=unit(3.8,'mm'),legend.spacing.y=unit(3,'pt'),legend.margin=margin(0,0,0,1),
 legend.box.spacing=unit(3,'pt'))+guides(fill=guide_legend(ncol=1))
}
c <- donut(cl,cols,'C  Viral class')
dd <- donut(fa,fcols,'D  Viral family')
ph<-unique(h[order(h$display_order),c('host_phylum','display_order')])$host_phylum
h$host_phylum<-factor(h$host_phylum,levels=ph)
h$group<-factor(h$group,c('All AVG','AMG','APG','AReG'))
hcols<-c(Bacillota_A='#DB9D85',Bacteroidota='#64BC8C',Bacillota='#9CA9E2',
 Pseudomonadota='#E09A98',Bacillota_C='#80B87C',Actinomycetota='#80ADDC',Other='#C5C9CD')
host<-ggplot(h,aes(group,fill=host_phylum))+th+scale_fill_manual(values=hcols,name='Predicted host phylum',drop=FALSE)
htot<-unique(h[c('group','total_vOTUs')])
e<-host+geom_col(aes(y=vOTUs),width=.65,position=position_stack(reverse=TRUE))+
 geom_text(data=htot,aes(group,total_vOTUs,label=format(total_vOTUs,big.mark=',',trim=TRUE)),
 inherit.aes=FALSE,vjust=-.5,size=mm_text(12))+
 scale_y_continuous(labels=scales::label_comma(),breaks=c(0,1000,2000,3000),limits=c(0,3250),expand=c(0,0))+
 labs(title='E  Predicted host counts',y='Number of vOTUs')
f<-host+geom_col(aes(y=percent),width=.65,position=position_stack(reverse=TRUE))+
 scale_y_continuous(breaks=seq(0,100,25),labels=function(x)paste0(x,'%'),expand=c(0,0))+
 labs(title='F  Predicted host composition',y='Proportion')
leg_host<-extract_legend(f,4)
# Lifestyle counts: unique AVG-positive vOTUs in the eligible lifestyle set.
life_raw <- read.delim('vOTU_class_counts.tsv')
life_positive <- unique(life_raw[life_raw$Copies > 0,c('Contig','Lifestyle')])
stopifnot(!anyDuplicated(life_positive$Contig), nrow(life_positive)==4860)
life_counts <- as.data.frame(table(factor(life_positive$Lifestyle,levels=c('Virulent','Temperate'))))
names(life_counts) <- c('Lifestyle','vOTUs')
stopifnot(identical(as.integer(life_counts$vOTUs),c(2811L,2049L)))
write.table(life_counts,'lifestyle_counts.tsv',sep='\t',quote=FALSE,row.names=FALSE)
write.table(life_positive,'lifestyle_vOTUs.tsv',sep='\t',quote=FALSE,row.names=FALSE)
life_counts$Lifestyle <- factor(life_counts$Lifestyle,levels=c('Temperate','Virulent'))
g <- ggplot(life_counts,aes(x=Lifestyle,y=vOTUs,fill=Lifestyle))+
 geom_col(width=.55)+
 geom_text(aes(label=format(vOTUs,big.mark=',',trim=TRUE)),hjust=-.2,size=mm_text(12))+
 scale_fill_manual(values=c(Virulent='#CC6B2C',Temperate='#2C7F91'))+
 scale_y_continuous(limits=c(0,3200),breaks=c(0,1000,2000,3000),labels=scales::label_comma(),expand=c(0,0))+
 coord_flip()+th+
 labs(title='G  Predicted viral lifestyle',x=NULL,y='Number of AVG-positive vOTUs')+
 theme(axis.title.x=element_text(size=fs(13.5),margin=margin(t=6)),axis.title.y=element_blank())
ggsave('lifestyle_panel_R.pdf',g,width=12.8,height=2.3,device=cairo_pdf)
ggsave('lifestyle_panel_R.png',g,width=12.8,height=2.3,dpi=300,bg='white')
gs<-lapply(list(a,b,c,dd,e,f,g),ggplotGrob)
for(pair in list(c(1,2),c(5,6))) {
 w<-unit.pmax(gs[[pair[1]]]$widths,gs[[pair[2]]]$widths)
 for(j in pair)gs[[j]]$widths<-w
}
draw<-function(){
 grid.newpage()
 pushViewport(viewport(x=.5,y=.5,width=.985,height=.985,
 layout=grid.layout(6,2,widths=c(1,1),heights=unit(c(3.6,.7,4.1,3.6,1.0,2.3),'in'))))
 for(i in seq_len(6)) {
  row<-c(1,1,3,3,4,4)[i];col<-ifelse(i%%2==1,1,2)
  pushViewport(viewport(layout.pos.row=row,layout.pos.col=col));grid.draw(gs[[i]]);popViewport()
 }
 pushViewport(viewport(layout.pos.row=2,layout.pos.col=1:2));grid.draw(leg_conf);popViewport()
 pushViewport(viewport(layout.pos.row=5,layout.pos.col=1:2));grid.draw(leg_host);popViewport()
 pushViewport(viewport(layout.pos.row=6,layout.pos.col=1:2));grid.draw(gs[[7]]);popViewport(2)
}
cairo_pdf('AVG_overview_with_lifestyle_R.pdf',width=12.8,height=15.6);draw();dev.off()
png('AVG_overview_with_lifestyle_R.png',width=12.8,height=15.6,units='in',res=300,type='cairo',bg='white');draw();dev.off()
svglite::svglite('AVG_overview_with_lifestyle_R.svg',width=12.8,height=15.6);draw();dev.off()
writeLines(capture.output(sessionInfo()),'R_sessionInfo.txt')
