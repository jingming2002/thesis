args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else '.'
library(ggplot2)
library(grid)
d <- read.delim(file.path(out,'inputs/functional_mapping.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(d$Protein),nrow(d)==11725)
classes <- c('AMG','APG','AReG')
count <- function(x) length(unique(x))
m <- d[!is.na(d$Display_module)&nzchar(d$Display_module),]
pc <- aggregate(Protein~AVG_category+Display_module,m,count)
vc <- aggregate(Contig~AVG_category+Display_module,m,count)
names(pc)[3] <- 'Proteins';names(vc)[3] <- 'vOTUs'
tab <- merge(pc,vc)
names(tab)[1:2] <- c('Class','Module')
tab <- tab[order(match(tab$Class,classes),-tab$Proteins,tab$Module),]
tab$Label <- tab$Module
tab$Label[tab$Module=='Anaerobic RNR activation'] <- 'Anaerobic RNR activation (NrdG)'
tab$Label[tab$Module=='Cobalamin biosynthesis'] <- 'Cobalamin biosynthesis-related'
tab$Label[tab$Module=='Toxin-antitoxin-associated functions'] <- 'Toxin-antitoxin-associated candidates'
write.table(tab,file.path(out,'functional_module_protein_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
summary <- do.call(rbind,lapply(classes,function(cl){
 x <- d[d$AVG_category==cl,]; mapped <- sum(!is.na(x$Display_module)&nzchar(x$Display_module))
 data.frame(Class=cl,Total_proteins=nrow(x),Proteins_in_displayed_modules=mapped,Other_or_unresolved=nrow(x)-mapped)
}))
write.table(summary,file.path(out,'class_count_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
wrap <- function(x) vapply(x,function(s)paste(strwrap(s,width=30),collapse='\n'),character(1))
plots <- lapply(classes,function(cl){
 x <- tab[tab$Class==cl,]
 x$Row <- factor(x$Module,levels=rev(x$Module))
 lim <- c(AMG=760,APG=640,AReG=1490)[cl]
 br <- list(AMG=c(0,200,400,600),APG=c(0,200,400,600),AReG=c(0,400,800,1200))[[cl]]
 ggplot(x,aes(x=Proteins,y=Row))+
  geom_col(width=.59,fill='#377F92')+
  geom_text(aes(label=format(Proteins,big.mark=',',trim=TRUE)),hjust=-.18,size=3.3,colour='#24343E')+
  scale_y_discrete(labels=setNames(wrap(x$Label),x$Module),expand=expansion(add=.6))+
  scale_x_continuous(limits=c(0,lim),breaks=br,expand=expansion(mult=c(0,0)),labels=function(v)format(v,big.mark=',',trim=TRUE))+
  labs(title=cl,x='Number of proteins',y=NULL)+
  theme_minimal(base_size=11)+
  theme(plot.title=element_text(face='bold',size=15,colour='#22313D',margin=margin(b=16)),
   panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),panel.grid.major.x=element_line(colour='#E5EAED',linewidth=.35),
   axis.text.y=element_text(colour='#263640',size=10,lineheight=.95,margin=margin(r=8)),
   axis.text.x=element_text(colour='#61717D',size=9),axis.title.x=element_text(size=10,margin=margin(t=9)),
   plot.margin=margin(8,16,10,8))
})
# Keep the bar area equal in width across the three panels.
grobs <- lapply(plots,ggplotGrob)
w <- do.call(unit.pmax,lapply(grobs,function(g)g$widths))
for(i in seq_along(grobs)) grobs[[i]]$widths <- w
draw <- function(){
 grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
 grid.text('Functional composition of auxiliary viral genes',x=.017,y=.965,just='left',gp=gpar(fontsize=18,fontface='bold',col='#22313D'))
 for(i in seq_along(grobs)){
  pushViewport(viewport(x=(i-.5)/3,y=.465,width=1/3,height=.88))
  grid.draw(grobs[[i]]);popViewport()
 }
}
pdf(file.path(out,'AVG_functional_protein_counts.pdf'),width=17,height=8.2,useDingbats=FALSE);draw();dev.off()
png(file.path(out,'AVG_functional_protein_counts.png'),width=17,height=8.2,units='in',res=300,type='cairo');draw();dev.off()
png(file.path(out,'preview.png'),width=2040,height=984,res=120,type='cairo');draw();dev.off()
if(requireNamespace('svglite',quietly=TRUE)){svglite::svglite(file.path(out,'AVG_functional_protein_counts.svg'),width=17,height=8.2);draw();dev.off()}
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(summary)
