args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else '.'
library(grid)
read_tsv <- function(f) read.delim(file.path(out,f),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA',''))
write_tsv <- function(x,f) write.table(x,file.path(out,f),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
d <- read_tsv('inputs/original_mapping.tsv')
b <- read_tsv('inputs/eligible_background.tsv')
stopifnot(!anyDuplicated(b$Contig),!anyDuplicated(d$Protein))
d$Revision <- 'Original annotation-based module retained'
exclude_ta <- c('toxin','Bacterial toxin YdaS','Bacterial toxin ydaT',
 'hlg, luk; leukocidin/hemolysin toxin family protein',
 'Bacterial toxin homologue of phage lysozyme, C-term',
 'Papain fold toxin 1, glutamine deamidase',
 'berB; beta-exotoxin I transport system permease protein','ttsA; typhoid toxin secretion A')
bad <- which(d$Display_module=='Toxin-antitoxin-associated functions' & d$Function %in% exclude_ta)
d$Display_module[bad] <- NA
d$Revision[bad] <- 'Omitted from TA module: nonspecific toxin call or no explicit TA-family support'
tr <- which(d$Display_module=='Transcriptional repression' & d$Function=='translation repressor')
d$Display_module[tr] <- NA
d$Revision[tr] <- 'Translation repression omitted from transcriptional repression'
write_tsv(d,'revised_mapping_audit.tsv')
write_tsv(d[c(bad,tr),],'mapping_changes.tsv')
z <- merge(d[!is.na(d$Display_module),c('Contig','AVG_category','Display_module')],b[,c('Contig','PhaTYP_Pred')],by='Contig')
z <- unique(z)
write_tsv(z,'unique_module_carriers.tsv')
nv <- sum(b$PhaTYP_Pred=='virulent'); nt <- sum(b$PhaTYP_Pred=='temperate')
stopifnot(nv==6061,nt==3076)
mods <- unique(d[!is.na(d$Display_module),c('AVG_category','Display_module')])
res <- do.call(rbind,lapply(seq_len(nrow(mods)),function(i){
 x <- z[z$AVG_category==mods$AVG_category[i] & z$Display_module==mods$Display_module[i],]
 kv <- sum(x$PhaTYP_Pred=='virulent'); kt <- sum(x$PhaTYP_Pred=='temperate')
 test <- fisher.test(matrix(c(kv,nv-kv,kt,nt-kt),nrow=2,byrow=TRUE))
 data.frame(Class=mods$AVG_category[i],Module=mods$Display_module[i],Virulent_k=kv,Virulent_n=nv,
 Temperate_k=kt,Temperate_n=nt,Pooled_carriers=kv+kt,V_pct=100*kv/nv,T_pct=100*kt/nt,
 OR=unname(test$estimate),CI_lower=test$conf.int[1],CI_upper=test$conf.int[2],p=test$p.value)
}))
# BH family and display both contain all 32 mapped modules.
res$q <- p.adjust(res$p,'BH')
res$Displayed <- TRUE
classes <- c('AMG','APG','AReG')
stopifnot(nrow(res)==32)
res <- res[order(match(res$Class,classes),-res$Pooled_carriers,res$Module),]
write_tsv(res,'all_module_tests.tsv')
p <- res[res$Displayed,]
labels <- c('Anaerobic RNR activation'='Anaerobic RNR activation (NrdG)',
 'Folate-mediated one-carbon metabolism'='Folate-mediated one-carbon metabolism',
 'Cobalamin biosynthesis'='Cobalamin biosynthesis-related',
 'Sulfur mobilization / Fe-S assembly'='Sulfur mobilization / Fe-S assembly',
 'Toxin-antitoxin-associated functions'='Toxin-antitoxin-associated candidates',
 'Sporulation-associated functions'='Sporulation-associated functions',
 'Transcription initiation / sigma factors'='Transcription initiation / sigma factors',
 'Transcription regulation (unspecified)'='Transcription regulation (unspecified)')
p$Label <- p$Module
ii <- match(p$Module,names(labels)); p$Label[!is.na(ii)] <- labels[ii[!is.na(ii)]]
write_tsv(p,'plotted_data.tsv')
pal <- c(Virulent='#CC6B2C',Temperate='#2C7F91')
ink <- '#22313D'; muted <- '#64727D'; light <- '#E4E9ED'
# Fixed coordinates in a vector grid layout; all panels share both axis scales.
xl <- .35; xr <- .655; fl <- .735; fr <- .937; qx <- .978
xp <- function(v) xl+(xr-xl)*v/25
xf <- function(v) fl+(fr-fl)*(log2(v)+5)/9
txt <- function(label,x,y,size=10,col=ink,just='left',bold=FALSE){grid.text(label,x,y,just=just,gp=gpar(fontsize=size,col=col,fontface=if(bold)'bold' else 'plain',fontfamily='sans'))}
seg <- function(x0,y0,x1,y1,col=light,lwd=.7,lty=1){grid.segments(x0,y0,x1,y1,gp=gpar(col=col,lwd=lwd,lty=lty))}
pt <- function(x,y,fill,col=fill,r=1.5){grid.circle(x,y,r=unit(r,'mm'),gp=gpar(fill=fill,col=col,lwd=1))}
fmtq <- function(v) if(v<.001) '<0.001' else sprintf('%.3f',v)
draw <- function(){
 grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
 txt('Functional repertoires across viral lifestyles',.025,.972,17,bold=TRUE)
 pt(.028,.931,pal['Virulent']);txt(sprintf('Virulent  (n = %s)',format(nv,big.mark=',')),.040,.931,10)
 pt(.225,.931,pal['Temperate']);txt(sprintf('Temperate  (n = %s)',format(nt,big.mark=',')),.237,.931,10)
 txt('Function',.025,.887,10,bold=TRUE)
 txt('vOTUs carrying the function (%)',(xl+xr)/2,.887,10,just='centre',bold=TRUE)
 txt('Odds ratio (95% CI)',(fl+fr)/2,.897,10,just='centre',bold=TRUE)
 txt('Temperate-associated',fl,.876,8,col=pal['Temperate'])
 txt('Virulent-associated',fr,.876,8,col=pal['Virulent'],just='right')
 txt('BH q',qx,.887,9,just='right',bold=TRUE)
 rowstep <- .019
 counts <- as.numeric(table(factor(p$Class,levels=classes)))
 tops <- .837-c(0,cumsum(head(counts,-1)*rowstep+.063))
 titles <- c('AMG','APG','AReG')
 for(k in seq_along(classes)){
  dat <- p[p$Class==classes[k],]; top <- tops[k]
  grid.rect(x=.5,y=top,width=.95,height=.023,gp=gpar(fill='#F0F3F5',col=NA))
  txt(titles[k],.032,top,11,bold=TRUE)
  ys <- top-.030-(seq_len(nrow(dat))-1)*rowstep
  bottom <- min(ys)-.015; upper <- max(ys)+.014
  for(v in seq(0,25,5)) seg(xp(v),bottom,xp(v),upper)
  for(v in c(1/16,1/4,1,4,16)) seg(xf(v),bottom,xf(v),upper,col=if(v==1)'#7D8992' else light,lty=if(v==1)2 else 1)
  for(j in seq_len(nrow(dat))){
   a <- dat[j,];y <- ys[j]
   txt(a$Label,.032,y,9.5)
   seg(xp(a$T_pct),y,xp(a$V_pct),y,col='#BBC5CC',lwd=1.1)
   # Slight vertical offsets keep equal or nearly equal frequencies legible.
   pt(xp(a$V_pct),y+.0027,pal['Virulent'],r=1.35)
   pt(xp(a$T_pct),y-.0027,pal['Temperate'],r=1.35)
   colour <- if(a$OR>1)pal['Virulent'] else pal['Temperate']
   lo <- max(a$CI_lower,2^-5);hi <- min(a$CI_upper,16)
   seg(xf(lo),y,xf(hi),y,col=colour,lwd=1.3)
   seg(xf(lo),y-.0035,xf(lo),y+.0035,col=colour,lwd=1.1)
   seg(xf(hi),y-.0035,xf(hi),y+.0035,col=colour,lwd=1.1)
   if(is.finite(a$OR) && a$OR>=2^-5 && a$OR<=16){
    pt(xf(a$OR),y,if(a$q<.05)colour else 'white',col=colour,r=1.5)
   } else {
    txt(if(is.infinite(a$OR)) 'Inf' else format(a$OR,digits=2),fr-.008,y+.009,8,col=colour,just='right')
   }
   if(a$CI_upper>16) grid.segments(fr-.012,y,fr,y,arrow=arrow(length=unit(1.6,'mm')),gp=gpar(col=colour,lwd=1.3))
   if(a$CI_lower<2^-5) grid.segments(fl+.012,y,fl,y,arrow=arrow(length=unit(1.6,'mm')),gp=gpar(col=colour,lwd=1.3))
   txt(fmtq(a$q),qx,y,8.5,col=if(a$q<.05)ink else muted,just='right')
  }
  ay <- bottom-.012
  for(v in seq(0,25,5)) txt(as.character(v),xp(v),ay,8.5,col=muted,just='centre')
  for(v in c(1/16,1/4,1,4,16)) txt(c('1/16','1/4','1','4','16')[match(v,c(1/16,1/4,1,4,16))],xf(v),ay,8.5,col=muted,just='centre')
 }
}
pdf(file.path(out,'AVG_repertoire_redesign.pdf'),width=17,height=12.7,family='Helvetica',useDingbats=FALSE);draw();dev.off()
png(file.path(out,'AVG_repertoire_redesign.png'),width=17,height=12.7,units='in',res=300,type='cairo');draw();dev.off()
png(file.path(out,'preview.png'),width=2040,height=1524,res=120,type='cairo');draw();dev.off()
if(requireNamespace('svglite',quietly=TRUE)){svglite::svglite(file.path(out,'AVG_repertoire_redesign.svg'),width=17,height=12.7);draw();dev.off()}
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('Background:',nv,nt,'; modules:',nrow(res),'; displayed:',nrow(p),'\n')
cat('TA exclusions:',length(bad),'; translation-repression exclusions:',length(tr),'\n')
print(res[res$Module %in% c('Toxin-antitoxin-associated functions','Sporulation-associated functions'),])
