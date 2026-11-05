# Three key AVG phylogenies: Fhs, FolD and NrdG.
# Run: Rscript plot_AVG_trees.R
# Dependencies: install.packages(c("ape", "digest"))
# No Python, ggtree, external image files or network access is used by this script.
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(grid))
if (!requireNamespace("digest", quietly=TRUE)) stop("Install digest first")
args <- commandArgs(FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script)) {
  ROOT <- dirname(normalizePath(script[1]))
} else {
  ofiles <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  ROOT <- if (length(ofiles)) dirname(normalizePath(tail(ofiles,1)[[1]])) else getwd()
}
DATA <- file.path(ROOT, "data")
OUT <- file.path(ROOT, "figures")
AUDIT <- file.path(ROOT, "audit")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
dir.create(AUDIT, recursive=TRUE, showWarnings=FALSE)

# User-editable appearance. All font sizes are points at the exported canvas size.
COLORS <- c(XVIR="#C65C32", PIBAC="#40875D", NCBI="#377FA6", CUR="#655296")
FOCUS <- character(0) # Uniform display: no individually starred candidates.
LABEL_PT <- 11
PANEL_LABEL_PT <- 10.5
BRANCH_COL <- "#69757D"
STRONG_COL <- "#263D49"
SH_MIN <- 80
UF_MIN <- 95
PNG_DPI <- 180
# Overview hides individual NCBI/PIBAC text only; their tips remain drawn.
# Set TRUE to print every reference label in the overview as well.
OVERVIEW_ALL_LABELS <- FALSE

read_tsv <- function(p) read.delim(p, check.names=FALSE, stringsAsFactors=FALSE,
                                  quote="", comment.char="", na.strings=NULL,
                                  colClasses="character")
write_tsv <- function(x,p) write.table(x,p,sep="\t",row.names=FALSE,quote=FALSE,na="")
tree_structure <- function(tr) {
  nt <- Ntip(tr); nn <- nt+tr$Nnode
  ch <- vector("list", nn)
  for (i in seq_len(nrow(tr$edge))) ch[[tr$edge[i,1]]] <- c(ch[[tr$edge[i,1]]],tr$edge[i,2])
  root <- setdiff(tr$edge[,1],tr$edge[,2])[1]
  des <- vector("list",nn)
  get_desc <- function(v) {
    z <- if (v<=nt) v else unlist(lapply(ch[[v]],get_desc),use.names=FALSE)
    des[[v]] <<- z; z
  }
  get_desc(root)
  list(children=ch, descendants=des, root=root, nt=nt, nn=nn)
}
split_keys <- function(tr) {
  s <- tree_structure(tr); allnames <- sort(tr$tip.label)
  vapply(s$descendants,function(ids) {
    a <- sort(tr$tip.label[ids]); b <- setdiff(allnames,a)
    sa <- paste(a,collapse=";"); sb <- paste(b,collapse=";")
    if (length(a)<length(b)) sa else if (length(a)>length(b)) sb else min(sa,sb)
  },character(1))
}
support_matrix <- function(tr) {
  z <- matrix(NA_real_,Ntip(tr)+tr$Nnode,2,dimnames=list(NULL,c("SH_aLRT","UFBoot")))
  if (!is.null(tr$node.label)) for (j in seq_len(tr$Nnode)) {
    x <- suppressWarnings(as.numeric(strsplit(tr$node.label[j],"/",fixed=TRUE)[[1]]))
    if (length(x)==2) z[Ntip(tr)+j,] <- x
  }
  z
}
strong_nodes <- function(tr) {
  s <- support_matrix(tr)
  !is.na(s[,1]) & !is.na(s[,2]) & s[,1]>=SH_MIN & s[,2]>=UF_MIN
}

# Midpoint rooting directly on the weighted, undirected tree graph.
# Support labels are properties of splits: they are remapped by split identity.
# This avoids silently moving a support label to a different edge on rerooting.
midpoint_verified <- function(tr) {
  stopifnot(!is.null(tr$edge.length), all(is.finite(tr$edge.length)),all(tr$edge.length>=0),
            !anyDuplicated(tr$tip.label))
  old <- tr; st <- tree_structure(old); nt <- st$nt; nn <- st$nn
  oldkeys <- split_keys(old)
  oldsplits <- sort(unique(oldkeys[old$edge[,2]]))
  sm <- support_matrix(old)
  internal <- (nt+1):nn
  supported <- internal[!is.na(sm[internal,1]) & internal!=st$root]
  supportmap <- setNames(old$node.label[supported-nt],oldkeys[supported])
  d <- cophenetic.phylo(old)
  ends <- arrayInd(which.max(d),dim(d))[1,]
  half <- d[ends[1],ends[2]]/2
  adj <- vector("list",nn+1)
  for (i in seq_len(nrow(old$edge))) {
    a <- old$edge[i,1]; b <- old$edge[i,2]; w <- old$edge.length[i]
    adj[[a]] <- rbind(adj[[a]],c(b,w)); adj[[b]] <- rbind(adj[[b]],c(a,w))
  }
  prev <- rep(NA_integer_,nn+1); plen <- numeric(nn+1)
  prev[ends[1]] <- 0L; stack <- ends[1]
  while (length(stack)) {
    v <- tail(stack,1); stack <- head(stack,-1)
    if (v==ends[2]) break
    for (i in seq_len(nrow(adj[[v]]))) {
      u <- as.integer(adj[[v]][i,1])
      if (is.na(prev[u])) {prev[u]<-v; plen[u]<-adj[[v]][i,2]; stack<-c(stack,u)}
    }
  }
  path <- matrix(numeric(),ncol=3); v <- ends[2]
  while (prev[v]!=0) {path<-rbind(c(prev[v],v,plen[v]),path);v<-prev[v]}
  walked <- 0
  for (i in seq_len(nrow(path))) {
    a<-path[i,1];b<-path[i,2];w<-path[i,3]
    if (walked+w>=half) {
      left<-half-walked;right<-w-left
      if (abs(left)<1e-13) newroot<-a else if (abs(right)<1e-13) newroot<-b else {
        newroot<-nn+1L
        adj[[a]]<-adj[[a]][adj[[a]][,1]!=b,,drop=FALSE]
        adj[[b]]<-adj[[b]][adj[[b]][,1]!=a,,drop=FALSE]
        adj[[a]]<-rbind(adj[[a]],c(newroot,left))
        adj[[b]]<-rbind(adj[[b]],c(newroot,right))
        adj[[newroot]]<-rbind(c(a,left),c(b,right))
      }
      break
    }
    walked<-walked+w
  }
  ids<-integer(nn+1);ids[seq_len(nt)]<-seq_len(nt);counter<-nt
  edge<-matrix(integer(),ncol=2);lens<-numeric()
  orient <- function(v,parent=0L) {
    if (v>nt) {counter<<-counter+1L;ids[v]<<-counter}
    for (i in seq_len(nrow(adj[[v]]))) {
      u<-as.integer(adj[[v]][i,1]); if (u==parent) next
      orient(u,v)
      edge<<-rbind(edge,c(ids[v],ids[u]));lens<<-c(lens,adj[[v]][i,2])
    }
  }
  orient(newroot)
  tr<-list(edge=edge,edge.length=lens,tip.label=old$tip.label,Nnode=counter-nt)
  class(tr)<-"phylo";storage.mode(tr$edge)<-"integer"
  tr<-reorder.phylo(tr,"cladewise")
  newkeys<-split_keys(tr);tr$node.label<-rep("",tr$Nnode)
  for (v in (nt+1):(nt+tr$Nnode)) {
    if (v==nt+1) next
    val<-supportmap[newkeys[v]]
    if (length(val) && !is.na(val)) tr$node.label[v-nt]<-unname(val)
  }
  tr<-ladderize(tr,right=TRUE)
  nd<-cophenetic.phylo(tr)[rownames(d),colnames(d)]
  err<-max(abs(d-nd))
  stopifnot(err<1e-9,identical(oldsplits,sort(unique(split_keys(tr)[tr$edge[,2]]))))
  keys<-split_keys(tr);ns<-support_matrix(tr);nroot<-tree_structure(tr)$root
  for (v in (nt+1):(nt+tr$Nnode)) if (!is.na(ns[v,1]) && v!=nroot) {
    stopifnot(identical(unname(supportmap[keys[v]]),tr$node.label[v-nt]))
  }
  depth<-node.depth.edgelength(tr)
  stopifnot(all(abs(depth[ends]-half)<1e-9),
            all(names(supportmap)%in%keys[which(!is.na(ns[,1]))]))
  list(tree=tr,validation=data.frame(tips=nt,distance_pairs=nt*(nt-1)/2,
       max_distance_change=err,diameter=2*half,root_to_endpoint=half,
       topology_preserved=TRUE,support_splits_preserved=length(supportmap)))
}

build_model <- function(tr,metadata,title) {
  st<-tree_structure(tr); sm<-support_matrix(tr); strong<-strong_nodes(tr)
  meta<-metadata[match(tr$tip.label,metadata$tree_id),]
  stopifnot(!anyNA(meta$tree_id),all(meta$source%in%names(COLORS)))
  collapsed<-list();members<-list()
  visit<-function(v) {
    ids<-st$descendants[[v]];r<-meta[ids,,drop=FALSE]
    if (length(ids)>=2 && all(r$source%in%c("NCBI","PIBAC")) &&
        length(unique(r$phylum))==1 && !is.na(r$phylum[1]) && nzchar(r$phylum[1]) && strong[v]) {
      cid<-sprintf("C%02d",length(collapsed)+1)
      collapsed[[as.character(v)]]<<-list(label=paste0(r$phylum[1]," (n=",length(ids),")"),
                                        members=ids,id=cid)
      members[[cid]]<<-data.frame(clade_id=cid,tree_id=r$tree_id,source=r$source,
                          phylum=r$phylum,SH_aLRT=unname(sm[v,1]),UFBoot=unname(sm[v,2]))
    } else for (u in st$children[[v]]) visit(u)
  }
  visit(st$root)
  list(tree=tr,st=st,meta=meta,title=title,collapsed=collapsed,
       members=if(length(members)) do.call(rbind,members) else NULL,
       support=sm,strong=strong)
}
geometry <- function(m,full=FALSE,topology=FALSE) {
  cc<-if(full) list() else m$collapsed
  visible<-integer();tips<-integer();x<-y<-rep(NA_real_,m$st$nn);x[m$st$root]<-0
  el<-setNames(m$tree$edge.length,m$tree$edge[,2])
  walk<-function(v) {
    visible<<-c(visible,v)
    if (v<=m$st$nt || as.character(v)%in%names(cc)) {
      tips<<-c(tips,v);y[v]<<-length(tips)
    } else {
      for (u in m$st$children[[v]]) {x[u]<<-x[v]+if(topology) 1 else el[as.character(u)];walk(u)}
      yy<-y[m$st$children[[v]]];y[v]<<-mean(range(yy))
    }
  }
  walk(m$st$root)
  list(visible=visible,tips=tips,x=x,y=y,cc=cc)
}
node_label <- function(m,v,cc,overview=FALSE) {
  if(as.character(v)%in%names(cc)) return(cc[[as.character(v)]]$label)
  if (overview && !OVERVIEW_ALL_LABELS && !m$meta$source[v]%in%c("XVIR","CUR")) return("")
  label<-m$meta$label[v]
  if(overview && m$meta$source[v]=="CUR") label<-sub(" \\|.*$","",label)
  original<-sub("\\|.*$","",m$meta$original[v])
  if(original%in%FOCUS) label<-paste0("* ",label)
  label
}
node_color <- function(m,v) if(v<=m$st$nt) unname(COLORS[m$meta$source[v]]) else "#567A82"
ln <- function(x,y,col=BRANCH_COL,lwd=.6,lty=1) grid.lines(x,y,default.units="native",gp=gpar(col=col,lwd=lwd,lty=lty))
pt <- function(x,y,col,size=1) grid.points(x,y,default.units="native",pch=16,size=unit(size,"mm"),gp=gpar(col=col))

draw_circle <- function(m,overview=FALSE) {
  g<-geometry(m,topology=TRUE); n<-length(g$tips)
  a<-(8+344*(g$y-1)/max(1,n-1))*pi/180
  radius<-max(g$x,na.rm=TRUE)+1;r<-g$x/radius;lr<-1.06
  # A physical square viewport preserves circular geometry in every panel.
  lim<-if(overview) 1.65 else 1.95
  pushViewport(viewport(width=unit(1,"snpc"),height=unit(1,"snpc"),
                        xscale=c(-lim,lim),yscale=c(-lim,lim)))
  for(v in g$visible) {
    if(v%in%g$tips) {
      col<-node_color(m,v);end<-r[v]
      if(as.character(v)%in%names(g$cc)) {
        da<-(344/max(1,n-1))*pi/180*.3
        grid.polygon(c(r[v]*cos(a[v]),cos(a[v]-da),cos(a[v]+da)),
                     c(r[v]*sin(a[v]),sin(a[v]-da),sin(a[v]+da)),default.units="native",
                     gp=gpar(fill="#DCE8EA",col=col,lwd=.65));end<-1
      } else pt(end*cos(a[v]),end*sin(a[v]),col,if(overview) .85 else 1.2)
      label<-node_label(m,v,g$cc,overview)
      if(nzchar(label)) {
        ln(c(end,lr-.02)*cos(a[v]),c(end,lr-.02)*sin(a[v]),"#CCD2D6",.4,3)
        deg<-a[v]*180/pi;left<-deg>90 && deg<270
        bold<-v<=m$st$nt && m$meta$source[v]%in%c("XVIR","CUR")
        grid.text(label,lr*cos(a[v]),lr*sin(a[v]),default.units="native",
                  rot=if(left) deg+180 else deg,just=if(left) "right" else "left",
                  gp=gpar(col=col,fontsize=if(overview) PANEL_LABEL_PT else LABEL_PT,
                          fontface=if(bold) "bold" else "plain"))
      }
    } else {
      ch<-m$st$children[[v]];aa<-seq(min(a[ch]),max(a[ch]),length.out=60)
      ln(r[v]*cos(aa),r[v]*sin(aa))
      for(u in ch) ln(c(r[v],r[u])*cos(a[u]),c(r[v],r[u])*sin(a[u]))
      if(v!=m$st$root && !is.na(m$support[v,1])) {
        pt(r[v]*cos(a[v]),r[v]*sin(a[v]),if(m$strong[v]) STRONG_COL else "#CDD1D4",if(m$strong[v]) .9 else .55)
        if(!overview && m$strong[v] && !v%in%m$st$children[[m$st$root]]) {
          grid.text(paste(m$support[v,],collapse="/"),unit(r[v]*cos(a[v]),"native")+unit(1,"mm"),
                    unit(r[v]*sin(a[v]),"native")+unit(1,"mm"),just=c("left","bottom"),
                    gp=gpar(fontsize=7,col=STRONG_COL))
        }
      }
    }
  }
  popViewport()
}
draw_rect <- function(m,full=FALSE) {
  g<-geometry(m,full=full);n<-length(g$tips)
  depth<-node.depth.edgelength(m$tree);maxx<-max(depth);gap<-maxx*.02
  pushViewport(viewport(xscale=c(-.02*maxx,1.80*maxx),yscale=c(n+2,-1)))
  for(v in g$visible) {
    if(v%in%g$tips) {
      col<-node_color(m,v);end<-g$x[v]
      if(as.character(v)%in%names(g$cc)) {
        end<-max(depth[m$st$descendants[[v]]])
        grid.polygon(c(g$x[v],end,end),c(g$y[v],g$y[v]-.3,g$y[v]+.3),default.units="native",
                     gp=gpar(fill="#DCE8EA",col=col,lwd=.65))
      } else pt(end,g$y[v],col,1.1)
      grid.text(node_label(m,v,g$cc),end+gap,g$y[v],default.units="native",just="left",
                gp=gpar(fontsize=LABEL_PT,col=col,
                        fontface=if(v<=m$st$nt && m$meta$source[v]%in%c("XVIR","CUR")) "bold" else "plain"))
    } else {
      ch<-m$st$children[[v]];ln(rep(g$x[v],2),range(g$y[ch]),lwd=.7)
      for(u in ch) ln(c(g$x[v],g$x[u]),rep(g$y[u],2),lwd=.7)
    }
    if(v>m$st$nt && v!=m$st$root && !is.na(m$support[v,1])) {
      pt(g$x[v],g$y[v],if(m$strong[v]) STRONG_COL else "#CDD1D4",.7)
      if(m$strong[v]) grid.text(paste(m$support[v,],collapse="/"),
                         unit(g$x[v],"native")-unit(.7,"mm"),unit(g$y[v],"native")+unit(.8,"mm"),
                         just=c("right","bottom"),gp=gpar(fontsize=7.5,col=STRONG_COL))
    }
  }
  scale<-10^floor(log10(maxx*.15));if(maxx*.15/scale>=2)scale<-scale*2
  ln(c(0,scale),rep(n+1.05,2),"#333333",1)
  grid.text(paste(scale,"substitutions/site"),scale/2,n+1.65,default.units="native",gp=gpar(fontsize=9))
  popViewport()
}
draw_legend <- function() {
  labs<-c("Viral candidate","PIBAC","NCBI","Curated","Strong support")
  cols<-c(COLORS,STRONG_COL);xs<-seq(.12,.88,length.out=5)
  for(i in seq_along(labs)) {
    grid.points(xs[i]-.02,.5,pch=16,size=unit(1.6,"mm"),gp=gpar(col=cols[i]))
    grid.text(labs[i],xs[i],.5,just="left",gp=gpar(fontsize=12,col="#263238"))
  }
}
save_formats <- function(stem,width,height,draw,png=TRUE) {
  # Cairo embeds fonts and supports portable Unicode text and PDF output.
  # Render to a temporary file, close the device, then copy complete bytes.
  # This also avoids partially buffered files on synchronized filesystems.
  commit <- function(tmp,dest) {
    n<-file.info(tmp)$size
    blob<-readBin(tmp,"raw",n=n);stopifnot(length(blob)==n)
    writeBin(blob,dest);stopifnot(file.info(dest)$size==n);unlink(tmp)
  }
  tmp<-tempfile(fileext=".pdf")
  cairo_pdf(tmp,width=width,height=height,family="sans",bg="white")
  draw();dev.off()
  commit(tmp,paste0(stem,".pdf"))
  tmp<-tempfile(fileext=".svg")
  svg(tmp,width=width,height=height,family="sans",bg="white")
  draw();dev.off()
  commit(tmp,paste0(stem,".svg"))
  if(png) {
    tmp<-tempfile(fileext=".png")
    grDevices::png(tmp,width=width,height=height,units="in",res=PNG_DPI,type="cairo",bg="white")
    draw();dev.off()
    commit(tmp,paste0(stem,".png"))
  }
}
draw_single <- function(m,mode) {
  grid.newpage()
  grid.text(m$title,.035,.985,just=c("left","top"),gp=gpar(fontsize=16,fontface="bold"))
  pushViewport(viewport(x=.5,y=.963,width=.7,height=.025));draw_legend();popViewport()
  pushViewport(viewport(x=.5,y=.475,width=.94,height=.92))
  if(mode=="circular") draw_circle(m) else draw_rect(m,mode=="complete")
  popViewport()
}
draw_panels <- function(models,ncol=2) {
  nrow<-ceiling(length(models)/ncol)
  grid.newpage()
  pushViewport(viewport(x=.5,y=.51,width=.98,height=.965,
                        layout=grid.layout(nrow,ncol)))
  for(i in seq_along(models)) {
    pushViewport(viewport(layout.pos.row=ceiling(i/ncol),layout.pos.col=(i-1)%%ncol+1))
    grid.text(paste(LETTERS[i],models[[i]]$title),.025,.98,just=c("left","top"),
              gp=gpar(fontsize=19,fontface="bold",col="#202C34"))
    pushViewport(viewport(x=.5,y=.47,width=.96,height=.90))
    draw_circle(models[[i]],overview=TRUE);popViewport(2)
  }
  popViewport()
  pushViewport(viewport(x=.5,y=.014,width=.8,height=.02));draw_legend();popViewport()
}

# Load the three current trees (Fhs uses five high-quality viral candidates).
# Every output is reconstructed from raw Newick.
manifest<-read_tsv(file.path(DATA,"source","tree_manifest.tsv"))
order<-read_tsv(file.path(DATA,"tree_order.tsv"))
expected<-read_tsv(file.path(DATA,"expected_collapse_summary.tsv"))
expected_members<-read_tsv(file.path(DATA,"expected_collapsed_members.tsv"))
models<-list();validation<-list();summaries<-list();members<-list()
for(i in seq_len(nrow(order))) {
  name<-order$tree[i];section<-order$section[i]
  rec<-manifest[tools::file_path_sans_ext(basename(manifest$file))==name,,drop=FALSE]
  stopifnot(nrow(rec)==1)
  path<-file.path(DATA,"source",rec$file)
  stopifnot(digest::digest(path,algo="sha256",file=TRUE,serialize=FALSE)==rec$sha256)
  raw<-read.tree(path);stopifnot(Ntip(raw)==rec$n_tips)
  rooted<-midpoint_verified(raw)
  m<-build_model(rooted$tree,read_tsv(file.path(DATA,"metadata",paste0(name,".tsv"))),order$title[i])
  models[[name]]<-m
  hidden<-if(is.null(m$members)) 0 else nrow(m$members)
  summary<-data.frame(tree=name,section=section,original_tips=Ntip(raw),
           collapsed_clades=length(m$collapsed),hidden_reference_tips=hidden,
           displayed_units=Ntip(raw)-hidden+length(m$collapsed))
  exp<-expected[expected$tree==name,]
  stopifnot(summary$collapsed_clades==exp$collapsed_clades,
            summary$displayed_units==exp$displayed_units)
  # Compare exact member partitions, not just collapse counts, with previous output.
  partition<-function(df) if(is.null(df)||!nrow(df)) character() else
    unname(sort(vapply(split(df$tree_id,df$clade_id),function(x) paste(sort(x),collapse=";"),character(1))))
  stopifnot(identical(partition(m$members),partition(expected_members[expected_members$tree==name,])))
  validation[[name]]<-cbind(tree=name,rooted$validation)
  summaries[[name]]<-summary
  if(!is.null(m$members)) members[[name]]<-cbind(tree=name,m$members)
  dest<-file.path(OUT,section,name);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  write.tree(m$tree,file=file.path(dest,paste0(name,".midpoint.treefile")),digits=15)
  for(mode in c("collapsed","circular","complete")) {
    width<-if(mode=="circular") if(Ntip(m$tree)>130) 22 else 17 else 17
    height<-if(mode=="circular") width else max(7,length(geometry(m,mode=="complete")$tips)*.25+1.4)
    save_formats(file.path(dest,paste0(name,"_",mode)),width,height,
                 function() draw_single(m,mode),png=mode!="complete")
  }
  message(name,": ",Ntip(raw)," tips; ",length(m$collapsed)," collapsed clades; validated")
}
write_tsv(do.call(rbind,validation),file.path(AUDIT,"rooting_validation.tsv"))
write_tsv(do.call(rbind,summaries),file.path(AUDIT,"collapse_summary.tsv"))
write_tsv(do.call(rbind,members),file.path(AUDIT,"collapsed_members.tsv"))
main<-models[order$tree[order$main_panel==1]]
save_formats(file.path(OUT,"AVG3_midpoint_combined"),27,9.5,function() draw_panels(main,3))
capture.output(sessionInfo(),file=file.path(AUDIT,"R_sessionInfo.txt"))
message("Finished: ",OUT)
