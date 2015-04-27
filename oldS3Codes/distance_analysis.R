######################################################
#' Basic Distance Analysis and beta diversity
#' support all vegdist and unifrac
#' Mental like Test
#' Intra-group vs. Inter-group distance average
#' PCoA, NMDS, 2d + 3d
#' 
#########   Distance Analysis Fundemental   #########
distance_analysis <- function(physeq, distmatrix, axes=c(1,2), samtype=NULL, color=NULL, shape=NULL, label=NULL, normalize=TRUE, heatscale=c("white","steelblue"), max.dist=0.3, 
                              keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, trans=log_trans(4), title=NULL, flip=FALSE, td.axes=NULL, distname=NULL, name=NULL){
  dist <- list()
  group <- data.frame(sample_data(physeq))[,samtype]
  
  #######   Beta diversity Mental like Test   #######
  com <- combn(levels(group),2)
  beta.test <- matrix(0, ncol(com)+1, 3)
  rownames(beta.test) <- c("all", paste(com[1,], com[2,], sep=" vs "))
  colnames(beta.test) <- c("Anosim","MRPP","perMANOVA")
  
  anosim <- anosim(distmatrix, group, permutations=1000)
  mrpp <- mrpp(distmatrix, group, permutations=1000)
  permanova <- adonis(distmatrix~group, permutations=1000)
  beta.test[1,] <- c(anosim$signif, mrpp$Pvalue, permanova$aov.tab$"Pr(>F)"[1])
  
  for(i in 1:ncol(com)){
    sub_index <- group==com[1,i]|group==com[2,i]
    sub_dist <- as.dist(as.matrix(distmatrix)[sub_index, sub_index])
    sub_group <- factor(group[sub_index])
    beta.test[i+1,1] <- anosim(sub_dist, sub_group, permutations=1000)$signif
    beta.test[i+1,2] <- mrpp(sub_dist, sub_group, permutations=1000)$Pvalue
    beta.test[i+1,3] <- adonis(sub_dist~sub_group, permutations=1000)$aov.tab$"Pr(>F)"[1]
  }
  
  #######   Beta diversity Intra/Inter-Group Test   #######
  com <- cbind(rbind(levels(group),levels(group)), com)
  for(i in 1:ncol(com)){
    tmp <- as.vector(as.matrix(distmatrix)[group==com[1,i],group==com[2,i]])
    dist[[i]] <- tmp[tmp!=0]
  }
  distance <- unlist(dist)
  compare <- factor(rep(c(1:ncol(com)), sapply(dist,function(x){length(x)})))
  levels(compare) <- paste(com[1,],com[2,], sep=" vs. ")[as.numeric(levels(compare))]
  plotdata <- data.frame(compare, distance)
  Dsummary <- summaryBy(distance~compare, plotdata, FUN=c(mean, sd, length))
  signif <- try(pairwise.t.test(distance, compare, pool.sd=FALSE, p.adj="none"))
  
  #######   boxplot for Intra/Inter-Group Test   #######
  pbox <- ggplot(plotdata, aes(compare, distance)) + geom_boxplot() + theme_bw()
  if(flip==FALSE){
    pbox <- pbox + theme(axis.text.x=element_text(size=8, angle=90), axis.text.y=element_text(size=6),
                         panel.grid.major=element_blank(), panel.grid.minor=element_line(colour='grey', linetype='dashed'), 
                         legend.key = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())
  }
  if(flip==TRUE){
    pbox <- pbox + coord_flip() + theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=8),
                                        panel.grid.major=element_blank(), panel.grid.minor=element_line(colour='grey', linetype='dashed'), 
                                        legend.key = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())
  }
  
  #######   barplot   #######
  pbar <- ggplot(Dsummary, aes(compare, distance.mean)) + geom_bar(stat="identity", colour="black", fill="white") + theme_bw() + 
    geom_errorbar(aes(ymin=distance.mean-distance.sd/sqrt(distance.length), ymax=distance.mean+distance.sd/sqrt(distance.length)), width=.25, position=position_dodge(.9))
  
  if(flip==FALSE){
    pbar <- pbar + theme(axis.text.x=element_text(size=8, angle=90), axis.text.y=element_text(size=6),
                         legend.key = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())
  }
  if(flip==TRUE){
    pbar <- pbar + coord_flip() + theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=8),
                                        legend.key = element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank())
  }
  
  #######   PCoA and NMDS   #######
  PCoA <- ordinate(physeq, method="PCoA", distance=distmatrix)
  NMDS <- ordinate(physeq, method="NMDS", distance=distmatrix)
  ##   Notification: ordinate bug: line 239 and 243, failed to pass other par to metaMDS
  ##   return(metaMDS(ps.dist)) should be: return(metaMDS(ps.dist, ...))
  ##   If fixed, change code to: ordinate(physeq, method="NMDS", distance=distmatrix, k=3)
  if(!is.null(td.axes)){	NMDS.3d <- metaMDS(distmatrix, k=max(td.axes))	}
  
  #######   Plot Heatmap   #######
  #	heat.pcoa <- my.plot_heatmap(physeq, PCoA, label="ID", species.label="Taxa")
  #	heat.nmds <- my.plot_heatmap(physeq, method="NMDS", distmatrix, sample.label="ID", species.label="Taxa")
  
  #######   PCoA and NMDS scatter Plot   #######
  if(!is.null(keepOnlyTheseTaxa)){	physeq <- prune_species(keepOnlyTheseTaxa,   physeq)	}
  if(!is.null(keepOnlyTheseSample)){	physeq <- prune_samples(keepOnlyTheseSample, physeq)	}
  
  #######   2D scatter Plot   #######
  slot1 <- label%in%sample.variables(physeq)
  if(sum(slot1)){label <- label[slot1][1]}
  ppcoa <- plot_ordination(physeq, PCoA, type="samples", axes=axes, color=samtype, shape=shape) + geom_point(size=5) + 
    scale_color_manual(values=color) + geom_text(size=4, aes_string(label=label), vjust=2, na.rm=TRUE) + 
    #			scale_shape_manual(values=shape) + 
    theme_bw() + scale_x_continuous(name="PC1") + scale_y_continuous(name="PC2")
  pnmds <- plot_ordination(physeq, NMDS, type="samples", axes=axes, color=samtype, shape=shape) + geom_point(size=5) + 
    scale_color_manual(values=color) + geom_text(size=4, aes_string(label=label), vjust=2, na.rm=TRUE) + 
    #			scale_shape_manual(values=shape) + 
    theme_bw() + scale_x_continuous(name="MDS1") + scale_y_continuous(name="MDS2")
  
  #######   PCoA 3D Plot   #######
  if(!is.null(td.axes)){
    plotdata <- PCoA$vectors[,td.axes]
    if(!is.null(keepOnlyTheseSample)){
      plotdata <- plotdata[keepOnlyTheseSample,]
      group <- factor(as.vector(group[keepOnlyTheseSample]))
    }	
    cloud3d(plotdata, labels=group, filename=paste(name,"dist",distname,"pcoa3d","wrl",sep="."), type="vrml", pointstyle="s", 
            metalabels=rownames(plotdata), cols=color, scalefac=4, autoscale="independent", lab.axis=paste("PC",td.axes,sep=""), 
            col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
    plotdata <- NMDS.3d$points[,td.axes]
    if(!is.null(keepOnlyTheseSample)){
      plotdata <- plotdata[keepOnlyTheseSample,]
      group <- factor(as.vector(group[keepOnlyTheseSample]))
    }
    cloud3d(plotdata, labels=group, filename=paste(name,"dist",distname,"nmds3d","wrl",sep="."), type="vrml", pointstyle="s", 
            metalabels=rownames(plotdata), cols=color, scalefac=4, autoscale="independent", lab.axis=paste("MDS",td.axes,sep=""), 
            col.axis="darkblue", showaxis=TRUE, col.lab="darkgray", col.bg="white", cex.lab=1, htmlout=NULL, hwidth=1200, hheight=800, showlegend=TRUE)
  }
  
  #######   IGraph Network with max.dist.list  #######
  igraph <- list()
  if(!is.null(max.dist)){
    for(i in 1:length(max.dist)){
      dist <- max.dist[i]
      ig <- make_network(physeq, type="samples", distance=distmatrix, max.dist=dist, keep.isolates=TRUE)
      if(length(ig[[3]])){
        igraph[[paste("igraph", round(dist,2), sep="_")]] <- plot_network(ig, physeq, color=samtype, shape=shape, line_weight=0.4, label=NULL)
      }
    }
  }
  
  result <- list()
  result[["Dview"]] <- Dsummary
  result[["beta.test"]] <- data.frame(beta.test)
  result[["Signif"]] <- signif
  result[["PCoA"]] <- PCoA
  result[["NMDS"]] <- NMDS
  result[["boxplot"]] <- pbox
  result[["barplot"]] <- pbar
  result[["pcoaplot"]] <- ppcoa
  result[["nmdsplot"]] <- pnmds
  #	result[["heatmap.pcoa"]] <- heat.pcoa
  #	result[["heatmap.nmds"]] <- heat.nmds
  result[["igraph"]] <- igraph
  
  return(result)
}

####   Attach Plot3d from corrplot.r   ####
rgl_plot3d <- function(plotdata, group=NULL, color, lab=NULL, type="s", size=1){
  library(rgl)
  if(!is.null(group)){	color <- color[as.vector(group)]	}
  if(length(type)>1){	type <- type[as.vector(group)]	}
  if(length(size)>1){	size <- size[as.vector(group)]	}
  plot3d(plotdata, col=color, type=type, size=size, xlab=lab[1], ylab=lab[2], zlab=lab[3])
  detach("package:rgl")
}

