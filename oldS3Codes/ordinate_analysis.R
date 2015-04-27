######################################################
#' Basic Ordination Analysis and component composition, Utilize phyloseq ordination method
#' PCA, CCA, DCA, CA, RDA, split plot + biplot
#' PCA/CCA eigens, NMDS stress
#' 
#########   Ordination Analysis Fundemental   #########
ordinate_analysis <- function(physeq, ordination, axes=c(1,2), col.type=NULL, color=NULL, shape=NULL, label=NULL, normalize=TRUE, heatscale=c("white","steelblue"),
                              keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, trans=log_trans(4), title=NULL, td.axes=NULL, ordname, name){
  result <- list()
  
  #####   Plot Heatmap based on samtype and spetype   #####	
  #	result[["heatmap.ori"]] <- my.plot_heatmap(physeq, ps.ord=ordination, label=label, col.type=col.type, color=color, heatscale=heatscale, 
  #					keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, normalize=normalize, trans=NULL)
  #	result[[paste("heatmap",gsub("[^a-zA-Z0-9]","",trans$name),sep=".")]] <- my.plot_heatmap(physeq, ps.ord=ordination, label=label, col.type=col.type, color=color, heatscale=heatscale, 
  #					keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, normalize=normalize, trans=trans)
  
  #####  re-Organize physeq for color/shape/label #####
  if(length(col.type)==2){
    slot1 <- col.type%in%sample.variables(physeq)
    if(sum(slot1)){	sample_data(physeq)$Color <- data.frame(sample_data(physeq))[,col.type[slot1][1]]	}
    color <- rev(color)
    slot2 <- col.type%in%rank.names(physeq)
    if(sum(slot2)){	tax_table(physeq) <- cbind(tax_table(physeq),Color=as.vector(tax_table(physeq)[,col.type[slot2][1]]))	}
    col.type <- "Color"
  }
  if(length(shape)==2){
    slot1 <- shape%in%sample.variables(physeq)
    if(sum(slot1)){	sample_data(physeq)$Shape <- data.frame(sample_data(physeq))[,shape[slot1][1]]	}
    shape <- rev(shape)
    slot2 <- shape%in%rank.names(physeq)
    if(sum(slot2)){	tax_table(physeq) <- cbind(tax_table(physeq),Shape=as.vector(tax_table(physeq)[,shape[slot2][1]]))	}
    shape <- "Shape"
  }
  if(length(label)==2){
    slot1 <- label%in%sample.variables(physeq)
    if(sum(slot1)){	sample_data(physeq)$Label <- factor(data.frame(sample_data(physeq))[,label[slot1][1]])	}
    label <- rev(label)
    slot2 <- label%in%rank.names(physeq)
    if(sum(slot2)){	tax_table(physeq) <- cbind(tax_table(physeq), Label=as.vector(tax_table(physeq)[,label[slot2][1]]))	}
    label <- "Label"
  }
  
  #####   Subset Sample and Species  #####
  if(!is.null(keepOnlyTheseTaxa)){	physeq <- prune_species(keepOnlyTheseTaxa,   physeq)	}
  if(!is.null(keepOnlyTheseSample)){	physeq <- prune_samples(keepOnlyTheseSample, physeq)	}
  
  #####   Build Raw Ordination Plot DF   #####
  DF <- plot_ordination(physeq, ordination, type="split", axes=axes, color=col.type, shape=shape, label=label, justDF=TRUE)
  #	DF.3d <- plot_ordination(physeq, ordination, td.axes, type="split", color=col.type, shape=shape, justDF=TRUE)
  #####   Split Plot   #####
  result[["splitplot"]] <- ggplot(DF, aes_string(x=names(DF)[1], y=names(DF)[2], color=col.type, shape=shape, na.rm=TRUE)) + theme_bw() + 
    geom_point(aes(size=id.type), na.rm=TRUE) + scale_color_manual(values=color) + 
    scale_size_manual("type", values=c(samples=5, taxa=2)) + facet_wrap(~id.type, nrow=1)
  if(!is.null(label)){	result[["splitplot"]] <- result[["splitplot"]] + geom_text(size=4, aes_string(label=label), vjust=2, na.rm=TRUE)}
  #####  Biplot with Arrow  #####
  result[["biplot"]] <- ggplot(DF, aes_string(x=names(DF)[1], y=names(DF)[2], color=col.type, shape=shape, na.rm=TRUE)) + theme_bw() + 
    #					geom_point(subset = .(id.type=="samples"), aes(size=id.type), na.rm=TRUE) + 
    #					geom_segment(subset = .(id.type=="taxa"), aes_string(x=0, y=0, xend=names(DF)[1], yend=names(DF)[2], colour=col.type), size=1, arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) + 
    geom_point(data = subset(DF, id.type=="samples"), aes(size=id.type), na.rm=TRUE) + 
    geom_segment(data = subset(DF, id.type=="taxa"), aes_string(x=0, y=0, xend=names(DF)[1], yend=names(DF)[2], colour=col.type), size=1, arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) + 
    scale_color_manual(values=color) + scale_size_manual("type", values=c(samples=5, species=1)) + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  if(!is.null(label)){    result[["biplot"]] <- result[["biplot"]] + geom_text(size=4, aes_string(label=label), vjust=2, na.rm=TRUE)}
  #####   Stress/Eigen Plot   #####
  if(sum(class(ordination)%in%c("metaMDS","monoMDS"))){
    #####  NMDS Stress Info  #####
    k <- ncol(ordination$points)
    gof <- goodness(ordination)
    if(sum(gof)){cexp <- gof/mean(gof)*0.8} else{cexp <- 0.8}
    result[["stress"]] <- ordination$stress
    #####  Stress Plot  #####
    #		stressplot(ordination, p.col="blue", l.col="red", lwd=2, cex=0.8, ann=FALSE)
    #		mtext(side=1, text="Observed Dissimilarity", font=2, cex=0.9, line=1)
    #		mtext(side=2, text="Ordination Distance", font=2, cex=0.9, line=1)
    #		result[["stressplot"]] <- stressplot
    #####  Biplot  #####
    result[["biplot"]] <- ggplot(DF, aes_string(x=names(DF)[1], y=names(DF)[2], color=col.type, shape=shape, na.rm=TRUE)) + theme_bw() + 
      geom_point(aes(size=id.type), na.rm=TRUE) + scale_size_manual("type", values=c(samples=5, taxa=2)) + scale_color_manual(values=color)
    if(!is.null(label)){    result[["biplot"]] <- result[["biplot"]] + geom_text(size=4, aes_string(label=label), vjust=2, na.rm=TRUE)}
  } else if(sum(class(ordination)%in%"cca")){
    #####  Eigen Value  #####
    if(sum(ordination$CCA$eig)){eig <- ordination$CCA$eig} else{eig <- ordination$CA$eig}
    result[["eigen"]] <- eig
    #		result[["eigenplot"]] <- plot_eigen(vector=100*eig/sum(eig))
  } else if(class(ordination)=="dpcoa"){
    #####  Eigen Value  #####
    eig <- ordination$eig
    result[["eigen"]] <- eig
    #		result[["eigenplot"]] <- plot_eigen(vector=100*eig/sum(eig))
  } else{
    eigen <- ordination$evals
    decorana <- ordination$evals.decorana
    origin <- ordination$origin
    result[["eigen"]] <- eigen
    #		result[["eigenplot"]] <- plot_eigen
  }
  return(result)
}


###   vegan:::scores.default   ###
scores.dpcoa <- function(x, choices=NULL, display="sites", ...){
  ifelse(display=="species", coords <- x$l1, coords <- x$l2)
  if( is.null(choices) ){
    choices <- colnames(coords)
  }
  return( coords[, choices] )
}
