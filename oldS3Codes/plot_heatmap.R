######################################################
#' Modified Version of phyloseq:::plot_heatmap
my.plot_heatmap <- function(physeq, ps.ord=NULL, method="NMDS", distance="bray", label=NULL, col.type=NULL, color=NULL, normalize=TRUE, heatscale=c("#000033", "#66CCFF"), 
                            keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, na.value=NULL, trans=log_trans(4), breaks=c(1,2,4,8,16,32,60,80,100), max.label=250, ...){
  ##   Modified version of plot_heatmap function in phyloseq
  ##   1, Can take ordination-class as input if no method and distance are provided
  ##   2, Flexible Color/Label defination for axis-label by samtype and/or spetype with Color Legend
  ##   3, Customarize keepOnlyTheseTaxa, keepOnlyTheseSample, Normalization.
  
  ## Define Label and Color by Group
  sample.label <- species.label <- NULL
  sample.label <- label[label%in%sample.variables(physeq)][1]
  label <- rev(label)
  species.label <- label[label%in%rank.names(physeq)][1]
  
  sample.color <- species.color <- NULL
  samtype <- col.type[col.type%in%sample.variables(physeq)][1]
  col.type <- rev(col.type)
  spetype <- col.type[col.type%in%rank.names(physeq)][1]
  
  # Copy the approach from NeatMap and Phyloseq by doing ordination on samples, can take ordination subject for analysis
  # Get ps.ord ordination-class if ps.ord is not provided	
  if(is.null(ps.ord) & !is.null(method)){
    # Capture the NMDS iterations cat() output with capture.output
    junk <- capture.output( ps.ord <- ordinate(physeq, method, distance, ...), file=NULL)
  }
  
  # Enforce orientation and transfer otu table into matrix class
  if(!taxa_are_rows(physeq)){	physeq <- t(physeq)	}
  mot <- as(otu_table(physeq), "matrix")
  # Keep selected info and Define normalization or not. 
  if(is.numeric(normalize)){
    mot <- mot/normalize*100
  } else if(normalize){
    mot <- sapply(1:ncol(mot),function(x){mot[,x]/sum(mot[,x])}*100)
  } else {	mot <- mot	}
  colnames(mot) <- sample.names(physeq)
  
  # Group label by samtype and spetype and then trim Taxa and Sample
  if(!is.null(keepOnlyTheseTaxa)){	mot <- mot[keepOnlyTheseTaxa,]	}
  if(!is.null(keepOnlyTheseSample)){	mot <- mot[,keepOnlyTheseSample]	}
  if(!is.na(samtype)){
    tmp <- data.frame(sample_data(physeq)[colnames(mot), samtype])
    mot <- mot[, rownames(tmp)[order(tmp)]]
  }
  if(!is.na(spetype)){
    tmp <- data.frame(tax_table(physeq)[rownames(mot), spetype])
    mot <- mot[rownames(tmp)[order(tmp)], ]
  }
  
  # Initialize sample and species order vectors as NULL
  # Reorder by the angle in radial coordinates on the 2-axis plane.
  species.order <- rownames(mot)
  sample.order <- colnames(mot)
  if(!is.null(ps.ord)){
    reduction.result <- scores(ps.ord, choices=c(1, 2), display="sites")
    sample.order <- sample.names(physeq)[order(RadialCoords(reduction.result)[, "theta"])]
    sample.order <- sample.order[sample.order%in%colnames(mot)]
    mot <- mot[, sample.order]
    # Re-order species if possible
    test <- try(scores(ps.ord, choices=c(1, 2), display="species"), TRUE)
    if( class(test) != "try-error" & !is.null(test) ){
      species.reduct <- scores(ps.ord, choices=c(1, 2), display="species")
      species.order  <- species.names(physeq)[order(RadialCoords(species.reduct)[, "theta"])]
      species.order <- species.order[species.order%in%rownames(mot)]
      mot <- mot[species.order, ]
    }
  }
  
  # Alternative to melt that won't re-order your carefully-re-ordered stuff...
  adf <- data.frame(expand.grid(y=rownames(mot), x=colnames(mot)), value=as(mot, "vector"))
  
  # Now the plotting part
  p <- ggplot(adf, aes(x=as.factor(x), y=as.factor(y), fill=value)) + geom_tile() + theme(axis.text.x=element_blank(),axis.text.y=element_blank())
  p <- update_labels(p, list(fill="Abundance", y="", x=""))
  
  ## Define Color by Group
  if(!is.null(sample.label)){
    sample.color <- rep("black", length(sample.order))
    if(!is.null(samtype)){if(!is.na(samtype)){
      sam.group <- data.frame(sample_data(physeq))[sample.order, samtype]
      sample.color <- color[as.vector(sam.group)]
    }}
  }
  if(!is.null(species.label)){
    species.color <- rep("black", length(species.order))
    if(!is.null(spetype)){if(!is.na(spetype)){
      spe.group <- data.frame(tax_table(physeq))[species.order, spetype]
      species.color <- color[as.vector(spe.group)]
    }}
  }
  
  ## Axis Relabeling (Skipped if more than max.label):
  # Re-write sample-labels to some sample variable...
  if(!is.null(sample.label) & length(sample.order)<=max.label){if(!is.na(sample.label)){
    p <- p + theme(axis.text.x=element_text(size=2*treetextsize(0.1*length(sample.order)), colour=sample.color, angle = -90, hjust = 0))
    # Make a sample-named vector of the values for sample.label
    labvec <- as(getVariable(physeq, sample.label), "vector")
    names(labvec) <- sample.names(physeq)
    if(!is.null(sample.order)){
      # Re-order according to sample.order
      labvec <- labvec[sample.order]
    }
    # Add the sample.label re-labeling layer
    p <- p + scale_x_discrete(sample.label, labels=labvec)
  }}
  
  if(!is.null(species.label) & length(species.order)<=max.label){if(!is.na(species.label)){
    # Make a species-named vector of the values for species.label
    p <- p + theme(axis.text.y=element_text(size=2*treetextsize(0.1*length(species.order)), colour=species.color, hjust=1))
    labvec <- as(tax_table(physeq)[, species.label], "vector")
    names(labvec) <- species.names(physeq)
    if( !is.null(species.order) ){
      # Re-order according to species.order
      labvec <- labvec[species.order]
    }
    # Add the species.label re-labeling layer
    p <- p + scale_y_discrete(species.label, labels=labvec)
  }}
  
  ## Define na.value color
  if(is.null(na.value)){
    col.va <- round(col2rgb(heatscale[1])*1.2-col2rgb(heatscale[2])*0.2)
    col.va[col.va<0] <- 0
    col.va[col.va>255] <- 255
    na.value <- rgb(t(col.va), alpha=255, max=255)
  }
  
  # Color scale transformations
  if(!is.null(trans)){
    p <- p + scale_fill_gradient(low=heatscale[1], high=heatscale[2], trans=trans, breaks=breaks, na.value=na.value)
  } else {
    p <- p + scale_fill_gradient(low=heatscale[1], high=heatscale[2], na.value=na.value)
  }
  return(p)
}

####   Attach RadialCoords function  ####
RadialCoords <- function(pos){
  pos <- as.matrix(pos)
  n.points <- nrow(pos)
  radial <- matrix(0, n.points, 2)
  xc <- mean(pos[,1])
  yc <- mean(pos[,2])
  for(i in 1:n.points){
    radial[i,1] <- sqrt((pos[i,2]-yc)^2+(pos[i,1]-xc)^2)
    radial[i,2] <- atan2(pos[i,2]-yc, pos[i,1]-xc)
  }
  rownames(radial) <- rownames(pos)
  colnames(radial) <- c("r", "theta")
  return(radial)
}

####   Attach treetextsize function  ####
treetextsize <- function(n){
  # empirically chosen size-value calculator.
  s <- 6 * exp(-n/100)
  # enforce a floor.
  s <- ifelse(s > 0.5, s, 0.5)
  # enforce a max
  s <- ifelse(s < 4, s, 4)
  return(s)
}
