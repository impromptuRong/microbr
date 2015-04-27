######################################################
#' Basic Multivariate Analysis
#' refine otu table with smaller fraction random number for some 0 case
#' taxa_analysis.R
#' rarefaction curves.
#' distance_analysis.R
#' ordinate_analysis.R
#' ...
#' 
SA_uni <- function(physeq, samtype, spetype=NULL, col.sam=NULL, col.spe=NULL, shape=NULL, label=NULL, facet_formula=facet_formula, keepOnlyTheseTaxa=NULL, keepOnlyTheseSample=NULL, 
                   distlist=c("bray","jaccard"), max.dist=0.3, threshold=0, axes=c(1,2), seq.depth=100, heatscale=heatscale, trans=trans, name, flip=FALSE, td.axes=NULL, rarefy=FALSE, permu=100){
  #######    phyloseq detach   ########
  data <- data.frame(otu_table(physeq))
  tmp1 <- data[rowSums(data)!=0,]
  tmp2 <- data[rowSums(data)==0,]
  if(nrow(tmp2)){
    for(i in 1:nrow(tmp2)){
      n <- sample(2:10)
      if(ncol(tmp2)<10){n <- sample(2:ncol(tmp2))}
      tmp2[i,sample(1:ncol(tmp2))[1:n[1]]] <- 0.1/10^n[c(2:n[1],1)]
    }
    tmp2[is.na(tmp2)] <- 0
    data <- rbind(tmp1,tmp2)
  }
  otu_table(physeq) <- otu_table(data, taxa_are_rows=FALSE)
  # 	sample_data <- data.frame(sample_data(physeq))
  # 	taxonomy <- data.frame(tax_table(physeq))
  # 	ID <- sample.names(physeq)
  # 	feature <- species.names(physeq)
  # 	group <- factor(sample_data[,samtype])
  # 	type <- factor(taxonomy[,spetype])
  
  #####   Define scale for heatmap and Abundance plot  #####
  hei <- nspecies(physeq)
  wid <- nsamples(physeq)
  if(!is.null(keepOnlyTheseTaxa)){	hei <- sum(keepOnlyTheseTaxa %in% species.names(physeq))	}
  if(!is.null(keepOnlyTheseSample)){	wid <- sum(keepOnlyTheseSample %in% sample.names(physeq))	}
  
  #######    Taxanomy Analysis   #######
  taxa.info <- taxa_analysis(physeq, samtype=samtype, color=col.sam, facet_formula=facet_formula, threshold=threshold, normalize=seq.depth, 
                             keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, OTUpoints=FALSE, labelOTUs=FALSE, rarefy=rarefy, permu=permu)
  #####   Abundane plot   #####
  png(file=paste(name, "tax.boxplot", "png", sep="."), width=length(col.sam)*hei*50+500, height=2000, res=300)
  print(taxa.info$tax.boxplot)
  dev.off()
  png(file=paste(name, "tax.barplot", "png", sep="."), width=length(col.sam)*hei*50+500, height=2000, res=300)
  print(taxa.info$tax.barplot)
  dev.off()
  png(file=paste(name, "tax.richness.boxplot", "png", sep="."), width=5000, height=2000, res=300)
  print(taxa.info$tax.richness$boxplot)
  dev.off()
  png(file=paste(name, "tax.richness.barplot", "png", sep="."), width=5000, height=2000, res=300)
  print(taxa.info$tax.richness$barplot)
  dev.off()
  png(file=paste(name, "tax.alphaindex.boxplot", "png", sep="."), width=5000, height=2000, res=300)
  print(taxa.info$tax.alphaindex$boxplot)
  dev.off()
  png(file=paste(name, "tax.alphaindex.barplot", "png", sep="."), width=5000, height=2000, res=300)
  print(taxa.info$tax.alphaindex$barplot)
  dev.off()
  write.csv(taxa.info$tax.diversity, paste(name,".diversity.csv",sep=""))
  write.csv(taxa.info$tax.diver.test$comb_info, paste(name,".diver.test.csv",sep=""))
  
  if(rarefy){
    #####   Rarefaction Curve   #####
    png(file=paste(name, "tax.rarefy", samtype, "png", sep="."), width=6500, height=4000, res=300)
    print(taxa.info$tax.rarefy$group)
    dev.off()
    png(file=paste(name, "tax.sample.shannon", "png", sep="."), width=ceiling(sqrt(wid))*500+500, height=ceiling(sqrt(wid))*500+500, res=300)
    print(taxa.info$tax.rarefy$sample$shannon)
    dev.off()
    png(file=paste(name, "tax.sample.simpson", "png", sep="."), width=ceiling(sqrt(wid))*500+500, height=ceiling(sqrt(wid))*500+500, res=300)
    print(taxa.info$tax.rarefy$sample$simpson)
    dev.off()
    png(file=paste(name, "tax.sample.S.obs", "png", sep="."), width=ceiling(sqrt(wid))*500+500, height=ceiling(sqrt(wid))*500+500, res=300)
    print(taxa.info$tax.rarefy$sample$S.obs)
    dev.off()
    png(file=paste(name, "tax.sample.S.chao1", "png", sep="."), width=ceiling(sqrt(wid))*500+500, height=ceiling(sqrt(wid))*500+500, res=300)
    print(taxa.info$tax.rarefy$sample$S.chao1)
    dev.off()
    png(file=paste(name, "tax.sample.S.ACE", "png", sep="."), width=ceiling(sqrt(wid))*500+500, height=ceiling(sqrt(wid))*500+500, res=300)
    print(taxa.info$tax.rarefy$sample$S.ACE)
    dev.off()
  }
  
  #####   Group based raw heatmap   #####
  #	taxa.info[["raw.heatmap.ori"]] <- my.plot_heatmap(physeq, ps.ord=NULL, method=NULL, label=label, col.type=c(samtype, spetype), color=c(col.sam, col.spe), 
  #	                                                  heatscale=heatscale, keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, normalize=seq.depth, trans=NULL)
  #	taxa.info[[paste("raw.heatmap.",gsub("[^a-zA-Z0-9]","",trans$name),sep="")]] <- my.plot_heatmap(physeq, ps.ord=NULL, method=NULL, label=label, col.type=c(samtype, spetype), color=c(col.sam, col.spe), 
  #	                                                                                                heatscale=heatscale, keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, normalize=seq.depth, trans=trans)
  
  #	png(file=paste(name, "tax.rawheatmap.ori", "png", sep="."), width=wid*50+500, height=hei*50+200, res=300)
  #	    print(taxa.info$raw.heatmap.ori)
  #	dev.off()
  #	png(file=paste(name, "tax.rawheatmap", gsub("[^a-zA-Z0-9]","",trans$name), "png", sep="."), width=wid*50+500, height=hei*50+200, res=300)
  #	    print(taxa.info[[paste("raw.heatmap.",gsub("[^a-zA-Z0-9]","",trans$name),sep="")]])
  #	dev.off()
  
  #####   Correlation based heatmap   #####
  
  
  
  #######    Sample Distance Analysis   #######
  #######    Sample Distance Matrix    #######
  dist.all <- list()
  if(sum(getslots.phyloseq(physeq)%in%"phy_tree")){
    # registerDoParallel(cores=7)
    # dist.all[["unifrac.w.non_nor"]][["distmatrix"]] <- distance(physeq, method="unifrac", weighted=TRUE, normalized=FALSE, parallel=TRUE, fast=TRUE)
    # dist.all[["unifrac.w.nor"]][["distmatrix"]] <- distance(physeq, method="unifrac", weighted=TRUE, normalized=TRUE, parallel=TRUE, fast=TRUE)
    # dist.all[["unifrac.non_w"]][["distmatrix"]] <- distance(physeq, method="unifrac", weighted=FALSE, parallel=TRUE, fast=TRUE)
    # dist.all[["dpcoa"]][["distmatrix"]] <- distance(physeq, "dpcoa")
    # dist.all[["jsd"]][["distmatrix"]] <- distance(physeq, "jsd", parallel=TRUE)
    phy_unifrac <- fastUnifrac(physeq, seq.depth=seq.depth, parallel=FALSE, method=c("Edge","dPCoA","non_w","w.non_nor","w.nor"))
    dist.all[["unifrac.w.non_nor"]][["distmatrix"]] <- phy_unifrac$dist$w.non_nor
    dist.all[["unifrac.w.nor"]][["distmatrix"]] <- phy_unifrac$dist$w.nor
    dist.all[["unifrac.non_w"]][["distmatrix"]] <- phy_unifrac$dist$non_w
    dist.all[["dpcoa"]][["distmatrix"]] <- phy_unifrac$dist$dPCoA
  }
  for(i in 1:length(distlist)){	dist.all[[distlist[i]]][["distmatrix"]] <- distance(physeq, method=distlist[i])	}
  
  #######   Sample Distance Plot   #######
  beta.test <- list()
  for(i in 1:length(dist.all)){
    distname <- names(dist.all)[i]
    write.csv(as.matrix(dist.all[[distname]][["distmatrix"]]), paste(name,"dist",distname,"csv",sep="."))
    Dis_ana <- distance_analysis(physeq, distmatrix=dist.all[[distname]][["distmatrix"]], axes=axes, samtype=samtype, color=col.sam, shape=shape, 
                                 label=label, keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, trans=trans, max.dist=max.dist, 
                                 normalize=seq.depth, heatscale=heatscale, flip=flip, td.axes=td.axes, distname=distname, name=name)
    dist.all[[distname]][["dist.analysis"]] <- Dis_ana
    #####   Output Result   #####
    beta.test[[distname]] <- Dis_ana$beta.test
    if(class(Dis_ana$Signif)!="try-error"){
      write.csv(Dis_ana$Signif$p.value, paste(name, "dist", distname, "pvalue", "csv", sep="."))
    }
    png(file=paste(name, "dist", distname, "boxplot","png", sep="."), width=3200, height=2000, res=300)
    print(Dis_ana$boxplot)
    dev.off()
    png(file=paste(name, "dist", distname, "barplot","png", sep="."), width=3200, height=2000, res=300)
    print(Dis_ana$barplot)
    dev.off()
    png(file=paste(name, "dist", distname, "pcoaplot","png", sep="."), width=3200, height=3000, res=300)
    print(Dis_ana$pcoaplot)
    dev.off()
    png(file=paste(name, "dist", distname, "nmdsplot","png", sep="."), width=3200, height=3000, res=300)
    print(Dis_ana$nmdsplot)
    dev.off()
    if(length(names(Dis_ana$igraph))){
      for(j in 1:length(names(Dis_ana$igraph))){
        queue_name <- names(Dis_ana$igraph[j])
        png(file=paste(name, "dist", distname, queue_name, "png", sep="."), width=3200, height=3000, res=300)
        print(Dis_ana$igraph[[j]])
        dev.off()
      }
    }
  }
  
  write.csv(do.call(cbind, beta.test), paste(name, "dist", "beta.test", "csv", sep="."))
  
  #######   Ordination Analysis   #######
  ordinate.all <- list()
  #######   DPCoA: vegan scores function have bug for DPCoA now   #######
  # if(sum(getslots.phyloseq(physeq)%in%"phy_tree")){	ordinate.all[["DPCoA"]][["ordination"]] <- ordinate(physeq, method="DPCoA")	}
  #######   Correspondence Analysis   #######
  ordinate.all[["UCA"]][["ordination"]] <- ordinate(physeq, method="CCA")
  ordinate.all[["CCA"]][["ordination"]] <- ordinate(as.formula(paste("physeq", samtype, sep="~")), method="CCA")
  ordinate.all[["DCA"]][["ordination"]] <- ordinate(physeq, method="DCA")
  #######   Redundancy Analysis   #######
  ordinate.all[["PCA"]][["ordination"]] <- ordinate(physeq, method="RDA")
  ordinate.all[["RDA"]][["ordination"]] <- ordinate(as.formula(paste("physeq", samtype, sep="~")), method="RDA")
  #######   NMDS Analysis   #######
  for(i in 1:length(distlist)){	ordinate.all[[paste("NMDS",distlist[i],sep=".")]][["ordination"]] <- ordinate(physeq, method="NMDS", distance=distlist[i])	}
  
  #######   Ordination Plot   #######
  for(i in 1:length(ordinate.all)){	
    ordname <- names(ordinate.all)[i]
    Ord_ana <- ordinate_analysis(physeq, ordination=ordinate.all[[ordname]][["ordination"]], col.type=c(samtype, spetype), color=c(col.sam, col.spe), 
                                 shape=shape, label=label, normalize=seq.depth, heatscale=heatscale, trans=trans, axes=axes, td.axes=td.axes, 
                                 keepOnlyTheseTaxa=keepOnlyTheseTaxa, keepOnlyTheseSample=keepOnlyTheseSample, ordname=ordname, name=name)
    ordinate.all[[ordname]][["ordi.analysis"]] <- Ord_ana
    #####   Output Result   #####
    png(file=paste(name, "ordi", ordname, "splitplot", "png", sep="."), width=6400, height=3000, res=300)
    print(Ord_ana$splitplot)
    dev.off()
    png(file=paste(name, "ordi", ordname, "biplot", "png", sep="."), width=3400, height=3000, res=300)
    print(Ord_ana$biplot)
    dev.off()
    #		png(file=paste(name, "ordi", ordname, "supplot", "png", sep="."), width=3400, height=3000, res=300)
    #			print(Ord_ana$supplot+par.axis+scale_color_manual(values=c(col.sam, col.spe)))
    #		dev.off()
    #####   Define scale for heatmap  #####
    #		hei <- nspecies(physeq)
    #		wid <- nsamples(physeq)
    #		if(!is.null(keepOnlyTheseTaxa)){	hei <- sum(keepOnlyTheseTaxa%in%species.names(physeq))	}
    #		if(!is.null(keepOnlyTheseSample)){	hei <- sum(keepOnlyTheseSample%in%sample.names(physeq))	}
    #		png(file=paste(name, "ordi", ordname, "heatmap.ori", "png", sep="."), width=wid*50+500, height=hei*50+200, res=300)
    #			print(Ord_ana$heatmap.ori)
    #		dev.off()
    #		png(file=paste(name, "ordi", ordname, "heatmap", gsub("[^a-zA-Z0-9]","",trans$name), "png", sep="."), width=wid*50+500, height=hei*50+200, res=300)
    #			print(Ord_ana[[paste("heatmap",gsub("[^a-zA-Z0-9]","",trans$name),sep=".")]])
    #		dev.off()
  }
  
  #######   Return Result   #######
  result <- list()
  result[["beta.test"]] <- beta.test
  result[["taxonomy.info"]] <- taxa.info
  result[["distance.info"]] <- dist.all
  result[["ordination.info"]] <- ordinate.all
  return(result)
  
  #	SA_basic(data, group, ID, name, uni=TRUE, unidist, sleep)
}
