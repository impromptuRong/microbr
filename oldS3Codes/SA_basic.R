#########################################################
#' This is an old version of basic statistical analysis based on otu table and group
#' most of element have been merged into SA_uni.R
#' To do list: merging
#' hcluster plot
#' Lin's CCC heatmap
#' Discrimitive Analysis
#' 
SA_basic <- function(data, group, ID, name, uni=FALSE, unidist=0, sleep=20){
  bcdist <- vegdist(data,"bray")
  write.csv(as.matrix(bcdist),paste(name,".dist.bc.csv",sep=""))
  
  ##################   PCA, (D)PCoA, COA, NMDS   ###################
  result <- list()
  test <- matrix(0,5,2)
  colnames(test) <- c("Bray-Curtis","Unifrac")
  rownames(test) <- c("Anosim","MRPP","perManova","NMDS2","NMDS3")
  
  library(ade4)
  library(rgl)
  pca <- dudi.pca(data, scannf=FALSE, nf=3)
  coa <- dudi.coa(data, scannf=FALSE, nf=3)
  nmds2.bc <- metaMDS(bcdist, k=2)
  nmds3.bc <- metaMDS(bcdist, k=3)
  
  result[[1]] <- pca$eig[1:3]
  test[4,1] <- nmds2.bc$stress
  test[5,1] <- nmds3.bc$stress
  
  hclusterplot(bcdist, group, paste(name,".hcluster.bc.png", sep=""))
  nmdsplot(nmds2.bc, group, paste(name,".nmds2.bc.png",sep=""), cex=1, label=ID)
  nmdsplot(nmds3.bc, group, paste(name,".nmds3.bc.png",sep=""), cex=1, label=ID)
  Sys.sleep(sleep)
  
  if(uni==TRUE){
    nmds2.uni <- metaMDS(unidist, k=2)
    nmds3.uni <- metaMDS(unidist, k=3)
    test[4,2] <- nmds2.uni$stress
    test[5,2] <- nmds3.uni$stress
    hclusterplot(unidist, group, paste(name,".hcluster.uni.png",sep=""))
    nmdsplot(nmds2.uni, group, paste(name,".nmds2.uni.png",sep=""),cex=1, label=ID)
    nmdsplot(nmds3.uni, group, paste(name,".nmds3.uni.png",sep=""),cex=1, label=ID)
    Sys.sleep(sleep)
  }
  
  pcplot(pca, group, paste(name,".pca.png",sep=""), cbox=0.5, label=ID)
  Sys.sleep(sleep)
  pcplot(coa, group, paste(name,".coa.png",sep=""), cbox=0.5, label=ID)
  Sys.sleep(sleep)
  
  #########################################################
  con <- file(paste(name,".SA.log",sep=""), open="wt")
  sink(file=con, type="output")
  
  ########   Multi-test    #########
  #######  ANOSIM, MRPP, perMANOVA  #######
  anosim.bc <- anosim(bcdist, group, permutations=1000)
  mrpp.bc <- mrpp(bcdist, group, permutations=1000)
  per.bc <- adonis(bcdist~group, permutations=1000)
  #		perMANOVA.bc <- adonis(data ~ group, method="bray", permutations=1000)
  test[1,1] <- anosim.bc$signif
  test[2,1] <- mrpp.bc$Pvalue
  test[3,1] <- per.bc$aov.tab$"Pr(>F)"[1]
  
  print(anosim.bc)
  print(mrpp.bc)
  print(per.bc)
  
  if(uni==TRUE){
    anosim.uni <- anosim(unidist, group, permutations=1000)
    mrpp.uni <- mrpp(unidist, group, permutations=1000)
    per.uni <- adonis(unidist~group, permutations=1000)
    
    test[1,2] <- anosim.uni$signif
    test[2,2] <- mrpp.uni$Pvalue
    test[3,2] <- per.uni$aov.tab$"Pr(>F)"[1]
    
    print(anosim.uni)
    print(mrpp.uni)
    print(per.uni)
  }
  result[[2]] <- test
  
  #########################################################
  ##################  Diversity Index  ####################
  library(BiodiversityR)
  shannon <- diversity(data, index="shannon", MARGIN=1, base=exp(1))
  simpson <- diversity(data, index="simpson", MARGIN=1, base=exp(1))
  #		k <- sample(nrow(data), 9)
  #		rendiv <- renyi(data[k,])
  #		plot(rendiv)
  
  richness<- estimateR(round(data))
  #		Srar <- rarefy(Data_C, ceiling(min(rowSums(Data_C))))
  #		sac <- specaccum(Data)
  #		plot(sac, ci.type="polygon", ci.col="yellow")
  
  eco <- data.frame(t(rbind(shannon, simpson, richness)))
  eco$group <- group
  write.csv(eco, paste(name,".diver.csv",sep=""))
  
  eco_ttest <- para_t.test(eco[,c(1,2,3,4,6)], eco$group, paired=FALSE, alternative="two.sided")
  eco_npttest <- para_npt.test(eco[,c(1,2,3,4,6)], eco$group, paired=FALSE, alternative="two.sided")
  write.csv(t(cbind(eco_ttest,eco_npttest)), paste(name,".diver.test.csv",sep=""))
  
  #########################################################
  ##################         DA       #####################
  library(caret)
  library(lattice)
  library(rgl)
  library(MASS)
  ###########    LDA    ##########
  lda <- lda(group~., data, CV=FALSE)
  lda.hat <- predict(lda, decision.values=TRUE)
  tabTrain <- confusionMatrix(group, lda.hat$class)
  denplot(lda.hat$x, group, paste(name,".lda.density.png",sep=""))
  if(ncol(lda.hat$x)>1){
    ldaplot(lda.hat$x, group, paste(name,".lda.xyplot.png",sep=""))
  }
  Sys.sleep(sleep)
  
  ###########    PLSDA   #########
  useBayes <- plsda(data, group, ncomp=3, probMethod="Bayes")
  useSoftmax <- plsda(data, group, ncomp=3)
  useBayes.hat <- predict(useBayes, data)
  useSoftmax.hat <- predict(useSoftmax, data)
  print(confusionMatrix(predict(useBayes, data), group))
  print(confusionMatrix(predict(useSoftmax, data), group))
  
  #########################################################
  ##################  Lin's CCC heatmap  ##################
  #	library(ClassDiscovery)
  #	library(ClassComparison)
  #	library(epiR)
  #	library(proxy)
  #		lin <- function(x1, x2, ci="z-transform", conf.level=0.95){
  #			tmp <- epi.ccc(x1, x2, ci, conf.level)
  #			return(tmp$rho.c[,1])
  #		}
  #		pr_DB$set_entry(FUN=lin, names="Lin")
  #		method="Lin"
  #		heatmaplot(data, group, ID, method=method, paste(name,".",method,".heatmap.png",sep=""))
  #	detach("package:ClassDiscovery")
  #	detach("package:ClassComparison")
  #	detach("package:epiR")
  #	detach("package:proxy")
  
  sink()
  close(con)
  
  print(c("feature"=ncol(data),"sample"=nrow(data)))
  print(table(group))
  print(result[[1]])
  print(result[[2]])
}