######################################################
#' Enterotype analysis based on its official website
###########   Enterotype_analysis   ##############
enterotype_analysis <- function(phylo, name){
  require(ade4)
  require(cluster)
  require(clusterSim)
  
  if(class(phylo)=="phyloseq"){
    otu_table <- otu_table(phylo)@.Data
    if(otu_table(phylo)@taxa_are_rows){otu_table <- t(otu_table)}
  } else {
    otu_table <- as.matrix(phylo)
  }
  data <- otu_table    
  
  ##  Dissimilarity function and PAM clustering
  data.dist = dist.JSD(data)
  data.cluster = pam.clustering(data.dist, k=3)
  
  ##  Find Optimal number of clusters
  nclusters = NULL
  for (k in 1:20) {
    if (k==1){
      nclusters[k]=NA
    } else {
      data.cluster_temp=pam.clustering(data.dist, k)
      nclusters[k]=index.G1(t(data), data.cluster_temp,  d = data.dist, centrotypes = "medoids")
    }
  }
  
  png(file=paste(name, "kcluster", "png", sep="."), width=4000, height=2000, res=300)
  plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
  dev.off()
  kclust = which.max(nclusters)
  
  ##  This has shown that the optimal number of clusters for this particular dataset is 3 (k=3).
  data.cluster=pam.clustering(data.dist, k=kclust)
  ##  Cluster validation 
  obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
  
  #####  Graphical interpretation  #####
  ## Between-class analysis
  ## Remove noise
  # data.denoized = noise.removal(data, percent=0.01)
  
  ## PCA plot
  obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
  obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1)
  png(file=paste(name, "bca", "png", sep="."), width=3000, height=3000, res=300)
  s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)
  dev.off()
  
  ## PCoA plot
  obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
  png(file=paste(name, "pcoa", "png", sep="."), width=3000, height=3000, res=300)
  s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F)
  dev.off()
  
  ## Return Result
  result <- list()
  result[["k"]] <- kclust
  result[["silhouette"]] <- obs.silhouette
  result[["cluster"]] <- data.cluster
  
  return(result)
}

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) {
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]), as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}
##################################################