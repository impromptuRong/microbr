currpath <- "/Users/ruichenrong/Projects/DongPipeline/microbr"
library(phyloseq)
library(glmnet)
library(gcdnet)
library(penalizedSVM)
source("./R_code/z.source.r")

##############    Phylum Input   ###############
currpath <- "/Users/ruichenrong/Projects/DongPipeline/microbr"
setwd(currpath)
## Meta Data ##
metainfo <- read.csv("./0.raw/meta.info.csv", header=TRUE, row.names=1)
## Phylogenetic Tree ##
rawtree <- ape::read.tree("./0.raw/oral_phylum.ref.tre")
## Taxonomy File ##
taxonomy <- read.csv("./0.raw/dic.phylum.csv", header=TRUE, row.names=1)
## Abundance Table ##
rawdata <- read.csv("./0.raw/oral.summary.phylum.unrounded.csv", header=TRUE, row.names=1)

oral_raw
oral <- physet(otu_table=rawdata, sample_data=metainfo, tax_table=taxonomy, phy_tree=rawtree)

oral <- physet(otu_table=oral_raw$rawdata, sample_data=oral_raw$metainfo, tax_table=oral_raw$taxonomy)



## Phyloseq Class and Unifrac ##
raw.physeq <- phyloseq(otu_table(raw.data, taxa_are_rows=FALSE), phy_tree(raw.tree), 
                       tax_table(as.matrix(taxonomy)), sample_data(metainfo))
unilist <- fastUnifrac(raw.physeq, seq.depth=1000)
## unifs_cv1 Model ##
unifs_cv.w <- read.csv("./0.matlab.out/oral_phylum.unifs_cv1.weight.csv", header=TRUE, row.names=1)[colnames(unilist$edge_bool),]

alpha <- seq(0,1,0.1)
glmnet_cv_ori <- lapply(alpha, function(x){cv.glmnet(unilist$edge_matrix, metainfo$Periodontitis, nfolds=10, alpha=x, family="binomial", type.measure="class", standardize=FALSE)})
glmnet_cv_wei <- lapply(alpha, function(x){cv.glmnet(unilist$edge_matrix, metainfo$Periodontitis, weights=c(raw.tree$edge.length,0), nfolds=10, alpha=x, family="binomial", type.measure="class", standardize=FALSE)})
coef(glmnet_cv_wei[[11]], s="lambda.min")

################################################

########  Visualization Feature Tree  ##########
#####    Label color based for Taxonomy    #####
lab.type <- factor(tax_table(raw.physeq)[,"Phylum"])
n <- nlevels(lab.type)
hues <- rainbow(n+1)
hues <- seq(15, 375, length=n+1)
col.lab <- c(hcl(h=hues, l=65, c=100)[1:n],"#000000")
names(col.lab) <- c(levels(lab.type), "unclassified")
col.label <- col.lab[lab.type]
barplot(rep(1,n+1), yaxt="n", col=col.lab)

#####   Branch color and width for weight  #####
weight <- apply(abs(unifs_cv.w), 2, function(x){x/max(x)})
weight.adj <- log(1000*weight+1)/3
edge.width <- weight*2

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
barplot(rep(1,100), yaxt="n", col=jet.colors(100))
barplot(rep(1,30), yaxt="n", col=topo.colors(30))
par(mfrow=c(3,4))
for(k in 1:ncol(weight)){
  edge.color <- rev(topo.colors(30))[as.numeric(cut(weight[,k], breaks=30))]
  plot(phy_tree(raw.physeq), edge.color=edge.color, edge.width=edge.width, tip.color=col.label, font=1, cex=1, direction="downwards")
}
barplot(rep(1,30), yaxt="n", col=topo.colors(30))
################################################

########     NMDS Plot Comparison     ##########
all1 <- unilist$dist$dPCoA
all2 <- unilist$dist$w.non_nor
# c1 <- as.matrix(UniFrac(raw.physeq, weighted=FALSE, normalized=FALSE))
# c2 <- as.matrix(UniFrac(raw.physeq, weighted=TRUE, normalized=FALSE))
# c3 <- as.matrix(UniFrac(raw.physeq, weighted=TRUE, normalized=TRUE))
# c4 <- as.matrix(DPCoA(raw.physeq))
edgematrix <- unilist$edge_matrix
edgelength <- c(phy_tree(raw.physeq)$edge.length, 0)
for(k in 0:ncol(weight)){
  coref <- rep(TRUE, ncol(edgematrix))
  if(k){
    coref <- weight[,k]>10^-4
  }
  D1 <- matrix(0, nrow(edgematrix), nrow(edgematrix), dimnames=list(rownames(edgematrix), rownames(edgematrix)))
  D2 <- D1
  for(i in 1:nrow(edgematrix)){
    for(j in 1:nrow(edgematrix)){
      D1[i,j] <- (edgematrix[i,coref]-edgematrix[j,coref])^2 %*% edgelength[coref]
      D2[i,j] <- abs(edgematrix[i,coref]-edgematrix[j,coref]) %*% edgelength[coref]
    }
  }
  nmds1 <- plot_ordination(raw.physeq, ordinate(raw.physeq, method="NMDS", distance=as.dist(D1)), color="Group", shape="Periodontitis") + 
    geom_point(size=5) + scale_color_manual(values=c("blue","red","green","purple")) + theme_bw() + 
    scale_x_continuous(name="NMDS1") + scale_y_continuous(name="NMDS2")
  nmds2 <- plot_ordination(raw.physeq, ordinate(raw.physeq, method="NMDS", distance=as.dist(D2)), color="Group", shape="Periodontitis") + 
    geom_point(size=5) + scale_color_manual(values=c("blue","red","green","purple")) + theme_bw() + 
    scale_x_continuous(name="NMDS1") + scale_y_continuous(name="NMDS2")
  print(nmds1)
  print(nmds2)
}


##############################################################################
# tt <- cv.gcdnet(try.data, group, nfolds=10, method="hhsvm", pred.loss="misclass", standardize=FALSE, eps=1e-8, delta=0.1)
# gcdnet(try.data, group, nfolds=10, method="hhsvm", pred.loss="misclass", lambda2=0.5, standardize=FALSE, eps=1e-8, delta=0.1)
# ttt <- gcdnet(try.data, group, nlambda=100, method="hhsvm", lambda2=tt$lambda.min, standardize=FALSE, eps=1e-8, delta=0.1)
##############################################################################

od <- c("Acidobacteria","SR1","Firmicutes","Tenericutes","Fusobacteria","Cyanobacteria_Chloroplast","Chloroflexi",
        "Spirochaetes","Deinococcus_Thermus","Planctomycetes","Actinobacteria","Gemmatimonadetes","Synergistetes",
        "TM7","Bacteroidetes","Proteobacteria","Branch_1","Branch_2","Branch_3","Branch_4","Branch_5","Branch_6",
        "Branch_7","Branch_8","Branch_9","Branch_14","Branch_13","Branch_12","Branch_11","Branch_10","Root")

d = 0.008521655191070;
D = 1;
N = 31;
lambda2 = 0.002766614411988/D/N*d;
lambda1 = 0.002766614411988/D/N;

bw <- c(phy_tree(raw.physeq)$edge.length,0.00001)
data <- unilist$edge_matrix*100
group <- sam_data(raw.physeq)$Periodontitis
group <- c(1, -1)[as.numeric(group)]

result <- gcdnet(data, group, nlambda=1, method="hhsvm", lambda=lambda1, lambda2=lambda2, pf=1/bw, pf2=1/bw, standardize=FALSE, delta=0.0001)
as.matrix(result$beta[od,])
write.csv(as.matrix(result$beta[od,]), "pwd.csv")

gcdnet(x, y, nlambda = 100,
       method = c("hhsvm", "logit", "sqsvm", "ls"),
       lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
       lambda = NULL, lambda2 = 0,
       pf = rep(1, nvars), pf2 = rep(1, nvars), exclude,
       dfmax = nvars + 1, pmax = min(dfmax * 1.2, nvars), standardize = TRUE, eps = 1e-8, maxit = 1e6, delta = 2)

result <- gcdnet(data, group, method="hhsvm", lambda2=0, pf=1/bw, pf2=1/bw, standardize=FALSE, delta=0.0001)

