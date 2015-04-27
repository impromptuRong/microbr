######################################################
#' Apply basic statistical test to whole matrix with proper alpha correction

#######   Para T-test and Wilcoxon-test   ######
para_t.test <- function(Data, group, p.adj=c("none","bonferroni"), ...){
  result <- list()
  
  mean_info <- summaryBy(.~group, data=data.frame(Data, group), FUN=mean)
  rownames(mean_info) <- paste(mean_info$group, "mean", sep=".")
  colnames(mean_info) <- c("group", colnames(Data))
  sd_info <- summaryBy(.~group, data=data.frame(Data, group), FUN=sd)
  rownames(sd_info) <- paste(sd_info$group, "sd", sep=".")
  colnames(sd_info) <- c("group", colnames(Data))
  result$summary <- t(rbind(mean_info, sd_info)[,-1])
  
  ngroup <- nlevels(group)
  names <- combn(levels(group),2)
  names <- paste(names[1,], names[2,], sep="vs")
  comb_info <- result$summary
  for(k in 1:length(p.adj)){
    method <- p.adj[k]
    pvalue <- matrix(0, nrow=ncol(Data), ncol=length(names), dimnames=list(colnames(Data), names))
    ###  T-test  ###
    plist <- lapply(1:ncol(Data), function(x){try(pairwise.t.test(Data[,x], group, pool.sd=FALSE, p.adj=method, ...)$p.value)})
    for(i in 1:length(plist)){
      if(class(plist[[i]])=="matrix"){
        pvalue[i,] <- unlist(sapply(1:(ngroup-1), function(x){plist[[i]][x:(ngroup-1),x]}))
      } else {
        pvalue[i,] <- NA
      }
    }
    result[["Ttest"]][[method]] <- pvalue
    col_name <- c(colnames(comb_info), paste("Ttest", method, sep="."), colnames(pvalue))
    comb_info <- data.frame(comb_info, insert=rownames(pvalue), pvalue)
    colnames(comb_info) <- col_name
    ###  Wilcoxon-Test  ###
    plist <- lapply(1:ncol(Data), function(x){try(pairwise.wilcox.test(Data[,x], group, p.adj=method, ...)$p.value)})
    for(i in 1:length(plist)){
      if(class(plist[[i]])=="matrix"){
        pvalue[i,] <- unlist(sapply(1:(ngroup-1), function(x){plist[[i]][x:(ngroup-1),x]}))
      } else {
        pvalue[i,] <- NA
      }
    }
    result[["Wilcox"]][[method]] <- pvalue
    col_name <- c(colnames(comb_info), paste("Wilcox", method, sep="."), colnames(pvalue))
    comb_info <- data.frame(comb_info, insert=rownames(pvalue), pvalue)
    colnames(comb_info) <- col_name
  }
  result$comb_info <- comb_info
  return(result)
}


#######   Para-Chisquare-test   ######
para_chi.test <- function(Data, group, paired=FALSE, alternative="two.sided", conf.level=0.95){
  col <- combn(levels(group), 2)
  ngroup <- length(levels(group))
  if(length(paired)==1){
    paired <- rep(paired, ncol(col))
  }
  summary <- matrix(0, ncol(Data), ngroup, dimnames=list(colnames(Data), paste(levels(group), "count", sep=".")))
  chitest <- matrix(0, ncol(Data), ncol(col), dimnames=list(colnames(Data), paste(col[1,], col[2,], sep="vs")))
  for(i in 1:ncol(col)){
    Cdata1 <- Data[group==col[1,i],]
    Cdata2 <- Data[group==col[2,i],]
    for(j in 1:ncol(Data)){
      r1 <- Cdata1[,j]
      r2 <- Cdata2[,j]
      a <- length(r1[r1!=0])
      b <- length(r2[r2!=0])
      c <- length(r1)-a
      d <- length(r2)-b
      if(i < ngroup){
        summary[j,1] <- paste(a, length(r1),sep="|")
        summary[j,i+1] <- paste(b, length(r2),sep="|")
      }
      if(paired[i]==TRUE){
        chitest[j,i] <- mcnemar.test(matrix(c(a,b,c,d),2,2),correct=TRUE)$p.value
      }
      else{
        chitest[j,i] <- fisher.test(matrix(c(a,b,c,d),2,2),alternative=alternative)$p.value
      }
    }
  }
  
  result <- list()
  result$summary <- data.frame(summary)
  result$chitest <- chitest
  result$comb_info <- data.frame(summary, chitest)
  return(result)
}

#######   Para-Kolmogorov-test   ######
para_ks.test <- function(Data, group, mu=0, paired=FALSE, alternative="two.sided", conf.level=0.95){
}

#######   Para-Mann's NP multiple range test   #########
para_Mann <- function(alldata, group, conf.level=0.95){
  require(coin)
  cat <- levels(group)
  grc <- combn(cat,2)
  result <- matrix(0, ncol(alldata), ncol(grc)+2)
  rownames(result) <- colnames(alldata)
  colnames(result) <- c("kruskal","NDWD",paste(grc[1,],grc[2,],sep=" - "))
  for(i in 1:ncol(alldata)){
    subdata <- data.frame(feat=alldata[,i],site=group)
    kw <- kruskal_test(feat~group, data=subdata, distribution=approximate(B=9999))
    result[i,1] <- pvalue(kw)[1]
    if(require("multcomp")){
      NDWD <- oneway_test(feat~group, data=subdata,
                          ytrafo=function(data) trafo(data, numeric_trafo=rank),
                          xtrafo=function(data) trafo(data, factor_trafo=function(x) model.matrix(~x-1) %*% t(contrMat(table(x),"Tukey"))),
                          teststat="max", distribution=approximate(B=90000))
      ### global p-value
      result[i,2] <- pvalue(NDWD)[1]
      result[i,3:ncol(result)] <- as.vector(pvalue(NDWD, method = "single-step"))
    }
  }
  detach("package:coin")
  return(result)
}

####   Multiple Test Correction   ####
multi_correct <- function(data, ...){
  rowname <- rownames(data)
  colname <- colnames(data)
  return(matrix(p.adjust(c(as.matrix(data)), ...), dim(data), dimnames=list(rowname, colname)))    
}
