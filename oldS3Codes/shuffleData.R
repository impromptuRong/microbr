#################################################
#' OTU table random shuffling for equally sequencing depth
#' May cause bias for Richness estimation

################   Shuffle Data   ###############
shuffleData <- function(data, depth, permu=1000, taxa_AreRows=TRUE, rm.zero=TRUE){
  plotdata <- data.matrix(data)
  if(!taxa_AreRows){	plotdata <- t(plotdata)	}
  Bacname <- rownames(plotdata)
  
  result <- list()
  for(i in 1:ncol(plotdata)){
    weight <- as.vector(plotdata[,i])
    ss <- table(as.factor(plotdata[,i]))
    if(sum(weight)>depth){
      ss <- replicate(permu, sample(rep(Bacname[weight>0],weight[weight>0]),depth))
      #			ss <- replicate(permu, sample(Bacname[weight>0],depth,prob=weight)
      ss <- factor(as.vector(ss),levels=Bacname)
      result[[i]] <- table(ss)
    }
    else{ result[[i]] <- plotdata[,i]*permu	}
  }
  
  result <- do.call(cbind, result)/permu
  colnames(result) <- colnames(plotdata)
  if(rm.zero){	result <- result[rowSums(result)>0,]	}
  return(result)
}
