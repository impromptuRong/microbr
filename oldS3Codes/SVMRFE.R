######################################################
#' Implement SVM-RFE method for binary/multiple group feature selection
#' Support linear, RBF, ploynomial kernel for now
#' 2/3 straitified training: RFE + gridsearch for best par
#' evaluation: training CVacc + 1/3 testing acc

para_svmrfe <- function(data, group, ID, kernel="linear"){
  feature <- list()
  print(kernel)
  for(n in 1:100){
    traindata <- data.frame()
    testdata <- data.frame()
    cat <- levels(group)
    for(i in 1:length(cat)){
      c_data <- data[group==cat[i],]
      c_train <- sample(1:nrow(c_data),5)
      c_test <- setdiff(1:nrow(c_data),c_train)
      traindata <- rbind(traindata, data.frame(group=rep(cat[i],5),c_data[c_train,]))
      testdata <- rbind(testdata, data.frame(group=rep(cat[i],length(c_test)),c_data[c_test,]))
    }
    g_train <- traindata[,1]
    d_train <- traindata[,-1]
    d_train <- d_train[,colSums(d_train)!=0]
    traindata <- data.frame(group=g_train,d_train)
    testdata <- testdata[,colnames(traindata)]
    levels(testdata[,1]) <- levels(traindata[,1])
    
    svmrfe <- svmRFE(group~., traindata, testdata, kernel=kernel)
    #		svmRFEplot(svmrfe,filename="svmrfe.1.png")
    #		write.csv(svmrfe$featureRank,"svmrfe.1.csv")
    frank <- svmrfe$featureRank
    #		index <- sort((frank[,3]+frank[,4]),index.return=TRUE)$ix
    index <- sort((frank[,2]+frank[,4]),index.return=TRUE)$ix
    #		index <- sort((frank[,2]+frank[,3]+frank[,4]),index.return=TRUE)$ix
    
    feature[[n]] <- rownames(frank)[1:index[length(index)]]
  }
  return(table(unlist(feature)))
}

#######   svmRFE feature selection  #######
svmRFE <- function(formula, train, test, kernel="linear"){
  fstr <- toString(formula)
  type <- strsplit(fstr,"\\, ")[[1]][2]
  
  featureRankedList <- matrix(0,(ncol(train)-1),7)
  rownames(featureRankedList) <- c(1:(ncol(train)-1))
  colnames(featureRankedList) <- c("F_rank","Classification_Acc","Train_Acc","Test_Acc","Opt_Cost","Opt_Gamma","nSV")
  featuresort <- rep(0,(ncol(train)-1))
  
  rankedFeatureIndex <- ncol(train)-1
  survivingFeatures <- colnames(train)
  survivingFeatures <- c(survivingFeatures[survivingFeatures!=type],type)
  result <- list()
  
  while(length(survivingFeatures)>1){
    # Initialize Data
    train_data <- train[,survivingFeatures]
    test_data <- test[,survivingFeatures]
    if(kernel == "radial"){
      # Grid Search find a coarse Grid
      obj <- gridsearch(formula, train=train_data, test=test_data, cross=nrow(train_data), 
                        ranges=list(gamma=2^(-15:3),cost=2^(-5:15)),kernel="radial")
      coarse_Gamma <- log2(obj$best.model$gamma)
      coarse_C <- log2(obj$best.model$cost)
      range_Gamma <- seq(coarse_Gamma-2,coarse_Gamma+2,by=0.25)
      range_C <- seq(coarse_C-2,coarse_C+2,by=0.25)
      
      # Grid Search find best parameters
      obj_loo <- gridsearch(formula, train=train_data, test=test_data, cross=nrow(train_data), 
                            ranges=list(gamma=2^range_Gamma,cost=2^range_C), kernel="radial")
      
      Gamma_loo <- obj_loo$best.model$gamma
      C_loo <- obj_loo$best.model$cost
      print(obj_loo$best.performance)
      svmModel_loo <- obj_loo$best.model
      svmModel_loo <- svm(formula, data=train_data, cross=nrow(train_data), type="C-classification", 
                          kernel=kernel, cost=C_loo, gamma=Gamma_loo)
    }
    
    if(kernel == "linear"){
      svmModel_loo <- svm(formula, train_data, cross=nrow(train_data), kernel=kernel)
      #cachesize=500, 
    }
    
    # compute ranking criteria
    rankingCriteria <- svmweights(svmModel_loo)
    ranking <- sort(rankingCriteria, index.return=TRUE)$ix
    #		result$model[[rankedFeatureIndex]] <- svmModel_loo
    #		result$rankC[[rankedFeatureIndex]] <- rankingCriteria
    
    svmModel_loo_hat <- predict(svmModel_loo, train_data, decision.values=TRUE)
    svmModel_test_hat <- predict(svmModel_loo, test_data, decision.values=TRUE)
    tabTrain_loo <- confusion(train_data[,type], svmModel_loo_hat, printit=FALSE)
    tabTrain_test <- confusion(test_data[,type], svmModel_test_hat, printit=FALSE)
    
    # update feature ranked list
    featureRankedList[rankedFeatureIndex,1] <- rankedFeatureIndex
    featureRankedList[rankedFeatureIndex,2] <- svmModel_loo$tot.acc/100
    #		featureRankedList[rankedFeatureIndex,2] <- obj_loo$best.performance
    featureRankedList[rankedFeatureIndex,3] <- tabTrain_loo$acc
    featureRankedList[rankedFeatureIndex,4] <- tabTrain_test$acc
    featureRankedList[rankedFeatureIndex,5] <- svmModel_loo$cost
    featureRankedList[rankedFeatureIndex,6] <- svmModel_loo$gamma
    featureRankedList[rankedFeatureIndex,7] <- sum(svmModel_loo$nSV)
    featuresort[rankedFeatureIndex] <- survivingFeatures[ranking[1]]
    
    rankedFeatureIndex <- rankedFeatureIndex-1
    # eliminate the feature with smallest ranking criterion
    survivingFeatures <- survivingFeatures[-ranking[1]]
  }
  rownames(featureRankedList) <- featuresort
  result$featureRank <- featureRankedList
  
  return(result)
}

################################################
# weights and Criteria of the hiperplane
################################################
svmweights <- function(model){
  rankingCriteria <- rep(0,ncol(model$SV))
  
  #	rbf <- function(u, v, gamma){
  #		exp(-gamma*sum((u-v)^2))
  #	}
  #	class(rbf) <- "kernel"
  #	rbf <- rbfdot(sigma=model$gamma)
  
  if(model$nclasses==2){
    if(model$kernel==0){
      w <- t(model$coefs) %*% model$SV
      rankingCriteria <- w * w
    }
    if(model$kernel==1){
    }
    if(model$kernel==2){
      for(f in 1:ncol(model$SV)){
        KMat <- (model$coefs %*% t(model$coefs)) * kernelMatrix(rbfdot(sigma=model$gamma),model$SV[,-f],model$SV[,-f])
        rankingCriteria[f] <- -sum(KMat)
      }
    }
  }
  
  else{
    start <- c(1, cumsum(model$nSV)+1)
    start <- start[-length(start)]
    
    W <- matrix(0,ncol(model$SV),choose(model$nclasses,2))
    count <- 1
    for(i in 1:(model$nclasses-1)){
      for(j in (i+1):model$nclasses){
        ## ranges for class i and j:
        ri <- start[i]:(start[i] + model$nSV[i] - 1)
        rj <- start[j]:(start[j] + model$nSV[j] - 1)
        ## coefs and SV for (i,j):
        coefs <- c(model$coefs[ri, j-1], model$coefs[rj, i])
        SV <- data.matrix(model$SV[c(ri,rj),])
        if(model$kernel==0){
          w <- t(coefs) %*% SV
          W[,count] <- w * w
        }
        if(model$kernel==1){
        }
        if(model$kernel==2){
          for(nf in 1:ncol(model$SV)){
            KMat <- (coefs %*% t(coefs)) * kernelMatrix(rbfdot(sigma=model$gamma),SV[,-nf],SV[,-nf])
            W[nf,count] <- -sum(KMat)
          }
        }
        count <- count+1
      }
    }
    rankingCriteria <- rowMeans(W)
  }
  return(rankingCriteria)
}

gridsearch <- function(formula, train, test, cross=10, kernel="radial", ranges=list(gamma=2^(-15:3),cost=2^(-5:15))){
  fstr <- toString(formula)
  type <- strsplit(fstr,"\\, ")[[1]][2]
  
  para_grid <- expand.grid(ranges)
  para_list <- lapply(seq_len(nrow(para_grid)), function(i){para_grid[i,]})
  
  para_svm <- lapply(para_list, function(x){
    svm(formula, data=train, cross=cross, gamma=x$gamma, cost=x$cost, 
        cachesize=500, type="C-classification", kernel="radial")})
  
  para_train_predict <- lapply(para_svm, function(x){predict(x, train, decision.values=TRUE)})
  para_train_acc <- lapply(para_train_predict, function(x){confusion(train[,type], x, printit=FALSE)})
  
  para_test_predict <- lapply(para_svm, function(x){predict(x, test, decision.values=TRUE)})
  para_test_acc <- lapply(para_test_predict, function(x){confusion(test[,type], x, printit=FALSE)})
  
  performance <- data.frame(gamma=unlist(lapply(para_svm,function(x){x$gamma})),
                            cost=unlist(lapply(para_svm,function(x){x$cost})),
                            TrainError=unlist(lapply(para_svm,function(x){1-mean(x$acc)/100})),
                            TrainDispersion=unlist(lapply(para_svm,function(x){sd(x$acc)/100})),
                            TrainPredict=unlist(lapply(para_train_acc,function(x){1-x$acc})),
                            TestPredict=unlist(lapply(para_test_acc,function(x){1-x$acc}))
  )
  
  performance$SumError <- performance$TrainPredict+performance$TestPredict
  #	rank <- order(performance$SumError, performance$cost, performance$gamma)
  rank <- order(performance$TrainError, performance$cost, performance$gamma)
  result <- list(best.model=para_svm[[rank[1]]], best.performance=performance[rank[1],], performance=performance)
  return(result)
}

###### A function that calculates the confusion matrix and overall accuracy #####
confusion <- function(actual, predicted, names=NULL, printit=TRUE, prior=NULL){
  if(is.null(names)){	names <- levels(actual)	}
  result <- list()
  tab <- table(actual, predicted)
  acctab <- t(apply(tab, 1, function(x)x/sum(x)))
  dimnames(acctab) <- list(Actual=names,"Predicted (cv)"=names)
  result$tab <- acctab
  if(is.null(prior)){
    relnum <- table(actual)
    prior <- relnum/sum(relnum)
    acc <- sum(tab[row(tab)==col(tab)])/sum(tab)
    result$acc <- acc
  }
  else{
    acc <- sum(prior*diag(acctab))
    names(prior) <- names
    result$acc <- acc
  }
  if(printit){
    print(round(c("Overall accuracy"=acc,"Prior frequency"=prior),4))
  }
  if(printit){
    cat("\nConfusion matrix","\n")
    print(round(acctab,4))
  }
  return(result)
}
