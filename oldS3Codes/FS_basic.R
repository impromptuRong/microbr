#################################################
#' Basic supervised learning feature selection methods
#' Given otu-table and grouping/respond information
#' Statistical Test: T-test, Fisher exact test, Np T-test, Multi-group comparison
#' Random Forest Selection: Standard entropy method, Boruta selection
#' ENET Selection: ENET Logistic Regression, DrSVM, HHSVM
#' SVM-RFE: 2-group, multi-group
#' Indicator Species: labdsv
#' para_t.test
#' para_chi.test
#' para_Mann

FS_basic <- function(data, group, ID, name, method=c("StaTest", "RF", "ENet", "ISA", "SVMRFE"), kernel="linear"){
  require(doBy)
  require(ggplot2)
  method <- c("StaTest", "RF", "ENet", "ISA", "SVMRFE") %in% method
  if(!sum(method)){
    warning("Input Methods is not correct! Only show Univariate Test!")
    method = c(1,0,0,0,0)
  }
  feature <- sort(colnames(data))
  data <- data[, feature]
  result <- list()
  summary <- list()
  
  ########################################################
  ################  Basic test selection  ################
  if(method[1]){
    print("#### 1. Basic test selection ####")
    print("     1.1 Two-Group Comparison")
    parattest <- para_t.test(data, group, p.adj="none", paired=FALSE, alternative="two.sided")
    parachitest <- para_chi.test(data, group, paired=FALSE, alternative="two.sided")
    result[["BTest.pvalue"]] <- data.frame(parattest$Ttest$none, parattest$Wilcox$none, parachitest$chitest)
    colnames(result[["BTest.pvalue"]]) <- c("Ttest", "Wilcox", "FisherExact")
    result[["BTest.comb_info"]] <- cbind(parattest$comb_info, FisherExact=colnames(data), parachitest$comb_info)
    write.csv(result[["BTest.comb_info"]], paste(name, ".fs.Test.pairwise.csv", sep=""))
    print("#################################")
    summary[["BTest"]] <- result[["BTest.pvalue"]]<=0.05
    
    if(length(levels(group))>2){
      print("#################################")
      print("     1.2 Multi-Group Comparison")
      aov <- sapply(c(1:ncol(data)), function(x){summary(aov(data[,x]~group))[[1]][["Pr(>F)"]][1]})
      # fisher <- sapply(c(1:ncol(data)), function(x){fisher.test(table(factor(data[,x]==0,levels=c(TRUE,FALSE)), group),alternative="two.sided")$p.value})
      fisher <- sapply(c(1:ncol(data)), function(x){fisher.test(factor(data[,x]==0,levels=c(TRUE,FALSE)),group,alternative="two.sided")$p.value})
      MGC <- sapply(c(1:ncol(data)),function(x){TukeyHSD(aov(data[,x]~group))$group[,4]})
      colnames(MGC) <- colnames(data)
      np_MGC <- para_Mann(data, group, conf.level=0.95)
      result[["MTest"]] <- t(rbind(fisher, aov, MGC, t(np_MGC)))
      write.csv(result[["MTest"]], paste(name, ".fs.Test.group.csv", sep=""))
      print("#################################")
      # summary[["MTest"]] <- result[["BTest.pvalue"]]<=0.05
    }
  }
  
  ########################################################
  ##############   RandomForest Selection   ##############
  if(method[2]){
    print("##  2. RandomForest selection  ##")
    #######      Standard Method       #######
    print("    2.1 Standard Method")
    require(party)
    set.seed(326)
    data.cforest <- cforest(group ~ ., data=data.frame(group=group, data), controls=cforest_unbiased(ntree=1000, mtry=3))
    if(nlevels(group)==2){
      data.cforest.varimp <- varimpAUC(data.cforest, conditional=TRUE)
    } else{
      data.cforest.varimp <- varimp(data.cforest, conditional=TRUE)
    }
    criteria <- abs(min(data.cforest.varimp))
    plotdata <- data.frame(Importance=data.cforest.varimp, Taxa=factor(reorder(names(data.cforest.varimp), data.cforest.varimp)))
    write.csv(plotdata, paste(name, ".fs.RF.coef.csv",sep=""))
    result[["RF.coef"]] <- matrix(plotdata$Importance,,1,dimnames=list(rownames(plotdata),"RF.coef"))
    
    # rf <- randomForest(data, group, ntree=1000)
    # result[["RF"]] <- importance(rf)
    
    plotdata <- plotdata[order(plotdata$Importance,decreasing=TRUE),]
    plotscale <- min(nrow(plotdata), max(sum(plotdata$Importance>criteria), 50))
    result[["RF.varImp"]] <- ggplot(plotdata, aes(x=Importance, y=Taxa)) + geom_point(color="blue",size=1.5) + geom_vline(aes(xintercept=0),colour="blue") + 
      geom_vline(aes_string(xintercept=criteria),colour="red",linetype="longdash") + geom_vline(aes_string(xintercept=-criteria),colour="red",linetype="longdash") + 
      scale_y_discrete(limits=plotdata$Taxa[plotscale:1]) + theme_bw()
    
    png(file=paste(name,".fs.RF.tune.png",sep=""), width=2000, height=plotscale*50+100, res=300)
    print(result[["RF.varImp"]])
    dev.off()
    print("#################################")
    summary[["RF.coef"]] <- result[["RF.coef"]]>=abs(min(result[["RF.coef"]]))
    
    #######       Boruta Method        #######
    print("    2.2 Boruta Method")
    require(Boruta)
    boruta <- Boruta(group~., data=data.frame(group=group, data), doTrace=2, maxRuns=12, ntree=500) #  Default maxRuns=4
    result[["RF.Boruta"]] <- data.frame(Boruta=boruta$finalDecision)
    write.csv(result[["RF.Boruta"]], paste(name, ".fs.RF.Boruta.csv",sep=""))
    # borutaplot(boruta, paste(name, ".fs.RF.Boruta.png",sep=""))
    # tmp <- TentativeRoughFix(boruta)
    detach("package:Boruta")
    detach("package:party")
    detach("package:randomForest")
    print("#################################")
    summary[["RF.Boruta"]] <- result[["RF.Boruta"]]=="Confirmed"
  }
  
  ########################################################
  #############    ENET LogesicRegression    #############
  ################      HHSVM/DrSVM      #################
  if(method[3]){
    print("###  3. ENET Based selection  ###")
    #######   ENET LogesicRegression   #######
    require(glmnet)        
    print("     3.1 ENET LogesticRegression")
    alpha <- seq(0,1,0.1)
    family <- ifelse(nlevels(group)>2, "multinomial", "binomial")
    measure <- ifelse(nlevels(group)==2&&length(group)>=100, "auc", "class")
    cv.enet <- lapply(alpha, function(x){cv.glmnet(data.matrix(data), group, alpha=x, standardize=FALSE, family=family, type.measure=measure)})
    result[["cv.enet.par"]] <- rbind(alpha, sapply(cv.enet, function(x){c(min(x$cvm),x$lambda.min,x$lambda.1se)}))
    rownames(result[["cv.enet.par"]]) <- c("alpha", cv.enet[[1]]$name, "lambda.1se", "lambda.min")
    
    if(nlevels(group)>2){
      result[["cv.enet.coef.1se"]] <- data.matrix(do.call("cBind", lapply(cv.enet, function(x){do.call("cBind", coef(x,s="lambda.1se"))})))
      result[["cv.enet.coef.min"]] <- data.matrix(do.call("cBind", lapply(cv.enet, function(x){do.call("cBind", coef(x,s="lambda.min"))})))
    } else{
      cv.coef <- data.matrix(do.call("cBind", lapply(cv.enet, function(x){coef(x,s="lambda.1se")})))
      result[["cv.enet.coef.1se"]] <- cbind(cv.coef,-cv.coef)[,c(t(matrix(1:(length(alpha)*2),,2)))]
      cv.coef <- data.matrix(do.call("cBind", lapply(cv.enet, function(x){coef(x,s="lambda.min")})))
      result[["cv.enet.coef.min"]] <- cbind(cv.coef,-cv.coef)[,c(t(matrix(1:(length(alpha)*2),,2)))]
    }
    cname <- t(expand.grid(levels(group), alpha))
    cname <- paste(cname[1,],cname[2,],sep="_")
    colnames(result[["cv.enet.coef.1se"]]) <- colnames(result[["cv.enet.coef.min"]]) <- cname
    
    write.csv(result[["cv.enet.coef.1se"]], paste(name,".fs.ENet.coef.1se.csv",sep=""))
    write.csv(result[["cv.enet.coef.min"]], paste(name,".fs.ENet.coef.min.csv",sep=""))
    write.csv(result[["cv.enet.par"]], paste(name,".fs.ENet.par.csv",sep=""))
    
    coef.min <- data.frame(result[["cv.enet.coef.min"]][-1,]!=0)
    colnames(coef.min) <- cname
    tmp <- reshape(coef.min, direction="long", varying=cname, sep="_")
    tmp$coef <- rowSums(tmp[,-c(1,ncol(tmp))])>0
    coef.min <- reshape(tmp[,c("time","id","coef")], direction="wide", new.row.names=rownames(coef.min))
    coef.1se <- data.frame(result[["cv.enet.coef.1se"]][-1,]!=0)
    colnames(coef.1se) <- cname
    tmp <- reshape(coef.1se, direction="long", varying=cname, sep="_")
    tmp$coef <- rowSums(tmp[,-c(1,ncol(tmp))])>0
    coef.1se <- reshape(tmp[,c("time","id","coef")], direction="wide", new.row.names=rownames(coef.1se))
    summary[["cv.enet.coef"]] <- cbind(coef.min[,-1], coef.1se[,-1])
    cname <- t(expand.grid(alpha, c("enet.min","enet.1se")))
    cname <- paste(cname[2,],cname[1,],sep="_")
    colnames(summary[["cv.enet.coef"]]) <- cname
    print("#################################")
    
    if(nlevels(group)==2){
      print("     3.2 ENET SVM")
      #######      HHSVM      #######
      require(gcdnet)
      lambda2 <- c(seq(0,1,0.2), seq(2,10,2))
      cv.hhsvm <- lapply(lambda2, function(x){cv.gcdnet(data.matrix(data), group, nfolds=10, pred.loss="misclass", method="hhsvm", standardize=FALSE, lambda2=x)})
      index <- which.min(sapply(cv.hhsvm, function(x){min(x$cvm)}))
      #######   SCAD+L2 SVM   #######
      require(penalizedSVM)
      cv.drHSVM <- try(svm.fs(data.matrix(data), y=as.numeric(group)*2-3, fs.method="scad+L2", cross.outer=0, inner.val.method="gacv", show="none", cross.inner=10, verbose=FALSE), TRUE)
      while(class(cv.drHSVM)=="try-error"){
        cv.drHSVM <- try(svm.fs(data.matrix(data[,sample(colnames(data))]), y=as.numeric(group)*2-3, fs.method="scad+L2", cross.outer=0, inner.val.method="gacv", show="none", cross.inner=10, verbose=FALSE), TRUE)
      }
      cv.coef <- c(data.matrix(coef(cv.hhsvm[[index]],s="lambda.min")), cv.drHSVM$model$b, cv.drHSVM$model$w[colnames(data)])
      cv.coef[is.na(cv.coef)] <- 0
      result[["cv.DrSVM.coef"]] <- matrix(cv.coef, ncol=2, dimnames=list(c("(Intercept)",colnames(data)),c("hhsvm","drHSVM")))
      result[["cv.DrSVM.par"]] <- matrix(c(min(cv.hhsvm[[index]]$cvm), cv.hhsvm[[index]]$lambda.min, lambda2[index], 
                                           sum(cv.drHSVM$classes!=as.numeric(group)*2-3)/length(group), cv.drHSVM$lambda1, cv.drHSVM$lambda2), 
                                         nrow=3, ncol=2, dimnames=list(c("Misclassification Error","lambda1","lambda2"),c("hhsvm","drHSVM")))
      write.csv(result[["cv.DrSVM.coef"]], paste(name,".fs.DrSVM.coef.csv",sep=""))
      write.csv(result[["cv.DrSVM.par"]], paste(name,".fs.DrSVM.par.csv",sep=""))
      print("#################################")
      summary[["cv.DrSVM.coef"]] <- result[["cv.DrSVM.coef"]][-1,]!=0
    }
  }
  
  ########################################################
  ###########    Indicator Species Analysis    ###########
  if(method[4]){
    print("####   4. IndicatorSpecies   ####")        
    require(labdsv)
    isa <- indval(data, group, numitr=1000)
    # summary(isa, p=0.05, digits=2, show=p)
    result[["ISA"]] <- cbind(isa$indval, pval=isa$pval)
    write.csv(result[["ISA"]], paste(name, ".fs.ISA.csv", sep=""))
    detach("package:labdsv")
    print("#################################")
    summary[["ISA"]] <- matrix(result[["ISA"]]$pval<=0.05, ncol(data), 1, dimnames=list(colnames(data),"ISA"))
  }
  
  ########################################################
  ################    SVMRFE selection    ################
  if(method[5]){
    print("####   5. SVMRFE Selection   ####")
    require(e1071)
    require(kernlab)
    if(nlevels(group)==2){
      require(pathClass)
      cv.svmrfe <- crossval(data, group, theta.fit=fit.rfe, folds=10, repeats=10, scale="scale", stepsize=0.1, parallel=TRUE, DEBUG=FALSE)
      cv.choose <- extractFeatures(cv.svmrfe, toFile=FALSE)[colnames(data),]
      cv.choose[is.na(cv.choose)] <- 0
      result[["svmrfe"]] <- matrix(cv.choose, ncol(data), 1, dimnames=list(colnames(data),"svmrfe"))
      write.csv(result[["svmrfe"]], paste(name,".fs.SVMRFE.linear.csv",sep=""))
      # cv.auc <- data.frame(cv.svmrfe$auc)
      # plotdata <- reshape(cv.auc, direction="long", varying=colnames(cv.auc), v.names="AUC", timevar="Repeat", times=colnames(cv.auc), idvar="Fold", ids=rownames(cv.auc))
      # plotdata$Repeat <- factor(reorder(plotdata$Repeat, 1:nrow(plotdata)))
      # plotdata$Fold <- factor(reorder(plotdata$Fold, 1:nrow(plotdata)))
      # ggplot(plotdata, aes(x=Repeat, y=AUC)) + geom_boxplot() + labs(y="AUC", x="") + scale_y_continuous(limits=c(0,1)) + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1))
      
      # pred <- prediction(cv.svmrfe$cv, matrix(as.numeric(group)*2-3, nrow(data), ncol(cv.svmrfe$cv)))
      # perf <- performance(pred, measure="tpr", x.measure="fpr")
      # auc <- unlist(performance(pred, "auc")@y.values)
      # a <- data.frame(fpr=unlist(perf@x.values), tpr=unlist(perf@y.values), rep=rep(colnames(cv.svmrfe$cv), each=length()))
      # plot(perf, avg="horizontal", spread.estimate="boxplot", main=paste("mean AUC = ", round(mean(auc),digits=4), sep=""))
      print("#################################")
      summary[["svmrfe"]] <- result[["svmrfe"]]>=50
    }
    
    # ftable <- para_svmrfe(data, group, ID, kernel="linear")
    # result5 <- ftable
    # write.csv(ftable,paste(name, ".fs.svmrfe.",kernel,".csv",sep=""))
  }
  
  ########################################################
  #############    Result Summary/Combine    #############
  result$summary <- summary
  return(result)
}