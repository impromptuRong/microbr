#################################################
#' Old version of FS_basic:
#' Methods not in new pipeline:
#' SPLSDA
#' Multi-group SVM-RFE
FS_basic_old <- function(data, group, ID, name, method=c("StaTest", "ISA", "RF", "ENet", "SVMRFE", "SPLSDA"), kernel="linear"){
  method <- c("StaTest", "ISA", "RF", "ENet", "SVMRFE", "SPLSDA") %in% method;
  if(!sum(method)){
    warning("Input Methods is not correct! Only show Univariate Test!")
    method = c(1,0,0,0,0,0)
  }
  feature <- sort(colnames(data))
  data <- data[,feature]
  result <- list()
  
  ########################################################
  ################  basic test selection  ################
  if(method[1]){
    library(coin)
    parattest <- para_t.test(data, group, p.adj="none", paired=FALSE, alternative="two.sided")$comb_info
    parachitest <- para_chi.test(data, group, paired=FALSE,alternative="two.sided")
    result1 <- cbind(parattest, FisherExact=feature, parachitest)
    result[["BTest"]] <- result1
    write.csv(result1, paste(name, ".fs.Test.pairwise.csv",sep=""))
    
    if(length(levels(group))>2){
      aov <- sapply(c(1:ncol(data)),function(x){summary(aov(data[,x]~group))[[1]][["Pr(>F)"]][1]})
      # fisher <- sapply(c(1:ncol(data)),function(x){fisher.test(table(data[,x]==0,group),alternative=TRUE)$p.value})
      fisher <- sapply(c(1:ncol(data)),function(x){
        t<-table(data[,x]==0,group);
        p<-1;
        if(nrow(t)>1){p<-fisher.test(t,alternative=TRUE)$p.value}
        return(p);})
      MGC <- sapply(c(1:ncol(data)),function(x){TukeyHSD(aov(data[,x]~group))$group[,4]})
      colnames(MGC) <- colnames(data)
      np_MGC <- para_Mann(data,group,conf.level=0.95)
      comb_MGT <- t(rbind(fisher, aov, MGC, t(np_MGC)))
      result[["MTest"]] <- comb_MGT
      write.csv(comb_MGT, paste(name,".fs.Test.group.csv",sep=""))
    }
    detach("package:coin")
  }
  
  ################  Indicator Species Analysis  #############
  if(method[2]){
    library(labdsv)
    isa <- indval(data, group, numitr=1000)
    #		summary(isa)
    summary(isa, p=0.05, digits=2, show=p)
    result2 <- data.frame(isa$indcls,isa$maxcls,isa$pval)
    result2 <- result2[feature,]
    isa.core <- result2[isa$pval<=0.05,]
    write.csv(isa.core, paste(name,".fs.isa.csv",sep=""))
    detach("package:labdsv")
    result[["ISA"]] <- result2
  }
  
  ################  Random Forest Selection ############
  if(method[3]){
    library(randomForest)
    rf <- randomForest(data, group, ntree=1000)
    result3 <- importance(rf)
    result[["RF"]] <- result3
    write.csv(result3, paste(name, ".fs.rf.csv",sep=""))
  }
  
  #    ################  Random Forest Boruta  ##############
  #    if(method[3]){
  #        library(randomForest)
  #        library(Boruta)
  #        set.seed(125)
  #        boruta <- Boruta(group~., data=cbind(group,data), doTrace=2, maxRuns=12, ntree=500) #  Default maxRuns=4
  #        borutaplot(boruta,paste(name, ".fs.rf.png",sep=""))
  #        tmp <- TentativeRoughFix(boruta)
  #        result3 <- cbind(boruta$finalDecision,tmp$finalDecision)
  #        result3 <- result3[feature,]
  #        write.csv(data.frame(boruta$finalDecision,tmp$finalDecision), paste(name, ".fs.rf.csv",sep=""))
  #        getConfirmedFormula(clean.Boruta)
  #        getNonRejectedFormula(clean.Boruta)
  #        detach("package:Boruta")
  #        detach("package:randomForest")
  #        result[["RF"]] <- result3
  #    }
  
  ################   ELASTIC NET   ####################
  if(method[4]){
    library(glmnet)
    elastic_mod <- paraelastic(data, group, family="multinomial")
    ENfsplot(elastic_mod,paste(name, ".fs.EN.png",sep=""))
    result4 <- as.matrix(elastic_mod[[5]])[feature,]
    write.csv(as.matrix(elastic_mod[[5]]),paste(name,".fs.EN.coef.csv",sep=""))
    write.csv(elastic_mod[[6]],paste(name,".fs.EN.par.csv",sep=""))
    detach("package:glmnet")
    result[["ENet"]] <- result4
  }
  
  ################  SVM-RFE selection  ################
  if(method[5]){
    library(e1071)
    library(kernlab)
    ftable <- para_svmrfe(data, group, ID, kernel="linear")
    result5 <- ftable
    write.csv(ftable,paste(name, ".fs.svmrfe.",kernel,".csv",sep=""))
    detach("package:e1071")
    detach("package:kernlab")
    result[["SVMRFE"]] <- result5
  }
  
  ################      HHSVM      ####################
  
  ################   SPLSDA  Selection  ################
  if(method[6]){
    library(spls)
    library(nnet)
    library(MASS)
    cv.lda <- cv.splsda(data.matrix(data), as.numeric(group), classifier="lda", K=c(1:5), eta=seq(0.1,0.9,0.1), scale.x=FALSE, fold=5, n.core=4)
    cv.log <- cv.splsda(data.matrix(data), as.numeric(group), classifier="logistic", K=c(1:5), eta=seq(0.1,0.9,0.1), scale.x=FALSE, fold=5, n.core=4)
    
    splsda.lda <- splsda(data.matrix(data), as.numeric(group), classifier="lda", eta=cv.lda$eta.opt, K=cv.lda$K.opt, scale.x=FALSE)
    splsda.log <- splsda(data.matrix(data), as.numeric(group), classifier="logistic", eta=cv.log$eta.opt, K=cv.log$K.opt, scale.x=FALSE)
    print(splsda.lda)
    print(splsda.log)
    
    par <- matrix(c(splsda.lda$eta,splsda.lda$K,splsda.lda$kappa,splsda.log$eta,splsda.log$K,splsda.lda$kappa),3,2)
    colnames(par) <- c(splsda.lda$classifier, splsda.log$classifier)
    rownames(par) <- c("eta","K","kappa")
    write.csv(splsda.lda$W, paste(name, ".fs.splsda.lda.csv", sep=""))
    write.csv(splsda.log$W, paste(name, ".fs.splsda.log.csv", sep=""))
    write.csv(par, paste(name, ".fs.splsda.par.csv", sep=""))
    detach("package:spls")
    detach("package:nnet")
    detach("package:MASS")
    result[["SPLSDA"]][["lda"]] <- splsda.lda$w
    result[["SPLSDA"]][["log"]] <- splsda.log$w
    result[["SPLSDA"]][["par"]] <- par
  }
  
  ################    Combine Result Score  #################
  if(sum(method[1:5])==5){
    score <- matrix(0, length(feature), 6)
    rownames(score) <- feature
    colnames(score) <- c("ENET","ISA","RF","SVMRFE","Test","Score")
    score[,5] <- rowSums(cbind(as.numeric(result1[,3]),as.numeric(result1[,6]),as.numeric(result1[,7]))<=0.05)
    score[score>=2] <- 1	
    score[,2] <- result2[,3]<=0.05
    score[names(result5),4] <- result5
    #		result3[result3==3] <- 0
    #		score[,3] <- rowSums(result3/2)
    score[,3] <- result3[,4]/max(result3[,4])*2
    score[,1] <- rowSums(abs(result4)>0)/2
    score[,6] <- score[,1]/2+score[,2]/0.5+score[,3]+score[,4]/20+score[,5]/0.5
    write.csv(score, paste(name, ".fs.final.csv", sep=""))
    result[["score"]][["mat"]] <- score
    
    ##################   Tune Score cutoff   ##################
    cre <- seq(0,max(score[,6]),by=0.1)
    flist <- lapply(cre, function(x){feature[score[,6]>=x]})
    
    accuracy <- sapply(flist,function(x){CL_basic(data[x], group, ID, "all", kernel="linear", plot=FALSE)})
    unlink("all.CL.log")
    rownames(accuracy) <- c("PLSDAlda", "PLSDAbayes", "SVM", "RandomForest")
    tune <- t(rbind(cre, accuracy))
    write.csv(tune, paste(name, ".fs.tune.csv", sep=""))
    result[["score"]][["tune"]] <- tune
    
    png(file=paste(name,".fs.tune.png",sep=""), width=2000, height=1000, res=300)
    par(mar=c(2.5,2.5,2,1)+0.1, lwd=0.5, font.axis=1, cex.axis=1, pch=21, cex=1, mgp=c(2,0,0), tck=-0.01)
    plot(cre, accuracy[2,],main="",ylim=c(0,100),col="green",bg="green",ann=FALSE)
    #			points(cre, accuracy[2,],col="cyan",bg="cyan")
    points(cre, accuracy[3,],col="red",bg="red")
    points(cre, accuracy[4,],col="blue",bg="blue")
    mtext(side=1, text="score", font=2, line=1)
    mtext(side=2, text="Accuracy", font=2, line=1)
    legend("bottomright",legend=c("PLSDA","SVM","RandomForest"),col=c("green","red","blue"),pt.bg=c("green","red","blue"),pch=21,cex=1)
    dev.off()
    ##################################
  }
  
  ################   Return part of the model   #############
  #	return(elastic_mod)
  return(result)
}
