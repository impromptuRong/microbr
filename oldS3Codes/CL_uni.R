##############################################
#' Classification and Confusion Matrix for different method after feature selection
#' CL_uni is a wapper for CL_basic
CL_uni <- function(physeq, cat, name, kernel="linear", plot=TRUE){
  data <- data.frame(otu_table(physeq))
  sample_data <- data.frame(sample_data(physeq))
  # 	tree <- phy_tree(physeq)
  taxonomy <- data.frame(tax_table(physeq))
  
  ID <- sample.names(physeq)
  feature <- species.names(physeq)
  group <- sample_data[,cat]
  result <- CL_basic(data, group, ID, name, kernel=kernel, plot=plot)
  return(result)
}

CL_basic <- function(data, group, ID, name, kernel="linear", plot=TRUE){
  result <- rep(0,4)
  
  con <- file(paste(name,".CL.log",sep=""), open="wt")
  sink(file=con, type="output")
  
  library(caret)
  ##########  PDA  ##########
  if(ncol(data)>2){
    ###########    PLSDA   #########
    useBayes <- plsda(data, group, ncomp=3, probMethod="Bayes")
    useSoftmax <- plsda(data, group, ncomp=3)
    useBayes.hat <- predict(useBayes, data)
    useSoftmax.hat <- predict(useSoftmax, data)
    cof1 <- confusionMatrix(predict(useBayes, data), group)
    cof2 <- confusionMatrix(predict(useSoftmax, data), group)
    print(cof1)
    print(cof2)
    result[1] <- cof1$overall[1]*100
    result[2] <- cof2$overall[1]*100
    detach("package:klaR")
    detach("package:pls")
    detach("package:MASS")
  }
  ##########  SVM  ##########
  library(e1071)
  svm <- svm(group~., data, cross=nrow(data), kernel=kernel)
  svm.hat <- predict(svm,decision.values=TRUE)
  tabTrain <- confusionMatrix(group, svm.hat)
  cof3 <- summary(svm)
  print(cof3)
  result[3] <- cof3$tot.accuracy
  detach("package:e1071")
  ##########  RF  ##########
  library(randomForest)
  count <- 20
  while(count > 0){
    trf <- randomForest(group~., data=data, proximity=TRUE)
    acc <- 100-trf$err.rate[nrow(trf$err.rate),1]*100
    if(result[4] < acc){
      result[4] <- acc;
      count <- 20;
      rf <- trf;
    }
    count <- count-1;
  }
  print(rf)
  if(plot==TRUE){
    rf.cmd <- cmdscale(1-rf$proximity)
    gn <- nlevels(group)
    color <- c("blue","red","green","purple")[1:gn]
    png(file=paste(name,".rf.mds.png",sep=""), width=1000, height=1000, res=300)
    par(mar=c(1.5,1.5,1,1)+0.1, font.axis=1, mgp=c(2,0.15,0), cex=0.8, tck=-0.005)
    plot(rf.cmd, pch=21, col=color[as.numeric(group)], bg=color[as.numeric(group)], ann=FALSE)
    setoff <- max(rf.cmd[,2])/30-min(rf.cmd[,2])/30
    s.label(cbind(rf.cmd[,1], rf.cmd[,2]+setoff), lab=ID, clabel=0.5, add.plot=TRUE, boxes=F, grid=F)
    legend("topright",legend=levels(group),col=color,pt.bg=color,pch=21,cex=1)
    dev.off()
  }
  detach("package:randomForest")
  ################   Logestic Regression ##############
  if(ncol(data)<nrow(data)){
    glm.r <- glm(group~.,cbind(group,data), family=binomial(logit))
    print(summary(glm.r))
    print(predict(glm.r,type="response"))
  }
  detach("package:caret")
  
  sink()
  close(con)
  
  return(result)
}
