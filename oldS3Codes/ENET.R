##########################################################
#' Elastic Net wrapper for feature selection and tuning parameters
#' binomial/multinomial logistic regression
#' gaussian/poisson regression
#' cox-survival model 
#' 
paraelastic <- function(data, group, family="multinomial"){
  alpha <- seq(0,1,0.1)
  if(family=="multinomial"|family=="binomial"){
    paracvfit1 <- lapply(alpha, function(x){cv.glmnet(data.matrix(data),group,nfolds=5,alpha=x,family=family,type.measure="class")})
  }
  if(family=="gaussian"|family=="poisson"){
    paracvfit1 <- lapply(alpha, function(x){cv.glmnet(data.matrix(data),group,nfolds=5,alpha=x,family="gaussian",type.measure="deviance")})
  }
  paracvfit2 <- lapply(alpha, function(x){cv.glmnet(data.matrix(data),group,nfolds=5,alpha=x,family=family,type.measure="mse")})
  paracvfit1[[1]] <- NULL
  paracvfit2[[1]] <- NULL
  
  cvm1 <- matrix(unlist(lapply(paracvfit1,function(x){x$cvm})),nrow=(length(alpha)-1),byrow=TRUE)
  cvup1 <- matrix(unlist(lapply(paracvfit1,function(x){x$cvup})),nrow=(length(alpha)-1),byrow=TRUE)
  #cvsd1 <- matrix(unlist(lapply(paracvfit1,function(x){x$cvsd})),nrow=(length(alpha)-1),byrow=TRUE)
  #nzero1 <- matrix(unlist(lapply(paracvfit1,function(x){x$nzero})),nrow=(length(alpha)-1),byrow=TRUE)
  lamda1 <- matrix(unlist(lapply(paracvfit1,function(x){x$lambda})),nrow=(length(alpha)-1),byrow=TRUE)
  
  cvm2 <- matrix(unlist(lapply(paracvfit2,function(x){x$cvm})),nrow=(length(alpha)-1),byrow=TRUE)
  cvup2 <- matrix(unlist(lapply(paracvfit2,function(x){x$cvup})),nrow=(length(alpha)-1),byrow=TRUE)
  #cvsd2 <- matrix(unlist(lapply(paracvfit2,function(x){x$cvsd})),nrow=(length(alpha)-1),byrow=TRUE)
  #nzero2 <- matrix(unlist(lapply(paracvfit2,function(x){x$nzero})),nrow=(length(alpha)-1),byrow=TRUE)
  lamda2 <- matrix(unlist(lapply(paracvfit2,function(x){x$lambda})),nrow=(length(alpha)-1),byrow=TRUE)
  
  result1.cvm <- which(cvm1==min(cvm1), arr.ind=TRUE)
  result1.cvup<- which(cvup1==min(cvup1), arr.ind=TRUE)
  result2.cvm <- which(cvm2==min(cvm2), arr.ind=TRUE)
  result2.cvup<- which(cvup2==min(cvup2), arr.ind=TRUE)
  
  lam1.min.cvm <- unlist(lapply(c(1:nrow(result1.cvm)),function(x){lamda1[result1.cvm[x,1],result1.cvm[x,2]]}))
  lam1.index.cvm <- which(lam1.min.cvm==max(lam1.min.cvm), arr.ind=TRUE)
  mod1.cvm <- paracvfit1[[result1.cvm[lam1.index.cvm,1]]]
  
  lam1.min.cvup <- unlist(lapply(c(1:nrow(result1.cvup)),function(x){lamda1[result1.cvup[x,1],result1.cvup[x,2]]}))
  lam1.index.cvup <- which(lam1.min.cvup==max(lam1.min.cvup), arr.ind=TRUE)
  mod1.cvup <- paracvfit1[[result1.cvup[lam1.index.cvup,1]]]
  
  lam2.min.cvm <- unlist(lapply(c(1:nrow(result2.cvm)),function(x){lamda2[result2.cvm[x,1],result2.cvm[x,2]]}))
  lam2.index.cvm <- which(lam2.min.cvm==max(lam2.min.cvm), arr.ind=TRUE)
  mod2.cvm <- paracvfit2[[result2.cvm[lam2.index.cvm,1]]]
  
  lam2.min.cvup <- unlist(lapply(c(1:nrow(result2.cvup)),function(x){lamda2[result2.cvup[x,1],result2.cvup[x,2]]}))
  lam2.index.cvup <- which(lam2.min.cvup==max(lam2.min.cvup), arr.ind=TRUE)
  mod2.cvup <- paracvfit2[[result2.cvup[lam2.index.cvup,1]]]
  
  mod.par <- matrix(c(result1.cvm[lam1.index.cvm,1]*0.01,mod1.cvm$lambda.1se,mod1.cvm$lambda.min,
                      result1.cvup[lam1.index.cvup,1]*0.01,mod1.cvup$lambda.1se,mod1.cvup$lambda.min,
                      result2.cvm[lam2.index.cvm,1]*0.01,mod2.cvm$lambda.1se,mod2.cvm$lambda.min,
                      result2.cvup[lam2.index.cvup,1]*0.01,mod2.cvup$lambda.1se,mod2.cvup$lambda.min),
                    4,3,byrow=TRUE)
  colnames(mod.par) <- c("alpha","lambda.1se","lambda.min")
  
  if(family=="multinomial"|family=="binomial"){
    rownames(mod.par) <- c("cla.cvm","cla.cvup","mse.cvm","mse.cvup")
    
    mod1.f.cvm.1se  <- do.call("cBind",coef(mod1.cvm, s="lambda.1se"))
    mod1.f.cvm.min  <- do.call("cBind",coef(mod1.cvm, s="lambda.min"))
    mod1.f.cvup.1se <- do.call("cBind",coef(mod1.cvup,s="lambda.1se"))
    mod1.f.cvup.min <- do.call("cBind",coef(mod1.cvup,s="lambda.min"))
    mod2.f.cvm.1se  <- do.call("cBind",coef(mod2.cvm, s="lambda.1se"))
    mod2.f.cvm.min  <- do.call("cBind",coef(mod2.cvm, s="lambda.min"))
    mod2.f.cvup.1se <- do.call("cBind",coef(mod2.cvup,s="lambda.1se"))
    mod2.f.cvup.min <- do.call("cBind",coef(mod2.cvup,s="lambda.min"))
    
    mod.coef <- cBind(mod1.f.cvm.1se,mod1.f.cvm.min,mod1.f.cvup.1se,mod1.f.cvup.min,
                      mod2.f.cvm.1se,mod2.f.cvm.min,mod2.f.cvup.1se,mod2.f.cvup.min)
    colnames(mod.coef) <- rep(levels(group),8)
  }
  
  if(family=="gaussian"|family=="poisson"){
    rownames(mod.par) <- c("dev.cvm","dev.cvup","mse.cvm","mse.cvup")
    
    mod.coef <- cBind(coef(mod1.cvm,s="lambda.1se"),coef(mod1.cvm,s="lambda.min"),coef(mod1.cvup,s="lambda.1se"),coef(mod1.cvup,s="lambda.min"),
                      coef(mod2.cvm,s="lambda.1se"),coef(mod2.cvm,s="lambda.min"),coef(mod2.cvup,s="lambda.1se"),coef(mod2.cvup,s="lambda.min"))
    colnames(mod.coef) <- c("dev.cvm.1se","dev.cvm.min","dev.cvup.1se","dev.cvup.min","mse.cvm.1se","mse.cvm.min","mse.cvup.1se","mse.cvup.min")
  }
  
  result <- list()
  result[[1]] <- mod1.cvm
  result[[2]] <- mod1.cvup
  result[[3]] <- mod2.cvm
  result[[4]] <- mod2.cvup
  result[[5]] <- mod.coef
  result[[6]] <- mod.par
  return(result)
}