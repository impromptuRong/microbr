#######################################################
#' Regression Models for continuous respond case
#' Error estimation for regression model before/after feature selection
#' Things not included yet: possion regression, bayesian estimation
################   Basic Linear Regression   ###############
LR_uni <- function(physeq, cat, name, family=gaussian()){
  data <- data.frame(otu_table(physeq))
  sample_data <- data.frame(sample_data(physeq))
  tree <- phy_tree(physeq)
  taxonomy <- data.frame(tax_table(physeq))
  
  ID <- sample.names(physeq)
  feature <- species.names(physeq)
  RV <- as.vector(sample_data[,cat])
  LR_basic(data, RV, ID, name=name, family=family)
}

LR_basic <- function(data, RV, ID, name, family=gaussian()){
  fstr <- family$family
  linkstr <- family$link
  
  con <- file(paste(name,".LR.log",sep=""), open="wt")
  sink(file=con, type="output")
  
  #	##############    ELASTIC NET   #################
  #	if(fstr=="gaussian"|fstr=="bionomial"|fstr=="poisson"){
  #		library(glmnet)
  #		elastic_mod <- paraelastic(data, RV, family=fstr)
  #		ENfsplot(elastic_mod,paste(name, ".LR.EN.png",sep=""))
  #		write.csv(as.matrix(elastic_mod[[5]]),paste(name,".LR.EN.coef.csv",sep=""))
  #		write.csv(elastic_mod[[6]],paste(name,".LR.EN.par.csv",sep=""))
  #		detach("package:glmnet")
  #	}
  
  #############   AIC for Combined Model   ############
  if(ncol(data)<nrow(data)){
    library(MASS)
    m0 <- lm(RV~1, data=cbind(RV,data))
    m1 <- lm(RV~., data=cbind(RV,data))
    print(summary(m1))
    print(extractAIC(m1))
    coefs <- FScoefs(m0, m1, data=cbind(RV,data))
    AICplot(coefs, paste(name,".AIC.coef.png",sep=""))
    m0 <- glm(RV~1, data=cbind(RV,data), family=family)
    m1 <- glm(RV~., data=cbind(RV,data), family=family)
    print(summary(m1))
    result <- step(m0, scope=list(lower=formula(m0),upper=formula(m1)), direction="forward", trace=FALSE)
    print(result)
    result <- step(m1, scope=list(lower=formula(m0),upper=formula(m1)), direction="backward", trace=FALSE)
    print(result)
    detach("package:MASS")
  }
  
  ##############     Linear Models   ###############
  if(ncol(data)>=nrow(data)){
    #############   LM for each features   ############
    lm_sep <- matrix(0,ncol(data),4)
    colnames(lm_sep) <- c("abu","pvalue","R2","R2.adj")
    rownames(lm_sep) <- colnames(data)
    dir.create(paste("./LR_",name,sep=""))
    setwd(paste("./LR_",name,sep=""))
    order <- order(colSums(data),decreasing=TRUE)
    for(j in 1:length(order)){
      i <- order[j]
      bcname <- colnames(data)[i]
      abu <- sum(data[,i])/sum(data)*100
      mod <- lm(RV~.,data=cbind(RV,data)[,c(1,i+1)])
      p <- 1-pf(summary(mod)$fstatistic[1],summary(mod)$fstatistic[2],summary(mod)$fstatistic[3])
      R2 <- summary(mod)$r.squared
      R2.adj <- summary(mod)$adj.r.squared
      lm_sep[bcname,1] <- abu
      lm_sep[bcname,2] <- p
      lm_sep[bcname,3] <- R2
      lm_sep[bcname,4] <- R2.adj
      if(p <= 0.05){
        LRplot(mod, paste(j,round(abu,2),bcname,"lm.png",sep="_"))
      }
    }
    setwd("../")
    write.csv(lm_sep,paste(name,".LR.lm.sep.csv",sep=""))
  }
  
  sink()
  close(con)
}

################   AIC_FS_coefs for Regression   ###############
FScoefs <- function(m0, m1, data, trace=FALSE, direction="forward") {
  keepCoef <- function(m, aic){
    all <- names(coef(m1))
    new <- names(coef(m))
    ans <- rep(0, length(all))
    ans[match(new, all)] <- coef(m)
    ans
  }
  out <- with(data, stepAIC(m0, scope=list(lower=formula(m0), upper=formula(m1)), k=0, trace=trace, keep=keepCoef, direction=direction))
  rownames(out$keep) <- names(coef(m1))
  out$keep
}

