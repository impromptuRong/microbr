################################################################################
# What's in this wrapper-methods.R file: 
# alpha diversity, richness and rarefactions:  
# ordinate_analysis: 
# supervise_learning: 
# unsupervise_learning: 
# feature_selection: 
# enterotpe_analysis: 
################################################################################
#' @title Statistical Analysis for physet object
#' @description These methods are used to calculate univariate statistical 
#' analysis. Including T-test, Wilcoxon-test, Fisher-exact-test, Macnemar-test, 
#' Anova, Kruskal test, Tukey HSD and Mann's NP multiple group comparison. 
#' @param x A \code{matrix} or a \code{\link{physet-class}} objects.
#' @param formula A formula used for constrain analysis. 
#' @param data The data frame contains the grouping variable. 
#' @param ... other parameters. not used.
#' @return A large list with the following structure: 
#' 
#' @rdname statistics-method
#' @examples
#' data(oral)
#' phyTtest(oral, . ~ Group)
#' phyCtest(oral, . ~ Group)
#' phyMtest(oral, . ~ Group + Sex)
#' @export
phyTtest <- function(x, formula, data = sample_data(x), ...) {
  arglist <- c(list(x = quote(x), formula = formula, 
                    data = data), list(...))
  if (!"p.adj" %in% arglist)
    arglist$p.adj <- "none"
  if (inherits(x, "physet"))
    x <- as.matrix(otu_table(x))
  rhs <- all.vars(formula[-2])
  lhs <- all.vars(formula[-3])
  if (all(lhs != "."))
    x <- x[, lhs]
  group <- factor(data[rownames(x), rhs[1]])
  
  result <- list()
  mean_info <- doBy::summaryBy(. ~ group, data=data.frame(x, group), FUN=mean)
  rownames(mean_info) <- paste(mean_info$group, "mean", sep=".")
  colnames(mean_info) <- c("group", colnames(x))
  sd_info <- doBy::summaryBy(. ~ group, data = data.frame(x, group), FUN = sd)
  rownames(sd_info) <- paste(sd_info$group, "sd", sep=".")
  colnames(sd_info) <- c("group", colnames(x))
  result$summary <- t(rbind(mean_info, sd_info)[,-1])
  
  nG <- length(unique(group))
  names <- combn(unique(group), 2)
  names <- paste(names[1,], names[2,], sep = "vs")
  comb_info <- result$summary
  for (k in 1:length(arglist$p.adj)) {
    method <- arglist$p.adj[k]
    pvalue <- matrix(0, nrow = ncol(x), ncol = length(names), 
                     dimnames = list(colnames(x), names))
    
    ###  T-test  ###
    plist <- lapply(1:ncol(x), function(m)
      try(suppressWarnings(pairwise.t.test(x[, m], group, pool.sd=FALSE, 
                           p.adjust.method = method)$p.value)))
    for(i in 1:length(plist)){
      if(class(plist[[i]]) == "matrix"){
        pvalue[i,] <- unlist(sapply(1:(nG-1), function(m){plist[[i]][m:(nG-1), m]}))
      } else {
        pvalue[i,] <- NA
      }
    }
    result[["Ttest"]][[method]] <- pvalue
    col_name <- c(colnames(comb_info), paste("Ttest", method, sep="."), 
                  colnames(pvalue))
    comb_info <- data.frame(comb_info, insert = rownames(pvalue), pvalue)
    colnames(comb_info) <- col_name
    
    ###  Wilcoxon-Test  ###
    plist <- lapply(1:ncol(x), function(m) 
      try(suppressWarnings(pairwise.wilcox.test(x[, m], group, 
                               p.adjust.method = method, ...)$p.value)))
    for(i in 1:length(plist)){
      if(class(plist[[i]])=="matrix"){
        pvalue[i,] <- unlist(sapply(1:(nG-1), function(m){plist[[i]][m:(nG-1),m]}))
      } else {
        pvalue[i,] <- NA
      }
    }
    result[["Wilcox"]][[method]] <- pvalue
    col_name <- c(colnames(comb_info), paste("Wilcox", method, sep="."), 
                  colnames(pvalue))
    comb_info <- data.frame(comb_info, insert = rownames(pvalue), pvalue)
    colnames(comb_info) <- col_name
  }
  result$comb_info <- comb_info
  return(result)
}
################################################################################
#' @rdname statistics-method
#' @export
phyCtest <- function(x, formula, data = sample_data(x), ...) {
  arglist <- c(list(x = quote(x), formula = formula, data = data), list(...))
  if (!"paired" %in% names(arglist)) arglist$paired <- FALSE
  if (inherits(x, "physet"))
    x <- as.matrix(otu_table(x))
  rhs <- all.vars(formula[-2])
  lhs <- all.vars(formula[-3])
  if (all(lhs != "."))
    x <- x[, lhs]
  group <- factor(data[rownames(x), rhs[1]])
  
  col <- combn(unique(group), 2)
  nG <- length(unique(group))
  paired <- arglist$paired
  if (length(paired) == 1)
    paired <- rep(paired, ncol(col))
  summary <- matrix(0, nrow = ncol(x), ncol = nG, dimnames = 
                      list(colnames(x), paste(unique(group), "count", sep = ".")))
  chitest <- matrix(0, nrow = ncol(x), ncol = ncol(col), dimnames = 
                      list(colnames(x), paste(col[1,], col[2,], sep="vs")))
  for(i in 1:ncol(col)){
    Cdata1 <- x[group == col[1, i], ]
    Cdata2 <- x[group == col[2, i], ]
    for(j in 1:ncol(x)){
      r1 <- Cdata1[, j]
      r2 <- Cdata2[, j]
      a <- length(r1[r1 != 0])
      b <- length(r2[r2 != 0])
      c <- length(r1) - a
      d <- length(r2) - b
      if(i < nG){
        summary[j, 1] <- paste(a, length(r1), sep="|")
        summary[j, i+1] <- paste(b, length(r2), sep="|")
      }
      if (paired[i] == TRUE) {
        chitest[j, i] <- mcnemar.test(matrix(c(a,b,c,d), 2, 2))$p.value
      } else {
        chitest[j, i] <- fisher.test(matrix(c(a,b,c,d), 2, 2), ...)$p.value
      }
    }
  }
  
  result <- list()
  result$summary <- data.frame(summary)
  result$chitest <- chitest
  result$comb_info <- data.frame(summary, chitest)
  return(result)
}
################################################################################
#' @rdname statistics-method
#' @export
phyMtest <- function(x, formula, data = sample_data(x), ...) {
  arglist <- c(list(x = quote(x), formula = formula, data = data), list(...))
  if (inherits(x, "physet"))
    x <- as.matrix(otu_table(x))
  rhs <- all.vars(formula[-2])
  for (fm in rhs)
    data[, fm] <- factor(data[, fm])
  lhs <- all.vars(formula[-3])
  if (all(lhs != "."))
    x <- x[, lhs]
  data <- data[rownames(x), rhs, drop = FALSE]
  
  result <- sapply(1:ncol(x), function(i) {
    subdata <- data.frame(feat = x[, i], data)
    ft <- formula
    ft[2] <- expression(feat)
    subaov <- aov(ft, data = subdata)
#     fisher <- fisher.test(factor(x[, i]==0, levels=c(TRUE, FALSE)), data[,rhs[1]], alternative="two.sided")$p.value
    ft <- formula(paste0("feat ~ ",paste(all.vars(ft[-2]),collapse="|")))
    subkw <- coin::kruskal_test(ft, data = subdata, 
                                distribution = coin::approximate(B = 9999))
    c(anova = summary(subaov)[[1]][["Pr(>F)"]][1], 
      kruskal = coin::pvalue(subkw)[1], 
      TukeyHSD(subaov)$Group[,4])
#     if(require("multcomp")){
#       NDWD <- oneway_test(ft, data=subdata, 
#                           ytrafo=function(data) trafo(data, numeric_trafo=rank), 
#                           xtrafo=function(data) trafo(data, factor_trafo=function(x) model.matrix(~x-1) %*% t(contrMat(table(x),"Tukey"))), 
#                           teststat="max", distribution=approximate(B=90000))
#       ### global p-value
#       coin::pvalue(NDWD)[1]
#       as.vector(coin::pvalue(NDWD, method = "single-step"))
#     }
    })
  colnames(result) <- colnames(x)
  return(t(result))
}
################################################################################
#' @title Rarefaction Analysis
#' 
#' @param x A \code{\link{physet-class}} objects.
#' @param formula A formula used for constrain analysis. 
#' @param data The meta information contains RHS of formula. 
#' @param rarefy Whether using bootstrap random sampling to generate plotdata 
#' for rarefaction curves. 
#' @param permu Times for bootstrap, only used if \code{rarefy = TRUE}. 
#' @param ... other parameters. not used.
#' @return A large list with the following structure: 
#' 
#' @note
#' support shannon, simpson, pielou, chao1, ace \cr
#' Diversity statistical test \cr
#' 
#' @examples
#' data(oral)
#' rarefaction(oral, . ~ Group, rarefy = TRUE, permutation = 100)
#' @rdname rarefaction-methods
#' @export
#' 
rarefaction <- function(x, formula, data = sample_data(x), rarefy = FALSE, permu = 100, ...) {
  arglist <- c(list(x = quote(x), formula = . ~ 1, data = data, 
                    rarefy = rarefy, permu = permu), list(...))
  if (!missing(formula)) arglist$formula <- formula
  data <- eval(arglist$data)
  if (inherits(x, "physet"))
    x <- as.matrix(otu_table(x))
  
  #####   Nested bootstrap Rarefaction Curve function  #####
  sample_permute <- function(seq_list, dep_list, sample = NULL, permu = 100) {
    dep_eco <- lapply(dep_list, function(x) {
      sub_comm <- suppressWarnings(.alphaEst(t(as.matrix(table(sample(seq_list, x))))))
      replicate(permu, unlist(sub_comm), simplify = TRUE)})
    dep_mean <- sapply(dep_eco, function(x) {
      apply(x[c(1,2,4,6,7,8),], 1, mean)})
    rownames(dep_mean) <- paste(rownames(dep_mean), "mean", sep=".")
    dep_sd <- sapply(dep_eco, function(x) {
      apply(x[c(1,2,4,6,7,8),], 1, sd)/sqrt(permu)})
    rownames(dep_sd) <- paste(rownames(dep_sd), "se", sep=".")
    res <- data.frame(sample=sample, seq_dep=dep_list, 
                      t(rbind(dep_mean, dep_sd)))
    return(res)
  }
  
  res <- list()
  #####   Estimate alpha Diversity, Richness and Eveness   #####
  eco <- .alphaEst(x, split = TRUE)
  # eco_test <- para_anova_test(eco[, c(1,2,4,6,7,8)], formula, ...)
  res$alpha.diver <- eco
  # res$diver.test$grouping <- eco_test
  if (rarefy) {   ## Plot rarefaction curve for each sample
    rarefy_sample <- lapply(rownames(x), function(k) {
      seqlist <- rep(colnames(x), x[k,])
      deplist <- c(2 ^ unique(floor(log(c(1 : length(seqlist)), 2))), 
                   length(seqlist))
      sample_permute(seqlist, deplist, sample = k, permu = permu)
    })
    res$rarefy$sample <- do.call(rbind, rarefy_sample)
  }
  
  #####    Grouping analysis    #####
  rhs <- all.vars(arglist$formula[-2])
  if (length(rhs) != 0) {
    group <- data[[rhs[1]]]
    res$diver.test$pairwise <- phyTtest(as.matrix(eco[, c(1,2,4,6,7,8)]), 
                                        formula, data, p.adj = "none", ...)
    
    if (rarefy) {   ## Plot rarefaction curve for each group
      rarefy_group <- lapply(unique(group), function(k){
        seqlist <- rep(colnames(x), colSums(x[group==k, ]))
        deplist <- c(2 ^ unique(floor(log(c(1 : length(seqlist)), 2))), 
                     length(seqlist))
        sample_permute(seqlist, deplist, sample = k, permu = permu)
      })
      res$rarefy$grouping <- do.call(rbind, rarefy_group)
   }
  }
  class(res) <- "rarefy"
  return(res)
}
