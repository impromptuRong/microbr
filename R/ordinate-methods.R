################################################################################
#' @title Distance and Dissimilarity 
#' @description Calculate distance or dissimilarity for community data. 
#' 
#' @param x A \code{matrix} or any object that inherits this method: like 
#' \code{\link{physet-class}}
#' @param method the distance measure to be used. Default method is "bray". 
#' Supported distances include: unifrac.uw, unifrac.w.n, unifrac.w.un, JSD and 
#' all methods registered in \code{\link[vegan]{vegdist}}.
#' @param diag logical value indicating whether the diagonal of the distance 
#' matrix should be printed. 
#' @param upper logical value indicating whether the upper triangle of the 
#' distance matrix should be printed.
#' @param ... all other parameters passed to \code{vegdist} or \code{dist}.
#' @return A distance matrix of class \code{\link{dist}}. 
#' 
#' @seealso \code{\link[vegan]{vegdist}}, \code{\link[phyloseq]{distance}}
#' @examples
#' data(oral)
#' dist <- phydist(oral, method = "unifrac.w.un")
#' 
#' @rdname phydist
#' @import foreach
#' @export
phydist <- function(x, method = "bray", diag = FALSE, upper = FALSE, ...) 
  UseMethod("phydist")
#' @rdname phydist
#' @usage NULL
#' @export
phydist.default <- function(x, method = "bray", 
                            diag = FALSE, upper = FALSE, ...) {
  arglist <- list(...)
  if (any(grep("unifrac", method)) || method == "dpcoa") {
    if (!"bw" %in% names(arglist)) {
      warning("Branch length is missing.\nSet bw = 1.")
      arglist$bw <- rep(1, ncol(x))
    }
    if (length(arglist$bw) != ncol(x)) {
      warning("Branch length mismatch table.\nTrim bw to size: ncol(x)")
      arglist$bw <- arglist$bw[1:ncol(x)]
    }
  }
  
  if (method == "unifrac.uw") {
    idx <- combn(1:nrow(x), 2)
    x <- x[, arglist$bw>0]
    arglist$bw <- arglist$bw[arglist$bw>0]
    x <- x != 0
    D <- foreach(MICROBR_GLOBAL_i = idx[1,], MICROBR_GLOBAL_j = idx[2,], 
                 .combine = 'c', .packages = "microbr") %dopar% 
      sum(arglist$bw[xor(x[MICROBR_GLOBAL_i,], x[MICROBR_GLOBAL_j,])])/
      sum(arglist$bw[x[MICROBR_GLOBAL_i,] | x[MICROBR_GLOBAL_j,]])
  } else if (method == "unifrac.w.un") {
    x <- x[, arglist$bw>0]
    arglist$bw <- arglist$bw[arglist$bw>0]
    x <- x * rep(arglist$bw, each=nrow(x))
    D <- dist(x, method = "manhattan")
  } else if (method == "unifrac.w.n") {
    idx <- combn(1:nrow(x), 2)
    li <- ncol(x) %/% 2 + 1
    ## need to fix this part ##
    leaf_depth <- rep(1, li)
    D <- foreach(MICROBR_GLOBAL_i = idx[1,], MICROBR_GLOBAL_j = idx[2,], 
                 .combine = 'c', .packages = "microbr") %dopar% 
      (abs(x[MICROBR_GLOBAL_i,]-x[MICROBR_GLOBAL_j,]) %*% arglist$bw)/
      (x[MICROBR_GLOBAL_i,1:li]+x[MICROBR_GLOBAL_j,1:li] %*% leaf_depth)
  } else if (method == "dpcoa") {
    x <- x[, arglist$bw>0]
    arglist$bw <- arglist$bw[arglist$bw>0]
    x <- x * rep(sqrt(arglist$bw), each=nrow(x))
    D <- dist(x, method = "euclidean")
  } else if (method == "JSD") {
    idx <- combn(1:nrow(x), 2)
    D <- foreach(MICROBR_GLOBAL_i = idx[1,], MICROBR_GLOBAL_j = idx[2,], 
                 .combine = 'c', .packages = "microbr") %dopar% {
                   u <- x[MICROBR_GLOBAL_i,]
                   v <- x[MICROBR_GLOBAL_j,]
                   m <- (u + v) / 2
                   return(sum(u*log(u/m), v * log(v/m), na.rm = TRUE)/2)}
  } else {
    D <- c(vegan::vegdist(as.matrix(x), method = method, ...))
  }
  attributes(D) <- list(Size = nrow(x), Labels = rownames(x), 
                        Diag = diag, Upper = upper, method = method, 
                        call = match.call(), class = "dist")
  return(D)
}
#' @rdname phydist
#' @usage NULL
#' @export
phydist.data.frame <- function(x, method = "bray", 
                               diag = FALSE, upper = FALSE, ...) {
  return(phydist.default(data.matrix(x), method = method, 
                         diag = diag, upper = upper, ...))
}
#' @rdname phydist
#' @usage NULL
#' @export
phydist.physet <- function(x, method = "bray", 
                           diag = FALSE, upper = FALSE, ...) {
  if (any(grep("unifrac", method)) || method == "dpcoa")
    return(phydist.default(edge_mat(x)/seqdep(x), method = method, 
                           diag = diag, upper = upper, bw = edge_len(x), ...))
  else
    return(phydist.default(otu_table(x)/seqdep(x), method = method, 
                           diag = diag, upper = upper, ...))
}
#' @rdname phydist
#' @usage NULL
#' @export
setGeneric("phydist")
#' @rdname phydist
setMethod("phydist", "physet", phydist.physet)
################################################################################
#' @title Ordination Analysis for Communities
#' @description This is a wrapper for supported ordination analysis. 
#' 
#' @section Usage:
#' \code{ordinate(x, apply, formula, dfun = vegdist, method = "bray", ...)} \cr
#' \code{phycca(x, formula, data, ...)} \cr
#' \code{phyrda(x, formula, data, ...)} \cr
#' \code{phycap(x, formula, data, comm, method, ...)} \cr
#' \code{phymds(x, formula, comm, method, ...)} \cr
#' \code{phydca(x, ...)} \cr
#' \code{phymca(x, ...)} \cr
#' 
#' @param x A \code{matrix}, \code{data.frame} or any object inherits these 
#' methods, like \code{\link{physet-class}}. \code{x} can also be a \code{dist} 
#' matrix if ordination method use dissimilarity matrix. 
#' @param apply The ordination analysis method. Current supported methods: CA, 
#' CCA, DCA, MCA, RDA, PCA, CAP, (N)MDS, PCoA, DPCoA. 
#' @param formula a model \code{formula} with the left hand side gives the 
#' community variables and the right hand side gives the constrains. The LHS 
#' part is matched to \code{x}, \code{comm} or \code{data} and The RHS part is 
#' matched to \code{x} or \code{data}. \dQuote{otu_table, edge_mat} can be 
#' used for \code{physet-class}. Conditioning variables can be given within a 
#' function \code{Condition}. See \code{\link[vegan]{cca}} for details. 
#' @param data A data frame contains variables on the RHS of the formula. 
#' If \code{data} is missing, the function will seek RHS in \code{x}. 
#' @param comm A (extra) community data used for calculate taxa scores based on 
#' LHS of the formula. If \code{comm} is missing, the function will seek LHS in 
#' \code{x}, then \code{data} if \code{x} is a \code{dist}. \code{comm} is not 
#' used to calculate distance but only used for taxa scores.
#' @param dfun The function used to calculate disimilarity matrix for analysis. 
#' Default method is \code{\link{phydist}} for \code{physet-class} and 
#' \code{\link{vegdist}} for default functions. The function parses \code{x} 
#' and LHS of \code{formula} to \code{dfun} if \code{x} is not a matrix. 
#' @param method The distance/dissimilarity method pass to \code{dfun}. 
#' @param ... Additional parameters pass to analysis method or dfun.
#' 
#' @return An ordination object. The class of the returned object depends upon 
#' the ordination method, the function and package that perform it. Also the 
#' object has a \code{data.frame} under \dQuote{plotdata} slot for ploting.
#' 
#' @details
#' [Partial] [Constrained] Correspondence Analysis and Redundancy Analysis \cr
#' [Partial] Constrained Analysis of Principal Coordinates \cr
#' Nonmetric Multidimensional Scaling \cr
#' Detrended Correspondence Analysis \cr
#' Multiple Correspondence Analysis \cr
#' 
#' @seealso \code{\link[phyloseq]{ordinate}}, \code{\link[vegan]{decorana}}, 
#' \code{\link[vegan]{cca}}/\code{\link[vegan]{rda}}, \code{\link[MASS]{mca}}, 
#' \code{\link[vegan]{metaMDS}}, \code{\link[vegan]{capscale}}.
#' 
#' @examples
#' data(oral)
#' ordinate(oral, "CCA", otu_table ~ Group)
#' ordinate(oral, "PCoA", method = "unifrac.uw")
#' ordinate(oral, "CAP", edge_mat ~ Group, distance = "unifrac.w.un")
#' 
#' @name ordinate
#' @rdname ordinate-method
#' @usage NULL
#' @export
ordinate <- function(x, ...) UseMethod("ordinate")
#' @rdname ordinate-method
#' @usage NULL
#' @export
ordinate.default <- function(x, apply, formula, ...) {
  if (missing(formula))
    formula <- . ~ 1
  if (apply == "NMDS" || apply == "MDS")  {
    ord <- phymds(x, formula, ...)
  } else if (apply == "PCoA") {
    formula[3] <- 1
    ord <- phycap(x, formula, ...)
  } else if (apply == "CAP") {
    ord <- phycap(x, formula, ...)
  } else if (apply == "CCA") {
    ord <- phycca(x, formula, ...)
  } else if (apply == "CA") {
    formula[3] <- 1
    ord <- phycca(x, formula, ...)
  } else if (apply == "DCA") {
    ord <- phydca(x, ...)
  } else if (apply == "MCA") {
    ord <- phymca(x, ...)
  } else if (apply == "RDA") {
    ord <- phyrda(x, formula, ...)
  } else if (apply == "PCA") {
    formula[3] <- 1
    ord <- phyrda(x, formula, ...)
  } else {
    stop("Invalid ordination method. ")
  }
  ord$coord <- list(sample = vegan::scores(ord, display = "sites"), 
                    taxa = vegan::scores(ord, display = "species"), 
                    sgof = try(vegan::goodness(ord, "sites", summarize=TRUE), 
                               silent = TRUE), 
                    tgof = try(vegan::goodness(ord, "species", summarize=TRUE), 
                               silent = TRUE),
                    eigen = try(vegan::eigenvals(ord), silent = TRUE))
  if (all(is.na(ord$coord$taxa))) ord$coord$taxa <- ord$coord$tgof <- NULL
  if (is(ord$coord$sgof, "try-error")) ord$coord$sgof <- NULL
  if (is(ord$coord$tgof, "try-error")) ord$coord$tgof <- NULL
  if (is(ord$coord$eigen, "try-error")) ord$coord$eigen <- NULL
  class(ord) <- append("ordination", class(ord))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
ordinate.physet <- function(x, apply, formula, ...){
  ordinate.default(x, apply, formula, ...)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("ordinate")
#' @rdname ordinate-method
#' @usage NULL
setMethod("ordinate", signature(x = "physet"), ordinate.physet)
################################################################################
#' @name phycca
#' @rdname ordinate-method
#' @usage NULL
#' @export
phycca <- function(x, formula, ...) UseMethod("phycca")
#' @rdname ordinate-method
#' @usage NULL
#' @export
phycca.default <- function(x, formula, ...) { 
  lhs <- all.vars(formula[-3])
  rhs <- all.vars(formula[-2])
  arglist <- c(formula = formula(paste(c("comm", as.character(formula[-2])), 
                                       collapse=" ")), list(...))
  if ("data" %in% names(arglist) || length(rhs) == 0 || 
        !all(rhs %in% colnames(x))) {
    comm <- as.matrix(x)
    if (all(lhs != "."))
      comm <- comm[, lhs, drop = FALSE]
  } else {
    data <- x[, rhs, drop = FALSE]
    arglist$data <- quote(data)
    comm <- as.matrix(x[, setdiff(colnames(x), rhs)])
    if (all(lhs != "."))
      comm <- comm[, lhs, drop = FALSE]
  }
  del_names <- setdiff(names(arglist), 
                       c("formula", "data", "na.action", "subset"))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::cca), arglist)))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phycca.physet <- function(x, formula, ...) {
  lhs <- all.vars(formula[-3])
  arglist <- c(formula = formula(paste(c("comm", as.character(formula[-2])), 
                                       collapse=" ")), list(...))
  if (lhs == "." || lhs == "otu_table") {
    comm <- otu_table(x)/seqdep(x)
  } else if (lhs == "edge_mat") {
    comm <- edge_mat(x)/seqdep(x)
  } else if (all(lhs %in% Tnames(x))) {
    comm <- otu_table(x)[, lhs, drop=FALSE]/seqdep(x)
  } else if (all(lhs %in% Enames(x))) {
    comm <- edge_mat(x)[, lhs, drop=FALSE]/seqdep(x)
  } else {
    stop("Some taxa are not in either otu_table or edge_mat.")
  }
  
  # comm <- eval(parse(text=paste(parlist[2], "(x)", sep="")))
  if (!"data" %in% names(arglist))
    arglist$data <- quote(sample_data(x))
  del_names <- setdiff(names(arglist), 
                       c("formula", "data", "na.action", "subset"))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::cca), arglist)))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("phycca")
#' @rdname ordinate-method
#' @usage NULL
setMethod("phycca", signature(x = "physet", formula = "formula"), phycca.physet)
################################################################################
#' @name phyrda
#' @rdname ordinate-method
#' @usage NULL
#' @export
phyrda <- function(x, formula, ...) UseMethod("phyrda")
#' @rdname ordinate-method
#' @usage NULL
#' @export
phyrda.default <- function(x, formula, ...) {
  lhs <- all.vars(formula[-3])
  rhs <- all.vars(formula[-2])
  arglist <- c(formula = formula(paste(c("comm", as.character(formula[-2])), 
                                       collapse=" ")), list(...))
  if ("data" %in% names(arglist) || length(rhs) == 0 
      || !all(rhs %in% colnames(x))) {
    comm <- as.matrix(x)
    if (all(lhs != "."))
      comm <- comm[, lhs, drop = FALSE]
  } else {
    data <- x[, rhs, drop = FALSE]
    arglist$data <- quote(data)
    comm <- as.matrix(x[, setdiff(colnames(x), rhs), drop = FALSE])
    if (all(lhs != "."))
      comm <- comm[, lhs, drop = FALSE]
  }
  del_names <- setdiff(names(arglist), 
                       c("formula", "data", "scale", "na.action", "subset"))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::rda), arglist)))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phyrda.physet <- function(x, formula, ...) {
  lhs <- all.vars(formula[-3])
  arglist <- c(formula = formula(paste(c("comm", as.character(formula[-2])), 
                                       collapse=" ")), list(...))
  if (lhs == "." || lhs == "otu_table") {
    comm <- otu_table(x)/seqdep(x)
  } else if (lhs == "edge_mat") {
    comm <- edge_mat(x)/seqdep(x)
  } else if (all(lhs %in% Tnames(x))) {
    comm <- otu_table(x)[, lhs, drop=FALSE]/seqdep(x)
  } else if (all(lhs %in% Enames(x))) {
    comm <- edge_mat(x)[, lhs, drop=FALSE]/seqdep(x)
  } else {
    stop("Some taxa are not in either otu_table or edge_mat.")
  }
  # comm <- eval(parse(text=paste(parlist[2], "(x)", sep="")))
  if (!"data" %in% names(arglist))
    arglist$data <- quote(sample_data(x))
  del_names <- setdiff(names(arglist), 
                       c("formula", "data", "scale", "na.action", "subset"))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::rda), arglist)))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("phyrda")
#' @rdname ordinate-method
#' @usage NULL
setMethod("phyrda", signature(x = "physet", formula = "formula"), phyrda.physet)
################################################################################
#' @name phycap
#' @rdname ordinate-method
#' @usage NULL
#' @export
phycap <- function(x, formula, ...) UseMethod("phycap")
#' @rdname ordinate-method
#' @usage NULL
#' @export
phycap.default <- function(x, formula, ...) {
  lhs <- all.vars(formula[-3])
  rhs <- all.vars(formula[-2])
  arglist <- c(formula = formula(paste(c("dis", as.character(formula[-2])), 
                                       collapse=" ")), list(...))
  if (!"data" %in% names(arglist) && length(rhs) != 0 && all(rhs %in% colnames(x))) {
    arglist$data <- x[, rhs, drop = FALSE]
    x <- as.matrix(x[, setdiff(colnames(x), rhs)])  
  }
  if (!"comm" %in% names(arglist))
    arglist$comm <- as.matrix(x)
  if (!is.null(arglist$comm) && all(lhs != "."))
    arglist$comm <- arglist$comm[, lhs, drop = FALSE]
  data <- arglist$data
  comm <- arglist$comm
  
  dfun <- match.fun(phydist)
  if ("dfun" %in% names(arglist))
    dfun <- match.fun(arglist$dfun)
  dis <- dfun(x, ...)
  arglist$method <- arglist$dfun <- NULL
  if ("data" %in% names(arglist) && !is.null(arglist$data)){
    arglist$data <- quote(data)
    data <- data.frame(data)
  }
  if ("comm" %in% names(arglist) && !is.null(arglist$comm)){
    arglist$comm <- quote(comm)
    comm <- as.matrix(comm)
  }
  del_names <- setdiff(names(arglist),names(formals(vegan::capscale)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::capscale), arglist)))
  #  ord <- do.call(vegan::capscale, arglist)
  #  ord$call <- as.call(c(quote(vegan::capscale), arglist))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phycap.dist <- function(x, formula, ...) {
  lhs <- all.vars(formula[-3])
  rhs <- all.vars(formula[-2])
  arglist <- c(formula = formula(paste(c("dis", as.character(formula[-2])), 
                                       collapse=" ")), list(...))
  if (!"comm" %in% names(arglist) && "data" %in% names(arglist)) {
    comm <- arglist$data
    if (all(lhs != ".")) {
      comm <- comm[, lhs, drop = FALSE]
    } else {
      comm <- comm[, setdiff(colnames(x), rhs), drop = FALSE]
    }
    arglist$comm <- comm
  }
  data <- arglist$data
  comm <- arglist$comm
  
  dis <- x
  arglist$method <- arglist$dfun <- NULL
  if ("data" %in% names(arglist) && !is.null(arglist$data)){
    arglist$data <- quote(data)
    data <- data.frame(data)
  }
  if ("comm" %in% names(arglist) && !is.null(arglist$comm)){
    arglist$comm <- quote(comm)
    comm <- as.matrix(comm)
  }
  del_names <- setdiff(names(arglist),names(formals(vegan::capscale)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::capscale), arglist)))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phycap.physet <- function(x, formula, ...) {
  lhs <- all.vars(formula[-3])
  rhs <- all.vars(formula[-2])
  arglist <- c(formula = formula(paste(c("dis", as.character(formula[-2])), 
                                       collapse=" ")), list(...))
  if (!"data" %in% names(arglist))
    arglist$data <- sample_data(x)
  if (!"comm" %in% names(arglist)) {
    if (lhs == "." || lhs == "otu_table") {
      arglist$comm <- as.matrix(otu_table(x)/seqdep(x))
    } else if (lhs == "edge_mat") {
      arglist$comm <- as.matrix(edge_mat(x))/seqdep(x)
    } else if (all(lhs %in% Tnames(x))) {
      arglist$comm <- as.matrix(otu_table(x)[, lhs, drop=FALSE]/seqdep(x))
    } else if (all(lhs %in% Enames(x))) {
      arglist$comm <- as.matrix(edge_mat(x)[, lhs, drop=FALSE])/seqdep(x)
    } else {
      stop("Formula variables mismatch!")
    }
  }
  data <- arglist$data
  comm <- arglist$comm
  
  dfun <- match.fun(phydist)
  if ("dfun" %in% names(arglist))
    dfun <- match.fun(arglist$dfun)
  dis <- dfun(x, ...)
  lhs <- all.vars(formula[-3])
  arglist$method <- arglist$dfun <- NULL
  if ("data" %in% names(arglist) && !is.null(arglist$data)){
    arglist$data <- quote(data)
    data <- data.frame(data)
  }
  if ("comm" %in% names(arglist) && !is.null(arglist$comm)){
    arglist$comm <- quote(comm)
    comm <- as.matrix(comm)
  }
  del_names <- setdiff(names(arglist),names(formals(vegan::capscale)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::capscale), arglist)))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("phycap")
#' @rdname ordinate-method
#' @usage NULL
setMethod("phycap", signature(x = "physet", formula = "formula"), phycap.physet)
################################################################################
#' @name phymds
#' @rdname ordinate-method
#' @usage NULL
#' @export
phymds <- function(x, ...) UseMethod("phymds")
#' @rdname ordinate-method
#' @usage NULL
#' @export
phymds.default <- function(x, formula, ...) {
  # Retrieve community matrix
  arglist <- list(...)
  comm <- as.matrix(x)
  if ("comm" %in% names(arglist))
    comm <- arglist$comm
  if (!is.null(comm) && !missing(formula)){
    lhs <- formula
    if (inherits(formula, "formula"))
      lhs <- all.vars(formula[-3])
    if (all(lhs != "."))
      comm <- comm[, lhs, drop=FALSE]
  }
  
  # Calculate distance matrix, default is phydist method
  dfun <- match.fun(phydist)
  if ("dfun" %in% names(arglist))
    dfun <- match.fun(arglist$dfun)
  dis <- dfun(x, ...)
  
  # Calculate MDS ordination
  arglist$comm <- quote(dis)
  del_names <- setdiff(names(arglist),names(formals(vegan::metaMDS)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- do.call(vegan::metaMDS, arglist)
  ord$call <- match.call()
  # Calculate species scores
  if (!is.null(comm)) {
    expand <- TRUE
    if ("expand" %in% names(arglist))
      expand <- arglist$expand
    ord$species <- vegan::wascores(ord$points, w = comm, expand = expand)
  }
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phymds.dist <- function(x, formula, ...) {
  # Retrieve community matrix
  arglist <- list(...)
  comm <- arglist$comm
  if (!is.null(comm) && !missing(formula)) {
    lhs <- formula
    if (inherits(formula, "formula"))
      lhs <- all.vars(formula[-3])
    if (all(lhs != "."))
      comm <- comm[, lhs, drop = FALSE]
  }
  
  # Calculate MDS ordination
  arglist$comm <- quote(x)
  del_names <- setdiff(names(arglist),names(formals(vegan::metaMDS)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- do.call(vegan::metaMDS, arglist)
  ord$call <- match.call()
  # Calculate species scores
  if (!is.null(comm)) {
    expand <- TRUE
    if ("expand" %in% names(arglist))
      expand <- arglist$expand
    ord$species <- vegan::wascores(ord$points, w = comm, expand = expand)
  }
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phymds.physet <- function(x, formula, ...) {
  # Retrieve community matrix
  arglist <- list(...)
  if (!"comm" %in% names(arglist)) {
    lhs <- "."
    if (!missing(formula)){
      lhs <- formula
      if (inherits(formula, "formula"))
        lhs <- all.vars(formula[-3])
    }
    if (any(lhs == ".")) {
      comm <- as.matrix(otu_table(x)/seqdep(x))
    } else if (lhs == "edge_mat") {
      comm <- as.matrix(edge_mat(x))/seqdep(x)
    } else if (all(lhs %in% Tnames(x))) {
      comm <- comm[, lhs, drop=FALSE]
    } else if (all(lhs %in% Enames(x))) {
      comm <- as.matrix(edge_mat(x)[, lhs, drop=FALSE])/seqdep(x)
    } else {
      stop("Formula do not match input data! ")
    }
  } else {
    comm <- arglist$comm
    if (!is.null(comm) && !missing(formula)){
      lhs <- formula
      if (inherits(formula, "formula"))
        lhs <- all.vars(formula[-3])
      if (all(lhs != "."))
        comm <- comm[, lhs, drop=FALSE]
    }
  }
  
  # Calculate distance matrix, default is phydist
  dfun <- match.fun(phydist)
  if ("dfun" %in% names(arglist))
    dfun <- match.fun(arglist$dfun)
  dis <- dfun(x, ...)
  
  # Calculate MDS ordination
  arglist$comm <- quote(dis)
  del_names <- setdiff(names(arglist),names(formals(vegan::metaMDS)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- do.call(vegan::metaMDS, arglist)
  ord$call <- match.call()
  
  print(ord)
  # Calculate species scores
  if (!is.null(comm)) {
    expand <- TRUE
    if ("expand" %in% names(arglist))
      expand <- arglist$expand
    ord$species <- vegan::wascores(ord$points, w = comm, expand = expand)
  }
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("phymds")
#' @rdname ordinate-method
#' @usage NULL
setMethod("phymds", signature(x = "physet"), phymds.physet)
################################################################################
#' @name dca
#' @rdname ordinate-method
#' @usage NULL
#' @export
phydca <- function(x, ...) UseMethod("dca")
#' @rdname ordinate-method
#' @usage NULL
#' @export
phydca.default <- function(x, ...) {
  arglist <- c(veg = quote(x), list(...))
  del_names <- setdiff(names(arglist), names(formals(vegan::decorana)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::decorana), arglist)))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phydca.physet <- function(x, ...) {
  x <- as.matrix(otu_table(x)/seqdep(x))
  arglist <- c(veg = quote(x), list(...))
  del_names <- setdiff(names(arglist), names(formals(vegan::decorana)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::decorana), arglist)))
  return(ord)
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("phydca")
#' @rdname ordinate-method
#' @usage NULL
setMethod("phydca", "physet", phydca.physet)
################################################################################
#' @name phypcoa
#' @rdname ordinate-method
#' @usage NULL
#' @export
phypcoa <- function(x, ...) UseMethod("phypcoa")
#' @rdname ordinate-method
#' @usage NULL
#' @export
phypcoa.dist <- function(x, ...) {
  return(cmdscale(x, ...))
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phypcoa.default <- function(x, ...) {
  return(cmdscale(phydist(x, ...), ...))
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phypcoa.physet <- function(x, ...) {
  return(cmdscale(phydist(x, ...), ...))
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("phypcoa")
#' @rdname ordinate-method
#' @usage NULL
setMethod("phypcoa", "physet", phypcoa.physet)
################################################################################
#' @name phymca
#' @rdname ordinate-method
#' @usage NULL
#' @import MASS
#' @export
phymca <- function(x, ...) UseMethod("phymca")
#' @rdname ordinate-method
#' @usage NULL
#' @export
phymca.default <- function(x, ...) {
  return(MASS::mca(x, ...))
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
phymca.physet <- function(x, ...) {
  return(MASS::mca(as.matrix(otu_table(x)/seqdep(x)), ...))
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("phymca")
#' @rdname ordinate-method
#' @usage NULL
setMethod("phymca", "physet", phymca.physet)
################################################################################
