ordinate = function(physeq, method="DCA", distance="bray", formula=NULL, ...){
  # If `physeq` is a formula, post deprecated notice, attempt to convert and dispatch
  if( inherits(physeq, "formula") ){
    .Deprecated(msg=paste0("First argument, `physeq`, as formula is deprecated.\n",
                           "There is now an explicit `formula` argument.\n",
                           "Please revise method call accordingly."))
    # Create the new formula, RHS-only
    formchar = as.character(physeq)
    # Error if only RHS. Formula-first syntax required both sides.
    if(length(formchar) < 3){
      stop("Need both sides of formula in this deprecated syntax... Revisit ordinate() documentation / examples.")
    }
    # Replace with (presumed) phyloseq object.
    physeq <- get(as.character(physeq)[2])
    # Create the new formula, RHS-only. 
    newFormula = as.formula(paste0("~", formchar[length(formchar)]))  
    # Dispatch to (hopefully) ordinate,phyloseq
    return(ordinate(physeq, method=method, distance=distance, formula=newFormula, ...))
  }
  # Define table of currently-supported methods
  method_table <- c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
  # List supported method names to user, if requested.
  if( inherits(physeq, "character") ){
    if( physeq=="help" ){
      cat("Available arguments to methods:\n")
      print(c(method_table))
      cat("Please be exact, partial-matching not supported.\n")
      cat("Can alternatively provide a custom distance.\n")
      cat("See:\n help(\"distance\") \n")
      return()
    } else if( physeq=="list" ){
      return(c(method_table))
    } else {
      cat("physeq needs to be a phyloseq-class object, \n")
      cat("or a character string matching \"help\" or \"list\". \n")			
    }	
  }
  # Final check that `physeq` is a phyloseq or otu_table class
  if( !inherits(physeq, "phyloseq") & !inherits(physeq, "otu_table") ){
    stop("Expected a phyloseq object or otu_table object.")
  }
  # # Start with methods that don't require 
  # #  additional distance calculation. (distance argument ignored)
  # DCA
  if( method == "DCA" ){
    return( decorana(veganifyOTU(physeq), ...) )
  }
  # CCA / RDA
  if( method %in% c("CCA", "RDA") ){
    return(cca.phyloseq(physeq, formula, method, ...))
  }
  # CAP
  if( method == "CAP" ){
    # Call/return with do.call
    return(capscale.phyloseq(physeq, formula, distance, ...))
  }
  # DPCoA
  if( method == "DPCoA" ){
    return( DPCoA(physeq, ...) )
  }  
  # # Now resort to methods that do require a separate distance/dist-calc
  # Define ps.dist. Check the class of distance argument is character or dist
  if( inherits(distance, "dist") ){
    ps.dist <- distance
  } else if( class(distance) == "character" ){
    # There are some special options for NMDS/metaMDS if distance-method
    # is supported by vegdist, so check first. If not, just calculate distance	
    vegdist_methods <- c("manhattan", "euclidean", "canberra", "bray", 
                         "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", 
                         "mountford", "raup" , "binomial", "chao")			
    # NMDS with vegdist-method to include species		
    if(method == "NMDS" & distance %in% vegdist_methods){
      return(metaMDS(veganifyOTU(physeq), distance, ...))
    }
    # Calculate distance with handoff to distance()
    ps.dist <- distance(physeq, distance, ...)
  }
  # Vanilla MDS/PCoA
  if( method %in% c("PCoA", "MDS")){
    return(pcoa(ps.dist))
  }
  # NMDS with non-vegdist-method
  if(method == "NMDS"){
    return(metaMDS(ps.dist))
  }	
}

metaMDS <- function (comm, distance = "bray", k = 2, trymax = 20, 
                     engine = c("monoMDS", "isoMDS"), autotransform = TRUE, noshare = (engine == "isoMDS"), 
          wascores = TRUE, expand = TRUE, trace = 1, plot = FALSE, 
          previous.best, ...) 
{
  engine <- match.arg(engine)
  commname <- deparse(substitute(comm))
  if (any(autotransform, noshare > 0, wascores) && any(comm < 
                                                         0)) {
    warning("'comm' has negative data: 'autotransform', 'noshare' and 'wascores' set to FALSE")
    wascores <- FALSE
    autotransform <- FALSE
    noshare <- FALSE
  }
  if (inherits(comm, "dist")) {
    dis <- comm
    if (is.null(attr(dis, "method"))) 
      attr(dis, "method") <- "user supplied"
    wascores <- FALSE
  }
  else if (length(dim(comm) == 2) && ncol(comm) == nrow(comm) && 
             all(comm == t(comm))) {
    dis <- as.dist(comm)
    attr(dis, "method") <- "user supplied"
    wascores <- FALSE
  }
  else {
    dis <- metaMDSdist(comm, distance = distance, autotransform = autotransform, 
                       noshare = noshare, trace = trace, commname = commname, 
                       ...)
  }
  if (missing(previous.best)) 
    previous.best <- NULL
  out <- metaMDSiter(dis, k = k, trymax = trymax, trace = trace, 
                     plot = plot, previous.best = previous.best, engine = engine, 
                     ...)
  if (out$stress < 0.001) {
    warning("Stress is (nearly) zero - you may have insufficient data")
  }
  points <- postMDS(out$points, dis, plot = max(0, plot - 1), 
                    ...)
  if (is.null(rownames(points))) 
    rownames(points) <- rownames(comm)
  if (wascores) {
    comm <- eval.parent(parse(text = attr(dis, "commname")))
    wa <- wascores(points, comm, expand = expand)
  }
  else wa <- NA
  out$points <- points
  out$species <- wa
  out$call <- match.call()
  if (is.null(out$data)) 
    out$data <- commname
  class(out) <- c("metaMDS", class(out))
  out
}

metaMDSdist <- function (comm, distance = "bray", autotransform = TRUE, noshare = TRUE, 
          trace = 1, commname, zerodist = "ignore", distfun = phydist, 
          ...) 
{
  if (inherits(comm, "dist") || ncol(comm) == nrow(comm) && 
        all(comm == t(comm))) 
    return(comm)
  distname <- deparse(substitute(distfun))
  print(distname)
  distfun <- match.fun(distfun)
  zerodist <- match.arg(zerodist, c("fail", "add", "ignore"))
  if (missing(commname)) 
    commname <- deparse(substitute(comm))
  xam <- max(comm)
  if (autotransform && xam > 50) {
    comm <- sqrt(comm)
    commname <- paste("sqrt(", commname, ")", sep = "")
    if (trace) 
      cat("Square root transformation\n")
  }
  if (autotransform && xam > 9) {
    comm <- vegan::wisconsin(comm)
    commname <- paste("wisconsin(", commname, ")", sep = "")
    if (trace) 
      cat("Wisconsin double standardization\n")
  }
  print(distfun)
  dis <- distfun(comm, distance, ...)
  call <- attr(dis, "call")
  call[[1]] <- as.name(distname)
  attr(dis, "call") <- call
  if (zerodist != "ignore" && any(dis <= 0, na.rm = TRUE)) {
    if (zerodist == "fail") 
      stop("Zero dissimilarities are not allowed")
    else if (zerodist == "add") {
      zero <- min(dis[dis > 0], na.rm = TRUE)/2
      dis[dis <= 0] <- zero
      if (trace) 
        cat("Zero dissimilarities changed into ", zero, 
            "\n")
    }
  }
  maxdis <- abs(distfun(matrix(c(7, 0, 0, 3), 2, 2), method = distance, 
                        ...) - 1) < 1e-04
  if ((isTRUE(noshare) && any(tmp <- no.shared(comm))) || (!is.logical(noshare) && 
                                                             noshare >= 0 && sum(tmp <- vegan::no.shared(comm))/length(dis) > 
                                                             noshare)) {
    if (trace) 
      cat("Using step-across dissimilarities:\n")
    rn <- range(dis[tmp], na.rm = TRUE)
    if (rn[2]/rn[1] > 1.01) 
      warning("non-constant distances between points with nothing shared\n", 
              "  stepacross may be meaningless: consider argument 'noshare=0'")
    is.na(dis) <- tmp
    dis <- vegan::stepacross(dis, trace = trace, toolong = 0, ...)
    if (length(unique(vegan::distconnected(tmp, trace = trace))) > 
          1) 
      warning("Data are disconnected, results may be meaningless")
  }
  attr(dis, "maxdis") <- maxdis
  attr(dis, "commname") <- commname
  attr(dis, "function") <- distname
  dis
}
  
capscale <- function (formula, data, distance = "euclidean", sqrt.dist = FALSE, 
                      comm = NULL, add = FALSE, dfun = vegdist, metaMDSdist = FALSE, 
                      na.action = na.fail, subset = NULL, ...) 
{
  EPS <- sqrt(.Machine$double.eps)
  if (!inherits(formula, "formula")) 
    stop("Needs a model formula")
  if (missing(data)) {
    data <- parent.frame()
  }
  else {
    data <- ordiGetData(match.call(), environment(formula))
  }
  print("done!")
}
  
  formula <- formula(terms(formula, data = data))
  X <- eval(formula[[2]], envir = parent.frame(), enclos = environment(formula))
  if (!inherits(X, "dist")) {
    comm <- X
    dfun <- match.fun(dfun)
    if (metaMDSdist) {
      commname <- as.character(formula[[2]])
      X <- metaMDSdist(comm, distance = distance, zerodist = "ignore", 
                       commname = commname, distfun = dfun, ...)
      commname <- attr(X, "commname")
      comm <- eval.parent(parse(text = commname))
    }
    else {
      X <- dfun(X, distance)
    }
  }
  inertia <- attr(X, "method")
  if (is.null(inertia)) 
    inertia <- "unknown"
  inertia <- paste(toupper(substr(inertia, 1, 1)), substr(inertia, 
                                                          2, 256), sep = "")
  inertia <- paste(inertia, "distance")
  if (!sqrt.dist) 
    inertia <- paste("squared", inertia)
  if (add) 
    inertia <- paste(inertia, "(euclidified)")
  fla <- update(formula, X ~ .)
  environment(fla) <- environment()
  d <- ordiParseFormula(fla, if (is.data.frame(data) && !is.null(comm)) 
    cbind(data, comm)
    else data, envdepth = 1, na.action = na.action, subset = substitute(subset))
  if (!is.null(d$subset)) 
    d$X <- d$X[, d$subset, drop = FALSE]
  if (!is.null(d$na.action)) {
    d$X <- d$X[, -d$na.action, drop = FALSE]
  }
  X <- as.dist(d$X)
  k <- attr(X, "Size") - 1
  if (sqrt.dist) 
    X <- sqrt(X)
  if (max(X) >= 4 + .Machine$double.eps) {
    inertia <- paste("mean", inertia)
    adjust <- 1
  }
  else {
    adjust <- sqrt(k)
  }
  nm <- attr(X, "Labels")
  if (add) {
    X <- cmdscale(X, k = k, eig = TRUE, add = add)
    X$eig <- X$eig[X$eig > 0]
  }
  else X <- wcmdscale(X, eig = TRUE)
  if (is.null(rownames(X$points))) 
    rownames(X$points) <- nm
  X$points <- adjust * X$points
  if (adjust == 1) {
    X$eig <- X$eig/k
    if (!is.null(X$negaxes)) 
      X$negaxes <- X$negaxes/sqrt(k)
  }
  sol <- rda.default(X$points, d$Y, d$Z, ...)
  if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
    colnames(sol$CCA$u) <- colnames(sol$CCA$biplot) <- names(sol$CCA$eig) <- colnames(sol$CCA$wa) <- colnames(sol$CCA$v) <- paste("CAP", 
                                                                                                                                  1:ncol(sol$CCA$u), sep = "")
  }
  if (!is.null(sol$CA) && sol$CA$rank > 0) {
    colnames(sol$CA$u) <- names(sol$CA$eig) <- colnames(sol$CA$v) <- paste("MDS", 
                                                                           1:ncol(sol$CA$u), sep = "")
  }
  poseig <- length(sol$CA$eig)
  if (any(X$eig < 0)) {
    negax <- X$eig[X$eig < 0]
    sol$CA$imaginary.chi <- sum(negax)
    sol$tot.chi <- sol$tot.chi + sol$CA$imaginary.chi
    sol$CA$imaginary.rank <- length(negax)
    sol$CA$imaginary.u.eig <- X$negaxes
  }
  if (!is.null(comm)) {
    comm <- scale(comm, center = TRUE, scale = FALSE)
    sol$colsum <- apply(comm, 2, sd)
    if (!is.null(d$subset)) 
      comm <- comm[d$subset, , drop = FALSE]
    if (!is.null(d$na.action)) 
      comm <- comm[-d$na.action, , drop = FALSE]
    if (!is.null(sol$pCCA) && sol$pCCA$rank > 0) 
      comm <- qr.resid(sol$pCCA$QR, comm)
    if (!is.null(sol$CCA) && sol$CCA$rank > 0) {
      v.eig <- t(comm) %*% sol$CCA$u/sqrt(k)
      sol$CCA$v <- decostand(v.eig, "normalize", MARGIN = 2)
      comm <- qr.resid(sol$CCA$QR, comm)
    }
    if (!is.null(sol$CA) && sol$CA$rank > 0) {
      v.eig <- t(comm) %*% sol$CA$u/sqrt(k)
      sol$CA$v <- decostand(v.eig, "normalize", MARGIN = 2)
    }
  }
  else {
    sol$CA$v[] <- NA
    if (!is.null(sol$CCA)) 
      sol$CCA$v[] <- NA
    sol$colsum <- NA
  }
  if (!is.null(sol$CCA) && sol$CCA$rank > 0) 
    sol$CCA$centroids <- centroids.cca(sol$CCA$wa, d$modelframe)
  if (!is.null(sol$CCA$alias)) 
    sol$CCA$centroids <- unique(sol$CCA$centroids)
  if (!is.null(sol$CCA$centroids)) {
    rs <- rowSums(sol$CCA$centroids^2)
    sol$CCA$centroids <- sol$CCA$centroids[rs > 1e-04, , 
                                           drop = FALSE]
    if (nrow(sol$CCA$centroids) == 0) 
      sol$CCA$centroids <- NULL
  }
  sol$call <- match.call()
  sol$terms <- terms(formula, "Condition", data = data)
  sol$terminfo <- ordiTerminfo(d, data)
  sol$call$formula <- formula(d$terms, width.cutoff = 500)
  sol$call$formula[[2]] <- formula[[2]]
  sol$method <- "capscale"
  if (add) 
    sol$ac <- X$ac
  sol$adjust <- adjust
  sol$inertia <- inertia
  if (metaMDSdist) 
    sol$metaMDSdist <- commname
  sol$subset <- d$subset
  sol$na.action <- d$na.action
  class(sol) <- c("capscale", class(sol))
  if (!is.null(sol$na.action)) 
    sol <- ordiNAexclude(sol, d$excluded)
  sol
}

ordiGetData <- function (call, env) {
  call$scale <- call$distance <- call$comm <- call$add <- call$dfun <- call$sqrt.dist <- call$metaMDSdist <- call$subset <- NULL
  call$na.action <- na.pass
  call[[2]] <- NULL
  call[[1]] <- as.name("model.frame")
  eval(call, env)
}
  
metaMDSdist <- function (comm, distance = "bray", autotransform = TRUE, noshare = TRUE, 
          trace = 1, commname, zerodist = "ignore", distfun = vegdist, ...) 
  {
    if (inherits(comm, "dist") || ncol(comm) == nrow(comm) && 
          all(comm == t(comm))) 
      return(comm)
    distname <- deparse(substitute(distfun))
    distfun <- match.fun(distfun)
    zerodist <- match.arg(zerodist, c("fail", "add", "ignore"))
    formals(distfun) <- c(formals(distfun), alist(... = ))
    formals(stepacross) <- c(formals(stepacross), alist(... = ))
    if (missing(commname)) 
      commname <- deparse(substitute(comm))
    xam <- max(comm)
    if (autotransform && xam > 50) {
      comm <- sqrt(comm)
      commname <- paste("sqrt(", commname, ")", sep = "")
      if (trace) 
        cat("Square root transformation\n")
    }
    if (autotransform && xam > 9) {
      comm <- wisconsin(comm)
      commname <- paste("wisconsin(", commname, ")", sep = "")
      if (trace) 
        cat("Wisconsin double standardization\n")
    }
    dis <- distfun(comm, method = distance, ...)
    call <- attr(dis, "call")
    call[[1]] <- as.name(distname)
    attr(dis, "call") <- call
    if (zerodist != "ignore" && any(dis <= 0, na.rm = TRUE)) {
      if (zerodist == "fail") 
        stop("Zero dissimilarities are not allowed")
      else if (zerodist == "add") {
        zero <- min(dis[dis > 0], na.rm = TRUE)/2
        dis[dis <= 0] <- zero
        if (trace) 
          cat("Zero dissimilarities changed into ", zero, 
              "\n")
      }
    }
    maxdis <- abs(distfun(matrix(c(7, 0, 0, 3), 2, 2), method = distance, 
                          ...) - 1) < 1e-04
    if ((isTRUE(noshare) && any(tmp <- no.shared(comm))) || (!is.logical(noshare) && 
                                                               noshare >= 0 && sum(tmp <- no.shared(comm))/length(dis) > 
                                                               noshare)) {
      if (trace) 
        cat("Using step-across dissimilarities:\n")
      rn <- range(dis[tmp], na.rm = TRUE)
      if (rn[2]/rn[1] > 1.01) 
        warning("non-constant distances between points with nothing shared\n", 
                "  stepacross may be meaningless: consider argument 'noshare=0'")
      is.na(dis) <- tmp
      dis <- stepacross(dis, trace = trace, toolong = 0, ...)
      if (length(unique(distconnected(tmp, trace = trace))) > 
            1) 
        warning("Data are disconnected, results may be meaningless")
    }
    attr(dis, "maxdis") <- maxdis
    attr(dis, "commname") <- commname
    attr(dis, "function") <- distname
    dis
  }
<environment: namespace:vegan>


cap <- function(x, formula, ...) UseMethod("cap")
setGeneric("cap")
#' @export
cap.default <- function(x, formula, ...) {
  lhs <- all.vars(formula[-3])
  rhs <- all.vars(formula[-2])
  if ("data" %in% names(list(...)) || !all(rhs %in% colnames(x))) {
    comm <- as.matrix(x)
    if (all(lhs != "."))
      comm <- comm[, lhs, drop = FALSE]
    formula <- as.formula(paste(c("comm", as.character(formula[-2])), collapse=" "))
    return(vegan::capscale(formula, ...))
  } else {
    data <- x[, rhs, drop = FALSE]
    comm <- as.matrix(x[, setdiff(colnames(x), rhs)])
    if (all(lhs != "."))
      comm <- comm[, lhs, drop = FALSE]
    formula <- as.formula(paste(c("comm", as.character(formula[-2])), collapse=" "))
    return(vegan::capscale(formula, data = data, ...))
  }
}

cca.physet <- function(x, formula, ...) {
  lhs <- all.vars(formula[-3])
  if (lhs == "." || lhs == "otu_table") {
    comm <- otu_table(x)/seqdep(x)
  } else if (lhs == "edge_matrix") {
    comm <- edge_matrix(x)
  } else if (all(lhs %in% Tnames(x))) {
    comm <- otu_table(x)[, lhs, drop=FALSE]/seqdep(x)
  } else if (all(feature %in% Enames(x))) {
    comm <- edge_matrix(x)[, lhs, drop=FALSE]
  } else {
    stop("Some taxa are not in either otu_table or edge_matrix.")
  }
  # comm <- eval(parse(text=paste(parlist[2], "(x)", sep="")))
  formula <- as.formula(paste(c("comm", as.character(formula[-2])), collapse=" "))
  if ("data" %in% names(list(...)))
    return(vegan::cca(formula, ...))
  return(vegan::cca(formula, data = sample_data(x), ...))
}
setMethod("cca", signature(x = "physet", formula = "formula"), cca.physet)

f1 <- function(x, ...) {
  arglist <- list(...)
  dfun <- match.fun(phydist)
  if ("dfun" %in% names(arglist))
    dfun <- match.fun(arglist$dfun)
  dis <- dfun(x, ...)
  return(dis)
}

f2 <- function(x, ...) {
  mdscall <- match.call()
  mdscall[[1]] <- quote(phydist)
  if ("dfun" %in% names(mdscall))
    mdscall[[1]] <- mdscall$dfun
  dis <- eval(mdscall)
  return(dis)
}

f3 <- function(x, formula, ...) {
  arglist <- list(...)
  lhs <- all.vars(formula[-3])
  rhs <- all.vars(formula[-2])
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
  arglist$formula <- formula(paste(c("dis", as.character(formula[-2])), 
                                   collapse=" "))
  arglist$method <- arglist$dfun <- NULL
  if ("data" %in% names(arglist) && !is.null(arglist$data))
    arglist$data <- quote(data)
  if ("comm" %in% names(arglist) && !is.null(arglist$comm))
    arglist$comm <- quote(comm)
  del_names <- setdiff(names(arglist),names(formals(vegan::capscale)))
  for (i in del_names) arglist[[i]] <- NULL
  ord <- eval(as.call(c(quote(vegan::capscale), arglist)))
  #  ord <- do.call(vegan::capscale, arglist)
  #  ord$call <- as.call(c(quote(vegan::capscale), arglist))
  return(ord)
}


plot_ordination = function(physeq, ordination, type="samples", axes=1:2,
                           color=NULL, shape=NULL, label=NULL, title=NULL, justDF=FALSE){
  if(length(type) > 1){
    warning("`type` can only be a single option,
            but more than one provided. Using only the first.")
    type <- type[[1]]
  }
  if(length(color) > 1){
    warning("The `color` variable argument should have length equal to 1.",
            "Taking first value.")
    color = color[[1]][1]
  }
  if(length(shape) > 1){
    warning("The `shape` variable argument should have length equal to 1.",
            "Taking first value.")
    shape = shape[[1]][1]
  }
  if(length(label) > 1){
    warning("The `label` variable argument should have length equal to 1.",
            "Taking first value.")
    label = label[[1]][1]
  }
  official_types = c("sites", "species", "biplot", "split", "scree")
  if(!inherits(physeq, "phyloseq")){
    if(inherits(physeq, "character")){
      if(physeq=="list"){
        return(official_types)
      }
    } 
    warning("Full functionality requires `physeq` be phyloseq-class ",
            "with multiple components.")
  }
  # Catch typos and synonyms
  type = gsub("^.*site[s]*.*$", "sites", type, ignore.case=TRUE)
  type = gsub("^.*sample[s]*.*$", "sites", type, ignore.case=TRUE)
  type = gsub("^.*species.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*taxa.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*OTU[s]*.*$", "species", type, ignore.case=TRUE)
  type = gsub("^.*biplot[s]*.*$", "biplot", type, ignore.case=TRUE)
  type = gsub("^.*split[s]*.*$", "split", type, ignore.case=TRUE)
  type = gsub("^.*scree[s]*.*$", "scree", type, ignore.case=TRUE)
  # If type argument is not supported...
  if( !type %in% official_types ){
    warning("type argument not supported. `type` set to 'samples'.\n",
            "See `plot_ordination('list')`")
    type <- "sites"
  }
  if( type %in% c("scree") ){
    # Stop early by passing to plot_scree() if "scree" was chosen as a type
    return( plot_scree(ordination, title=title) )
  }
  # Define a function to check if a data.frame is empty
  is_empty = function(x){
    length(x) < 2 | suppressWarnings(all(is.na(x)))
  }
  # The plotting data frames.
  # Call scores to get coordinates.
  # Silently returns only the coordinate systems available.
  # e.g. sites-only, even if species requested.
  specDF = siteDF = NULL
  trash1 = try({siteDF <- scores(ordination, choices = axes, 
                                 display="sites", physeq=physeq)},
               silent = TRUE)
  trash2 = try({specDF <- scores(ordination, choices = axes, 
                                 display="species", physeq=physeq)},
               silent = TRUE)
  # Check that have assigned coordinates to the correct object
  siteSampIntx = length(intersect(rownames(siteDF), sample_names(physeq)))
  siteTaxaIntx = length(intersect(rownames(siteDF), taxa_names(physeq)))
  specSampIntx = length(intersect(rownames(specDF), sample_names(physeq)))
  specTaxaIntx = length(intersect(rownames(specDF), taxa_names(physeq)))
  if(siteSampIntx < specSampIntx & specTaxaIntx < siteTaxaIntx){
    # Double-swap
    co = specDF
    specDF <- siteDF
    siteDF <- co
    rm(co)
  } else {
    if(siteSampIntx < specSampIntx){
      # Single swap
      siteDF <- specDF
      specDF <- NULL
    }
    if(specTaxaIntx < siteTaxaIntx){
      # Single swap 
      specDF <- siteDF
      siteDF <- NULL
    }
  }
  # If both empty, warn and return NULL
  if(is_empty(siteDF) & is_empty(specDF)){
    warning("Could not obtain coordinates from the provided `ordination`. \n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
    return(NULL)
  }
  # If either is missing, do weighted average
  if(is_empty(specDF) & type != "sites"){
    message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
    specDF <- data.frame(wascores(siteDF, w = veganifyOTU(physeq)), stringsAsFactors=FALSE)
  }
  if(is_empty(siteDF) & type != "species"){ 
    message("Species coordinates not found directly in ordination object. Attempting weighted average (`vegan::wascores`)")
    siteDF <- data.frame(wascores(specDF, w = t(veganifyOTU(physeq))), stringsAsFactors=FALSE)
  }
  # Double-check that have assigned coordinates to the correct object
  specTaxaIntx <- siteSampIntx <- NULL
  siteSampIntx <- length(intersect(rownames(siteDF), sample_names(physeq)))
  specTaxaIntx <- length(intersect(rownames(specDF), taxa_names(physeq)))
  if(siteSampIntx < 1L & !is_empty(siteDF)){
    # If siteDF is not empty, but it doesn't intersect the sample_names in physeq, warn and set to NULL
    warning("`Ordination site/sample coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
    siteDF <- NULL
  }
  if(specTaxaIntx < 1L & !is_empty(specDF)){
    # If specDF is not empty, but it doesn't intersect the taxa_names in physeq, warn and set to NULL
    warning("`Ordination species/OTU/taxa coordinate indices did not match `physeq` index names. Setting corresponding coordinates to NULL.")
    specDF <- NULL
  }
  # If you made it this far and both NULL, return NULL and throw a warning
  if(is_empty(siteDF) & is_empty(specDF)){
    warning("Could not obtain coordinates from the provided `ordination`. \n",
            "Please check your ordination method, and whether it is supported by `scores` or listed by phyloseq-package.")
    return(NULL)
  }
  if(type %in% c("biplot", "split") & (is_empty(siteDF) | is_empty(specDF)) ){
    # biplot and split require both coordinates systems available. 
    # Both were attempted, or even evaluated by weighted average.
    # If still empty, warn and switch to relevant type.
    if(is_empty(siteDF)){
      warning("Could not access/evaluate site/sample coordinates. Switching type to 'species'")
      type <- "species"
    }
    if(is_empty(specDF)){
      warning("Could not access/evaluate species/taxa/OTU coordinates. Switching type to 'sites'")
      type <- "sites"
    }
  }
  if(type != "species"){
    # samples covariate data frame, `sdf`
    sdf = NULL
    sdf = data.frame(access(physeq, slot="sam_data"), stringsAsFactors=FALSE)
    if( !is_empty(sdf) & !is_empty(siteDF) ){
      # The first two axes should always be x and y, the ordination axes.
      siteDF <- cbind(siteDF, sdf[rownames(siteDF), ])
    }
  }
  if(type != "sites"){
    # taxonomy data frame `tdf`
    tdf = NULL
    tdf = data.frame(access(physeq, slot="tax_table"), stringsAsFactors=FALSE)
    if( !is_empty(tdf) & !is_empty(specDF) ){
      # The first two axes should always be x and y, the ordination axes.
      specDF = cbind(specDF, tdf[rownames(specDF), ])
    }
  }
  # In "naked" OTU-table cases, `siteDF` or `specDF` could be matrix.
  if(!inherits(siteDF, "data.frame")){
    #warning("Sample Co-variables apparently missing in provided `physeq` for this plot-type. Coercing coord matrix to data.frame.")
    siteDF <- as.data.frame(siteDF, stringsAsFactors = FALSE)
  }  
  if(!inherits(specDF, "data.frame")){
    #warning("Taxonomy apparently missing in provided `physeq` for this plot-type. Coercing coord matrix to data.frame.")
    specDF <- as.data.frame(specDF, stringsAsFactors = FALSE)
  }
  # Define the main plot data frame, `DF`
  DF = NULL
  DF <- switch(EXPR = type, sites = siteDF, species = specDF, {
    # Anything else. In practice, type should be "biplot" or "split" here.
    # Add id.type label
    specDF$id.type <- "Taxa"
    siteDF$id.type <- "Samples"
    # But what if the axis variables differ b/w them?
    # Coerce specDF to match samples (siteDF) axis names
    colnames(specDF)[1:2] <- colnames(siteDF)[1:2]
    # Merge the two data frames together for joint plotting.
    DF = merge(specDF, siteDF, all=TRUE)
    # Replace NA with "samples" or "taxa", where appropriate (factor/character)
    if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "Samples") }
    if(!is.null(shape)){ DF <- rp.joint.fill(DF, shape, "Taxa") }
    if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "Samples") }
    if(!is.null(color)){ DF <- rp.joint.fill(DF, color, "Taxa") }
    DF
  })
  # In case user wants the plot-DF for some other purpose, return early
  if(justDF){return(DF)}
  # Check variable availability before defining mapping.
  if(!is.null(color)){ 
    if(!color %in% names(DF)){
      warning("Color variable was not found in the available data you provided.",
              "No color mapped.")
      color <- NULL
    }
  }
  if(!is.null(shape)){ 
    if(!shape %in% names(DF)){
      warning("Shape variable was not found in the available data you provided.",
              "No shape mapped.")
      shape <- NULL
    }
  }
  if(!is.null(label)){ 
    if(!label %in% names(DF)){
      warning("Label variable was not found in the available data you provided.",
              "No label mapped.")
      label <- NULL
    }
  }
  # Grab the ordination axis names from the plot data frame (as strings)
  x = colnames(DF)[1]
  y = colnames(DF)[2]   
  # Mapping section
  if( ncol(DF) <= 2){
    # If there is nothing to map, enforce simple mapping.
    message("No available covariate data to map on the points for this plot `type`")
    ord_map = aes_string(x=x, y=y)
  } else if( type %in% c("sites", "species", "split") ){
    ord_map = aes_string(x=x, y=y, color=color, shape=shape, na.rm=TRUE)
  } else if(type=="biplot"){
    # biplot, `id.type` should try to map to color and size. Only size if color specified.
    if( is.null(color) ){
      ord_map = aes_string(x=x, y=y, size="id.type", color="id.type", shape=shape, na.rm=TRUE)
    } else {
      ord_map = aes_string(x=x, y=y, size="id.type", color=color, shape=shape, na.rm=TRUE)
    }
  }
  # Plot-building section
  p <- ggplot(DF, ord_map) + geom_point(na.rm=TRUE)
  # split/facet color and shape can be anything in one or other.
  if( type=="split" ){
    # split-option requires a facet_wrap
    p <- p + facet_wrap(~id.type, nrow=1)
  }
  # If biplot, adjust scales
  if( type=="biplot" ){  
    if( is.null(color) ){
      # Rename color title in legend.
      p <- update_labels(p, list(colour="Ordination Type")) 
    } 
    # Adjust size so that samples are bigger than taxa by default.
    p <- p + scale_size_manual("type", values=c(Samples=5, Taxa=2))
  }
  # Add text labels to points
  if( !is.null(label) ){
    label_map <- aes_string(x=x, y=y, label=label, na.rm=TRUE)
    p = p + geom_text(label_map, data=rm.na.phyloseq(DF, label),
                      size=2, vjust=1.5, na.rm=TRUE)
  }
  # Optionally add a title to the plot
  if( !is.null(title) ){
    p = p + ggtitle(title)
  }
  # Add fraction variability to axis labels, if available
  if( length(extract_eigenvalue(ordination)[axes]) > 0 ){
    # Only attempt to add fraction variability
    # if extract_eigenvalue returns something
    eigvec = extract_eigenvalue(ordination)
    # Fraction variability, fracvar
    fracvar = eigvec[axes] / sum(eigvec)
    # Percent variability, percvar
    percvar = round(100*fracvar, 1)
    # The string to add to each axis label, strivar
    # Start with the curent axis labels in the plot
    strivar = as(c(p$label$x, p$label$y), "character")
    # paste the percent variability string at the end
    strivar = paste0(strivar, "   [", percvar, "%]")
    # Update the x-label and y-label
    p = p + xlab(strivar[1]) + ylab(strivar[2])
  }
  # Return the ggplot object
  return(p)
}

plot_scree = function(ordination, title=NULL){
  # Use get_eigenvalue method dispatch. It always returns a numeric vector.
  x = extract_eigenvalue(ordination)
  # Were eigenvalues found? If not, return NULL
  if( is.null(x) ){
    cat("No eigenvalues found in ordination\n")
    return(NULL)
  } else {
    # If no names, add them arbitrarily "axis1, axis2, ..., axisN"
    if( is.null(names(x)) ) names(x) = 1:length(x)
    # For scree plot, want to show the fraction of total eigenvalues
    x = x/sum(x)
    # Set negative values to zero
    x[x <= 0.0] = 0.0		
    # Create the ggplot2 data.frame, and basic ggplot2 plot
    gdf = data.frame(axis=names(x), eigenvalue = x)
    p = ggplot(gdf, aes(x=axis, y=eigenvalue)) + geom_bar(stat="identity")
    # Force the order to be same as original in x
    p = p + scale_x_discrete(limits = names(x))
    # Orient the x-labels for space.
    p = p + theme(axis.text.x=element_text(angle=90, vjust=0.5))
    # Optionally add a title to the plot
    if( !is.null(title) ){
      p <- p + ggtitle(title)
    }		
    return(p)
  }
}

extract_eigenvalue = function(ordination) UseMethod("extract_eigenvalue", ordination)
# Default is to return NULL (e.g. for NMDS, or non-supported ordinations/classes).
extract_eigenvalue.default = function(ordination) NULL
# for pcoa objects
extract_eigenvalue.pcoa = function(ordination) ordination$values$Relative_eig
# for CCA objects
extract_eigenvalue.cca = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for RDA objects
extract_eigenvalue.rda = function(ordination) c(ordination$CCA$eig, ordination$CA$eig)
# for dpcoa objects
extract_eigenvalue.dpcoa = function(ordination) ordination$eig
# for decorana (dca) objects
extract_eigenvalue.decorana = function(ordination) ordination$evals


#'    \item{}{
#'        \itemize{ordination: CA, CCA, DCA, PCA, RDA information. } {
#'        \item{ordinate}{ The ordination object generated by function 
#'                          \code{\link{ordinate}}. }
#'        \item{figures}{ The figures for the ordination object. }
#'    }
#'    \itemize{methods}{ The distance methods applied to analysis. 
#'        \item{CAP}{}
#'        \item{PCoA}{}
#'        \item{NMDS}{}
#'    }
#' }

cap.default <- function(x, formula, ...) {
  arglist <- list(...)
  lhs <- all.vars(formula[-3])
  rhs <- all.vars(formula[-2])
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
  arglist$formula <- formula(paste(c("dis", as.character(formula[-2])), 
                                   collapse=" "))
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

