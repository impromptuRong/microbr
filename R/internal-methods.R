################################################################################
# Generate Tree Node abundance matrix
# Random shuffling otu_table (Not recommended).
# Transfer objects to data.frame for easy calculation and plot
################################################################################
#' @title Generate Tree Node abundance matrix
#' @description A internal function used to generate edge structure of tree and 
#' edge abundance matrix for further analysis like unifrac.
#' 
#' @usage .cumNode(tree, ...)
#' @param tree A \code{\link[ape]{phylo}} object, or a \code{phy_tree} slot.
#' @param ... unused parameters. 
#' @return A sparse boolean \code{Matrix} stored taxa composition for nodes 
#' and leafs. 
#' 
#' @details The function is automatically called in 
#' \code{\link{physet-constructor}} when \code{phy_tree} is provided. 
#' Fast-Unifrac algorithm is used for parallel computation. 
#' 
#' @rdname cumNode
#' @import foreach
#' @export
#' @keywords internal
.cumNode <- function(tree, ...) UseMethod(".cumNode")
#' @rdname cumNode
#' @usage NULL
#' @export
.cumNode.phylo <- function(tree, ...) {
  Ntaxa <- length(tree$tip.label)
  Nedge <- nrow(tree$edge)
  edge_len <- tree$edge.length
  
  ###   extending function   ###
  extend_par <- function(edge, tree=tree, Ntaxa=Ntaxa) {
    int.nodes <- tree$edge[edge, 2]
    if (int.nodes <= Ntaxa)
      return(int.nodes)
    sons <- c()
    repeat {
      int.nodes <- tree$edge[which(tree$edge[, 1] %in% int.nodes), 2]
      if (length(int.nodes) == 0)
        break
      sons <- c(sons, int.nodes)
    }
    sons <- sons[sons <= Ntaxa]
    return(sons)
  }
  extend <- function(edge) extend_par(edge, tree, Ntaxa)
  
  ###   Non-parallel and foreach   ###
  edge_list <- foreach(MICROBR_GLOBAL_x = 1:Nedge, 
                       .packages="microbr") %dopar% extend(MICROBR_GLOBAL_x)
  #   if (isTRUE(all.equal(class(edge_list), "try-error")))
  #     system.time(edge_list <- lapply(1:Nedge, extend))
  edge_name <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]
  #   if (is.rooted(tree)) {
  #     rootlabel <- tree$node.label[which(table(c(tree$edge))==2) - Ntaxa]
  #     edge_list[[Nedge+1]] <- c(1:Nedge)
  #     edge_name <- c(edge_name, rootlabel)
  #   }
  #   names(edge_len) <- edge_name
  
  return(Matrix(sapply(edge_list, function(x){1 : Ntaxa %in% x}), 
                dimnames = list(tree$tip.label, edge_name)))
}
#' @rdname cumNode
#' @usage NULL
#' @export
setGeneric(".cumNode")
#' @rdname cumNode
setMethod(".cumNode", signature(tree="phylo"), .cumNode.phylo)
################################################################################
#' @title Random Sampling and Shuffling
#' @description Random subsampling a matrix, or a given object with abundance or 
#' count table. This command is a generic function that treat different input 
#' differently. Read \strong{Details} for more information. 
#' 
#' @param x A \code{numeric}, \code{character} \code{vector} or \code{matrix}. 
#' It can also be any object that inherits this method, like 
#' \code{\link{physet-class}}.
#' @param size The sampling pool size, sequencing depth. If \code{size} is 
#' larger than population, the input \code{x} will be returned. 
#' @param rep Repeat sampling \code{rep} times and take the average.
#' @param ... all other parameters. Not available for now.
#' @return A \code{vector}, \code{matrix}, or other object.
#' 
#' @details This function performs differently based on different input. If a 
#' \code{character vector} is given, it will random pick \code{size} elements 
#' without replacement. If a \code{integer vector} is given, the function will 
#' first create a vector for population pool based on the \code{names} and 
#' \code{count} and then random pick \code{size} elements without replacement. 
#' If a \code{double vector} is given, the function will do sampling WITH 
#' replacement by using relative abundance as probabilities. The two different
#' implementations for \code{integer} and \code{double} are generally the same 
#' when \code{rep} is sufficiently large. But when both population size and 
#' time of \code{rep} are extremely small, \code{double} method may generate 
#' odd result. Other inputs will utilize above three methods. 
#' 
#' @note Some objects like \code{\link[Matrix]{Matrix}} will automatically 
#' transfer \code{integer} to \code{double}. So shuffle the matrix first before 
#' assign it to object who use \code{Matrix} in \code{physet} if you really 
#' feel uncomfortable with using probability for \code{double}. Note that 
#' shuffle only 1 time DO increase bias while shuffled data may cause bias for 
#' richness estimation. To solve these problem, \code{physet} will normalize 
#' the otu_table with seqdep for most of the statistical analysis and will use 
#' raw otu_table for richness rarefy analysis. Also as shuffle large \code{rep} 
#' times is extremly time consuming for high-dimensional data and the result 
#' is just as same as normalizing the data by sequencing depth. So it is not 
#' recommended to use this method unless it is really necessary. 
#' 
#' @references phyloseq-paper
#' @seealso \code{\link{sample}}, \code{\link[permute]{shuffle}}
#' @examples
#' data(oral)
#' sub <- .subsampling(oral, 100, 100)
#' sub <- .subsampling(as.matrix(otu_table(oral)), 100, 10)
#' 
#' @rdname subsampling
#' @export
#' @keywords internal
.subsampling <- function(x, ...) UseMethod(".subsampling")
#' @rdname subsampling
#' @export
.subsampling.integer <- function(x, size, rep, ...) {
  if (size > sum(x))
    return(x)
  result <- tabulate(replicate(rep, sample(rep(1:length(x), x), size)), 
                     length(x))
  names(result) <- names(x)
  return(result/rep)
}
#' @rdname subsampling
#' @export
.subsampling.double <- function(x, size, rep, ...) {
  if (size > sum(x))
    return(x)
  result <- tabulate(replicate(rep, sample(1:length(x), size, 
                                           replace=TRUE, x)), length(x))
  names(result) <- names(x)
  return(result/rep)
}
#' @rdname subsampling
#' @export
.subsampling.character <- function(x, size, rep, ...) {
  if (size > length(x))
    return(c(table(x)))
  return(c(table(replicate(rep, sample(x, size, ...))))/rep)
}
#' @rdname subsampling
#' @export
.subsampling.matrix <- function(x, size, rep, ...) {
  result <- t(sapply(1:nrow(x), function(i) {
    .subsampling(x[i,], size, rep, ...)}))
  rownames(result) <- rownames(x)
  return(result)
}
#' @rdname subsampling
#' @usage NULL
#' @export
.subsampling.physet <- function(x, size, rep, ...) {
  otu_table(x) <- .subsampling(as.matrix(otu_table(x)), size, rep, ...)
  return(x)
}
#' @rdname subsampling
#' @usage NULL
#' @export
setGeneric(".subsampling")
#' @rdname subsampling
setMethod(".subsampling", signature(x="physet"), .subsampling.physet)
################################################################################
#' @title Convert Objects to Data.frame
#' @description Internal function to convert matrix, data.frame, physet objects 
#' to long-shape data.frame (sparse-Matrix format) for ploting and testing. 
#' 
#' @param x A \code{list} or a object inherits this method: \code{physet-class}. 
#' @param i,j,... The indices specifying elements in object. Indices can be 
#' +/- \code{numeric}, \code{character} and \code{logical} \code{vectors} or 
#' empty (missing). 
#' @param threshold A \code{character} string indicates threshold like: 
#' \code{"> 0"}, \code{"== 10"}, \code{"<= 0.4"}. 
#' @return A \code{data.frame}, each line indicates one record plus its row 
#' and column information in x. 
#' 
#' @seealso \code{\link{reshape}}
#' @examples
#' data(oral)
#' .m2df(as.matrix(otu_table(oral)), 1:20, -c(1:10), "> 100")
#' .m2df(oral, 1:20, -c(1:10), "> 100")
#' 
#' @rdname m2df
#' @export
#' @keywords internal
.m2df <- function(x, ...) UseMethod(".m2df")
#' @rdname m2df
#' @export
.m2df.default <- function(x, ...) data.frame(x, stringAsFactors=FALSE, ...)
#' @rdname m2df
#' @export
.m2df.matrix <- function(x, i, j, threshold = NULL, ...) {
  if (!missing(i)) x <- x[i, , drop = FALSE]
  if (!missing(j)) x <- x[, j, drop = FALSE]
  res <- data.frame(expand.grid(Row.names=rownames(x), Col.names=colnames(x), 
                                stringsAsFactors = FALSE), x=c(x))
  if (is.character(threshold))
    res <- res[eval(parse(text=paste("res$x", threshold))), ]
  return(res)
}
#' @rdname m2df
#' @usage NULL
#' @export
.m2df.physet <- function(x, i, j, threshold = NULL, ...) {
  res <- .m2df.matrix(as.matrix(otu_table(x)), i, j, ...)
  RN <- rownames(res)
  if (!is.null(sample_data(x)))
    res <- cbind(res, sample_data(x)[res$Row.names, , drop=FALSE])
  if (!is.null(tax_table(x)))
    res <- cbind(res, tax_table(x)[res$Col.names, , drop=FALSE])
  rownames(res) <- RN
  if (is.character(threshold))
    res <- res[eval(parse(text=paste("res$x", threshold))), ]
  return(res)
}
#' @rdname m2df
#' @usage NULL
#' @export
setGeneric(".m2df")
#' @rdname m2df
setMethod(".m2df", "physet", .m2df.physet)
################################################################################
#' @title Estimate alpha diversity, richness and eveness 
#' @description A function used to estimate richness for a given object. 
#' 
#' @param x A \code{matrix} or \code{physet} object. 
#' @param split treat samples separately (\code{TRUE}) or as one community. 
#' @param ... not used
#' 
#' @seealso 
#' \code{\link[vegan]{diversity}}, \code{\link[phyloseq]{estimate_richness}}. 
#' @examples
#' data(oral)
#' .alphaEst(oral)
#' @rdname alphaEst
#' @export
#' 
.alphaEst <- function(x, split = TRUE, ...){
  arglist <- c(x = quote(x), split = split, list(...))
  if(inherits(x, "physet"))
    x <- as.matrix(otu_table(x))
  if (!split)
    otu_table <- t(as.matrix(colSums(otu_table)))
 
  # Check for Input, singletons, and fraction Data.
  wtag <- FALSE
  if (!any(x == 1)) {
    warning("The input does not have singletons.")
    wtag <- TRUE
  }
  if (sum(ceiling(x)) > sum(x)) {
    warning("Ceiling Data is used to calculate richness.")
    wtag <- TRUE
  }
  if (wtag)
    warning("Richness Estimations are unreliable.
            It is suggested to use raw count table.")
  
  # Calculate standard alpha diversity, Richness and Eveness
  shannon <- vegan::diversity(x, index="shannon", MARGIN=1, base=exp(1))
  simpson <- vegan::diversity(x, index="simpson", MARGIN=1, base=exp(1))
  richness <- vegan::estimateR(ceiling(x))
  pielou <- ifelse(vegan::specnumber(x)>1, shannon/log(vegan::specnumber(x)), 0)
  res <- data.frame(t(rbind(richness, shannon, simpson, pielou)))
  return(res)
}
################################################################################
#' @title Multiple Test Correction
#' @param p probability matrix
#' @param ... other parameters passed to \code{\link{p.adjust}}
#' @rdname pcorrect
#' @export
#' @keywords internal
.pcorrect <- function(p, ...){
  return(matrix(p.adjust(c(as.matrix(p)), ...), dim(p), 
                dimnames = list(rownames(p), colnames(p))))    
}
################################################################################
#' @title Generate colors and shapes for a given vector
#' @param x A list of 2 vectors or factors, each unique element will be assigned 
#' with a unique color/shape.
#' @param value The real colors used for plotting. The function will generate 
#' a sequence of color/shape if \code{value} is missing. 
#' @param ... other parameters, not used
#' @rdname colors-method
#' @export
#' @keywords internal
.phycolors <- function(x, value, ...){
  if (!missing(value)) {
    lv <- unique(c(x[[1]], x[[2]]))
    value <- rep_len(value, length.out = length(lv))
    names(value) <- lv
  } else {
    # samples
    xs <- x[[1]]
    lvs <- setdiff(unique(xs), c("unclassified", "Unclassified", "NA"))
    hues <- c("red", "green", "blue", "purple", "cyan")
    vals <- c(hues[1:length(lvs)])
    names(vals) <- lvs
    # taxa
    xt <- x[[2]]
    lvt <- setdiff(unique(xt), c("unclassified", "Unclassified", "NA"))
    # huet <- rainbow(n+1)
    huet <- hcl(h = seq(15, 375, length = length(lvt)), l = 65, c = 100)
    valt <- c(huet[1:length(lvt)], "black", "black", "lightgrey")
    names(valt) <- c(lvt, "unclassified", "Unclassified", "NA")
    value <- c(vals, valt)
  }
  return(value)
}
################################################################################
#' @rdname colors-method
#' @export
#' @keywords internal
.physhapes <- function(x, value, ...){
  if (!missing(value)) {
    lv <- unique(c(x[[1]], x[[2]]))
    value <- rep_len(value, length.out = length(lv))
    names(value) <- lv
  } else {
    # samples
    xs <- x[[1]]
    lvs <- setdiff(unique(xs), c("unclassified", "Unclassified", "NA"))
    hues <- c(15:18, 21:25)
    vals <- c(hues[1:length(lvs)])
    names(vals) <- lvs
    # taxa
    xt <- x[[2]]
    lvt <- setdiff(unique(xt), c("unclassified", "Unclassified", "NA"))
    huet <- c(3, 7:12, 14)
    valt <- c(huet[1:length(lvt)], 4, 4, 13)
    names(valt) <- c(lvt, "unclassified", "Unclassified", "NA")
    value <- c(vals, valt)
  }
  return(value)
}
################################################################################