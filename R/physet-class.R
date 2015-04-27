################################################################################
# The microbr package is a general purpose wrapper and pipeline for some routine 
# microbiome and metagenomics analysis. Similar package (like phyloseq) has been 
# developed and well maintained in Github and Bioconductor. 
# This package is used in Dong's Lab only or for personal learning purpose. Some
# of the codes are copied or highly similar to other contributions. Please DO
# CONTACT original authors for permission if codes are used for other purpose.
# 
# What's in this physet-class.R file: 
# Define virtual classes inherited from other packages
# Define a S4 class physet with 4 primary slots and 4 additional slots.
# Primary slots: otu_table, sam_table, tax_table, seqdep
# Additional slots: phy_tree, edge_com, edge_len, edge_mat
# Initialization method, show methods, assign methods, replace methods
################################################################################
#' @title Inherited and Virtual Classes
#' @description Define virtual classes and inherit \code{ape::phylo-class}
#' @usage NULL
#' @details S3 Class \code{\link[ape]{phylo}} is not formally defined in ape 
#' package. In order to inherit all the methods in \code{\link{ape}} package, 
#' an empty list is registered as \dQuote{phylo} and extended as a S4 virtual 
#' class allow \code{NULL}. 
#' 
#' @seealso \code{\link[ape]{phylo}}, \code{\link{setOldClass}}
#' @name phylo-class
#' @rdname virtual-class
#' @format NULL
#' @import methods
#' @import Matrix
#' @exportClass phylo
#' @keywords internal
# @importClassesFrom ape phylo
phylo <- structure(list(), class = "phylo")
setOldClass("phylo")
setClassUnion("phyloOrNULL", c("phylo", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
# setClassUnion("numericOrNULL", c("numeric", "NULL"))
# setClassUnion("MatrixOrNULL", c("Matrix", "NULL"))
################################################################################
#' @title The S4 class for microbiome analysis
#'
#' @description The \code{physet-class} is used to store most of the microbiome 
#' components in a single instance. The physet class allows easy and flexible 
#' manipulation on sample/taxa and it can facilitate downstream analyses. 
#' 
#' The \code{\link{physet-constructor}} ensures that different data components 
#' are compatible. Initially slots are set to \code{NULL}. The 4 primary slots 
#' can be updated with corresponding assignment methods while additional slots 
#' can be updated with phylogenetic tree. Also \code{\link{extract}} method can 
#' be used to exclude unwanted sample/taxa from downstream analyses.
#' 
#' The \code{physet-class} contains 4 primary slots. They are the fundemental 
#' components for most of the downstream analyses. 
#' Primary slots:
#' @slot otu_table A [\code{Nsample * Ntaxa}] \code{\link[Matrix]{Matrix}} for 
#'        count table or abundance table. The otu_table will keep raw count or 
#'        abundance and it will be normalized by \code{seqdep} automatically
#'        for most statistical analysis. The otu_table could also be randomly 
#'        subsampled with \code{\link{.subsampling}} (not recommanded, see 
#'        reference). 
#'        New otu_table can be assigned with \code{\link{otu_table<-}}. 
#' @slot sample_data A \code{\link{data.frame}} for meta data information. This 
#'        slot is matched accroding to \code{rownames(otu_table)}. So samples 
#'        only occured in \code{otu_table} will have \code{NA} attributes while 
#'        samples not occured in \code{otu_table} will be eliminated. One can 
#'        use \code{\link{extract}} method to manipulate unwanted samples. New 
#'        \code{sample_data} can be assigned by \code{\link{sample_data<-}}.
#' @slot tax_table A \code{\link{data.frame}} for taxonomy information. This 
#'        slot is matched accroding to \code{colnames(otu_table)}. So taxa only 
#'        occured in \code{otu_table} will have \code{NA} attributes while taxa 
#'        not occured in \code{otu_table} will be eliminated. One can use 
#'        \code{\link{extract}} method to eliminate unwanted taxa. New taxonomy 
#'        information can be assigned with \code{\link{tax_table<-}}.
#' @slot seqdep Sequencing depth of each sample for normalization. In order to 
#'        keep analysis under same scope, this slot will not be changed when 
#'        you eliminate unwanted taxa. So use \code{\link{seqdep<-}} to modify 
#'        it if you want to re-normalize \code{otu_table}, \code{edge_mat}. 
#' 
#' When phylogenetic tree is provided, 4 additional slots are added for 
#' phylogeny based analyses like unifrac and dpcoa. 
#' Additional slots:
#' @slot phy_tree A \code{\link[ape]{phylo}}-class from the ape package for 
#'        a rooted phylogenetic tree. Assign value to \code{phy_tree} will 
#'        not modify the 6 primary slot. This approach guarantee that tree 
#'        is kept without trimming and it will not influence non-phylogeny 
#'        analysis. \code{physet} will not modify \code{phy_tree}, please 
#'        check \pkg{\link[ape]{ape}} package for tree manipulation.
#' @slot edge_len A named vector stored the branch length for each edge in 
#'        \code{phy_tree}. 
#' @slot edge_com A boolean \code{\link[Matrix]{Matrix}} stored hierarchcal 
#'        structure of input tree. It can simplify some analysis and it is 
#'        useful to generate supporting files for iTOL plot. 
#' @slot edge_mat A [\code{Nsample * Nedge}] \code{\link{Matrix}} records the 
#'        total taxa abundance under each node for each sample. Any change in 
#'        \code{otu_table} and \code{phy_tree} will modify the value. The 
#'        matrix is used to speed up unifrac, dpcoa diversity calculation and 
#'        is used for phylogeny based feature selection. 
#' 
#' @seealso
#'  The constructor \code{\link{physet-constructor}}, 
#'  accessors \code{\link{otu_table}}, \code{\link{sample_data}}, 
#'  \code{\link{tax_table}}, \code{\link{phy_tree}}.
#' 
#' @examples
#' getSlots("physet")
#' slotNames("physet")
#' @name physet-class
#' @rdname physet-class
#' @exportClass physet
physet <- setClass(
  Class = "physet",
  slots = list(
    otu_table = "Matrix",
    sample_data = "listOrNULL",
    tax_table = "listOrNULL",
    phy_tree = "phyloOrNULL",
    edge_len = "numeric", 
    edge_com = "Matrix", 
    edge_mat = "Matrix", 
    seqdep = "numeric"
  ),
  prototype = list(
    otu_table = NULL,
    sample_data = NULL,
    tax_table = NULL,
    phy_tree = NULL, 
    edge_len = NULL, 
    edge_com = NULL, 
    edge_mat = NULL, 
    seqdep = 1
  )
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
)
################################################################################
#' Build a physet-class from its components
#'
#' The \code{physet-constructor} method is the initialization method to create 
#' a valid \code{\link{physet-class}}. This function firstly manage 4 major 
#' components defined in \code{\link{physet-class}}. Then if a phylogenetic 
#' tree is given, the internal function \code{\link{.cumNode}} is called to set 
#' up 3 additional components. Unprovided slots will be kept as \code{NULL} and 
#' can be assigned later. Also any physet instance can be re-build with the 
#' constructor with new slots.
#' 
#' @section Usage: 
#' physet(...) \cr
#' ## S4 method for signature 'physet' \cr
#' initialize(.Object, ...) \cr
#' 
#' @param .Object The physet object to create. 
#' @param ... Components to be passed to the constructor, 4 major components, 
#' phylogenetic tree and another \code{physet} class are allowed as input. 
#' The constructor has no order for inputs. Specify input parameter with its 
#' slot name in physet: "otu_table=yourdata, phy_tree=yourtree". 
#'  \itemize{
#'    \item{physet}{ a \code{\link{physet}} object. Further arguments will 
#'          be used to update this object. Use \code{physet = }. }
#'    \item{otu_table}{ a \code{matrix}, \code{Matrix} or \code{data.frame} 
#'          for taxa table. Use \code{otu_table = }. }
#'    \item{sample_data}{ a \code{\link{data.frame}} for sample information. 
#'          All \code{factor} will be flatterned as \code{atomic} and will be 
#'          transfered back when needed. This will avoid null level problem.
#'          If \code{physet} is provided, the \code{sample_data} will expand 
#'          not replace the original \code{sample_data} in \code{physet}. Only 
#'          conflict attributes will be overwritten. Use \code{sample_data = }. }
#'    \item{tax_table}{ a \code{\link{data.frame}} for taxonomy information. 
#'          All \code{factor} will be flatterned as \code{atomic} and will be 
#'          transfered back when needed. This will avoid null level problem.
#'          If \code{physet} is provided, the \code{sample_data} will expand 
#'          not replace the original \code{sample_data} in \code{physet}. Only 
#'          conflict attributes will be overwritten. Use \code{tax_table = }. }
#'    \item{seqdep}{ a \code{numeric} value or a \code{numeric vector} with 
#'          length \code{Nsample}. By default, the constructor will use 
#'          \code{rowSums(otu_table)} as \code{seqdep} if not specified. Use 
#'          \code{seqdep = }. }
#'  }
#'
#' @return A physet object is returned. All slots are \code{NULL} if nothing is 
#'  passed. Otherwise returned object will at least contains 4 primary slots. 
#'
#' @seealso \code{\link{access}}, \code{\link{assign}}, \code{\link{extract}}, 
#' @examples
#' require(foreach)
#' data(oral_raw)
#' oral <- physet(otu_table=oral_raw$rawdata, sample_data=oral_raw$metainfo, 
#'                tax_table=oral_raw$taxonomy, phy_tree=oral_raw$rawtree)
#' oral
#' slotNames(oral)
#' @aliases physet-constructor
#' @rdname physet-constructor
#' @usage NULL
#' @export physet
setMethod("initialize", "physet", function(.Object, ...) {
  arglist <- list(...)
  cat("New physet object:")
  if ("physet" %in% names(arglist)) {
    .Object@otu_table <- arglist$physet@otu_table
    .Object@sample_data <- arglist$physet@sample_data
    .Object@tax_table <- arglist$physet@tax_table
    .Object@seqdep <- arglist$physet@seqdep
  }
  data <- .Object@otu_table
  rowN <- rownames(data)
  colN <- colnames(data)
  
  if ("otu_table" %in% names(arglist)) {
    cat("\nRegister otu_table ... ")
    data <- Matrix(data.matrix(arglist$otu_table))
    rownames(data) <- gsub("[^a-zA-Z0-9]", "_", rownames(data))
    colnames(data) <- gsub("[^a-zA-Z0-9]", "_", colnames(data))
    rowN <- rownames(data)
    colN <- colnames(data)
    .Object@otu_table <- data
    .Object@seqdep <- rowSums(data) + 1e-20
    if (!is.null(.Object@sample_data) && !"sample_data" %in% names(arglist))
      arglist$sample_data <- .Object@sample_data
    .Object@sample_data <- data.frame(Sam_ID = rowN, row.names = rowN, 
                                      stringsAsFactors = FALSE)
    if (!is.null(.Object@tax_table) && !"tax_table" %in% names(arglist))
      arglist$tax_table <- .Object@tax_table
    .Object@tax_table <- data.frame(Tax_ID = colN, row.names = colN, 
                                    stringsAsFactors = FALSE)
    cat("Done!")
  }
  
  if ("seqdep" %in% names(arglist))
    .Object@seqdep <- rep(c(arglist$seqdep), length(rowN)/length(arglist$seqdep))
  
  if ("sample_data" %in% names(arglist) && !is.null(arglist$sample_data)) {
    cat("\nRegister sample_data ... ")
    arglist$sample_data <- data.frame(lapply(arglist$sample_data, as.vector), 
                                      row.names=rownames(arglist$sample_data), 
                                      stringsAsFactors = FALSE)
    rownames(arglist$sample_data) <- gsub("[^a-zA-Z0-9]", "_", 
                                          rownames(arglist$sample_data))
    arglist$sample_data <- merge(.Object@sample_data, arglist$sample_data, 
                                 by="row.names", all.x=TRUE, sort=FALSE, 
                                 suffixes = c("._rm",""))
    rownames(arglist$sample_data) <- arglist$sample_data$Row.names
    exc <- c(1, grep("\\._rm$", names(arglist$sample_data)))
    .Object@sample_data <- arglist$sample_data[rowN, -exc, drop=FALSE]
    cat("Done!")
  }
  
  if ("tax_table" %in% names(arglist) && !is.null(arglist$tax_table)) {
    cat("\nRegister tax_table ... ")
    arglist$tax_table <- data.frame(lapply(arglist$tax_table, as.vector), 
                                    row.names=rownames(arglist$tax_table), 
                                    stringsAsFactors = FALSE)
    rownames(arglist$tax_table) <- gsub("[^a-zA-Z0-9]", "_", 
                                        rownames(arglist$tax_table))
    arglist$tax_table <- merge(.Object@tax_table, arglist$tax_table, 
                               by="row.names", all.x=TRUE, sort=FALSE, 
                               suffixes = c("._rm",""))
    rownames(arglist$tax_table) <- arglist$tax_table$Row.names
    exc <- c(1, grep("\\._rm$", names(arglist$tax_table)))
    .Object@tax_table <- arglist$tax_table[colN, -exc, drop=FALSE]
    cat("Done!")
  }
  
  update_edge_mat <- FALSE
  if ("phy_tree" %in% names(arglist)) {
    cat("\nRegister phy_tree ... ")
    .Object@phy_tree <- arglist$phy_tree
    if (!is.null(.Object@phy_tree)) {
      # is.rooted(raw.tree)
      if (is.null(.Object@phy_tree$node.label)) 
        .Object@phy_tree$node.label <- paste("Branch", 1:.Object@phy_tree$Nnode, 
                                             sep="_")
      .Object@phy_tree$tip.label <- gsub("[^a-zA-Z0-9]", "_", 
                                         .Object@phy_tree$tip.label)
      .Object@edge_com <- .cumNode(.Object@phy_tree)
      update_edge_mat <- TRUE
    }
    cat("Done!")
  } else {
    if ("physet" %in% names(arglist) && !is.null(arglist$physet@phy_tree)) {
      .Object@phy_tree <- arglist$physet@phy_tree
      .Object@edge_com <- arglist$physet@edge_com
      .Object@edge_len <- arglist$physet@edge_len
      .Object@edge_mat <- arglist$physet@edge_mat
      if ("otu_table" %in% arglist)
        update_edge_mat <- TRUE
    }
  }
  
  # update edge_len and edge_mat
  if (update_edge_mat) {
    cat("\nRegister edge_mat ... ")
    # Select common taxa and common node
    tree <- .Object@phy_tree
    c_taxa <- intersect(tree$tip.label, colN)
    if (length(c_taxa) < length(tree$tip.label)) {
      warning("Not all taxa in tree are included in otu_table")
      tree <- ape::drop.tip(tree, setdiff(tree$tip.label, c_taxa))
    }
    if (length(c_taxa) < length(colN))
      warning("Not all taxa in otu_table are included in tree")
    c_taxa <- tree$tip.label
    c_node <- c(tree$tip.label, tree$node.label)[tree$edge[,2]]
    print(c_node)
    print(c_taxa)
    .Object@edge_len <- tree$edge.length
    .Object@edge_mat <- .Object@otu_table[, c_taxa, drop=FALSE] %*% 
      .Object@edge_com[c_taxa, c_node, drop=FALSE]
    names(.Object@edge_len) <- c_node
    cat("Done!")
  }
  return(.Object)
})
################################################################################
#' @title Show method for physet object
#' @description The function to show or print a physet object.
#' @param object A physet object
#' @seealso \code{\link{summary}}
#' @examples
#' data(oral)
#' oral
#' print(oral)
#' @rdname display-method
#' @export
#' @keywords internal
setMethod("show", "physet", function(object) {
  if (is.null(object@otu_table)) {
    cat("Empty physet object")
    return(NULL)
  }
  
  cat("physet-class object\nData info: ", 
      Nsample(object), " samples; ", Ntaxa(object), " taxa", sep="")
  if (!is.null(object@phy_tree))
    cat("; ", Nedge(object), " edges", sep="") 
  cat("\n")
  
  if (!is.null(object@sample_data))
    cat("Meta info: ", ncol(object@sample_data), " attributes\n", sep="")
  if (!is.null(object@tax_table))
    cat("Taxa info: ", ncol(object@tax_table), " attributes\n", sep="")
  if (!is.null(object@phy_tree))
    cat("Tree info: ", length(object@phy_tree$tip.label), " leafs; ", 
        object@phy_tree$Nnode, " internal nodes\n           ",
        ifelse(ape::is.rooted(object@phy_tree), "Rooted; ", "Unrooted; "),
        ifelse(is.null(object@phy_tree$edge.length), 
               "no branch lengths", 
               "includes branch lengths"),
        sep="")
})
################################################################################
#' @title Access slots or attributes in physet object. 
#' @description Access slots in physet object by slot name. Or get attributes in 
#' sample_data and tax_table like name and size. 
#' 
#' @section Usage:
#' access(x, slot_name) \cr
#' slot_name(x) \cr
#' x$name \cr
#' 
#' @param x A physet object
#' @param slot_name A character indicating the slot name: \code{otu_table}, 
#' \code{sample_data}, \code{tax_table}, \code{seqdep}, \code{phy_tree}, 
#' \code{edge_com}, \code{edge_len}, \code{edge_mat}. The function will return 
#' \code{NULL} if \code{slot_name} does not exist. \code{Nsample}, \code{Ntaxa}, 
#' \code{Nedge} are allowed to access the count of samples, taxa and utilized  
#' edges. \code{Snames}, \code{Tnames} and \code{Enames} can return sample 
#' names, taxa names and utilized tree edge names in physet object. 
#' @param name A character specifying variables. \code{name} can be any value 
#' in \code{names(sample_data(x))} and \code{names(tax_table(x))}. If 
#' \code{name} exists in both \code{sample_data} and \code{tax_table}, this will 
#' return a vector contains both columns with this \code{name}. 
#' 
#' @return Corresponding object or information for \code{name} of \code{x}.
#' 
#' @seealso \code{\link{assign}}, \code{\link{extract}}
#' @examples
#' data(oral)
#' access(oral, "sample_data")
#' seqdep(oral)
#' oral$Group
#' Ntaxa(oral)
#' Snames(oral)
#' 
#' @name access
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("access", function(x, slot_name) standardGeneric("access"))
#' @rdname access-method
#' @usage NULL
setMethod("access", signature(x="physet"), function(x, slot_name) {
  if (slot_name == "Nsample") return(nrow(x@otu_table))
  if (slot_name == "Ntaxa") return(ncol(x@otu_table))
  if (slot_name == "Nedge") return(col(x@edge_mat))
  if (slot_name == "Snames") return(rownames(x@otu_table))
  if (slot_name == "Tnames") return(colnames(x@otu_table))
  if (slot_name == "Enames") return(colnames(x@edge_mat))
  return(slot(x, slot_name))
})
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("otu_table", function(x) standardGeneric("otu_table"))
#' @rdname access-method
#' @usage NULL
setMethod("otu_table", signature(x="physet"), function(x) x@otu_table)
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("sample_data", function(x) standardGeneric("sample_data"))
#' @rdname access-method
#' @usage NULL
setMethod("sample_data", signature(x="physet"), function(x) x@sample_data)
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("tax_table", function(x) standardGeneric("tax_table"))
#' @rdname access-method
#' @usage NULL
setMethod("tax_table", signature(x="physet"), function(x) x@tax_table)
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("phy_tree", function(x) standardGeneric("phy_tree"))
#' @rdname access-method
#' @usage NULL
setMethod("phy_tree", signature(x="physet"), function(x) x@phy_tree)
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("seqdep", function(x) standardGeneric("seqdep"))
#' @rdname access-method
#' @usage NULL
setMethod("seqdep", signature(x="physet"), function(x) x@seqdep)
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("edge_com", function(x) standardGeneric("edge_com"))
#' @rdname access-method
#' @usage NULL
setMethod("edge_com", signature(x="physet"), function(x) x@edge_com)
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("edge_len", function(x) standardGeneric("edge_len"))
#' @rdname access-method
#' @usage NULL
setMethod("edge_len", signature(x="physet"), function(x) x@edge_len)
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("edge_mat", function(x) standardGeneric("edge_mat"))
#' @rdname access-method
#' @usage NULL
setMethod("edge_mat", signature(x="physet"), function(x) x@edge_mat)
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("Nsample", function(x) standardGeneric("Nsample"))
#' @rdname access-method
#' @usage NULL
setMethod("Nsample", signature(x="physet"), function(x) nrow(x@otu_table))
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("Ntaxa", function(x) standardGeneric("Ntaxa"))
#' @rdname access-method
#' @usage NULL
setMethod("Ntaxa", signature(x="physet"), function(x) ncol(x@otu_table))
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("Nedge", function(x) standardGeneric("Nedge"))
#' @rdname access-method
#' @usage NULL
setMethod("Nedge", signature(x="physet"), function(x) ncol(x@edge_mat))
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("Snames", function(x) standardGeneric("Snames"))
#' @rdname access-method
#' @usage NULL
setMethod("Snames", signature(x="physet"), function(x) rownames(x@otu_table))
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("Tnames", function(x) standardGeneric("Tnames"))
#' @rdname access-method
#' @usage NULL
setMethod("Tnames", signature(x="physet"), function(x) colnames(x@otu_table))
#' @rdname access-method
#' @usage NULL
#' @export
setGeneric("Enames", function(x) standardGeneric("Enames"))
#' @rdname access-method
#' @usage NULL
setMethod("Enames", signature(x="physet"), function(x) colnames(x@edge_mat))
#' @rdname access-method
#' @usage NULL
#' @export
setMethod("$", signature(x="physet"), function(x, name) {
  result <- NULL
  if (name %in% names(sample_data(x)))
    result <- c(result, sample_data(x)[, name])
  if (name %in% names(tax_table(x)))
    result <- c(result, tax_table(x)[, name])
  return(result)
})
################################################################################
#' @title Update or slots or attributes in physet object. 
#' @description Assign slots or attributes with new information or object. The 
#' function will automatically update related slots based on changes. 
#' 
#' @section Usage:
#' name(x) <- value
#' x$name <- value
#' 
#' @param x A physet object. 
#' @param name Available slot name are: \code{otu_table}, \code{sample_data}, 
#' \code{tax_table}, \code{phy_tree}, \code{seqdep}. Any attributes in 
#' \code{names} for \code{sample_data(x)} and \code{tax_table(x)} can be used. 
#' If \code{name} is a new attributes for \code{x} object, specify the slot 
#' by using suffix \dQuote{._AddToS_} for \code{sample_data} and suffix
#' \dQuote{._AddToT_} for \code{tax_table} (the suffix will be removed). 
#' @param value A new object with supported class for \code{name}. 
#' @return A new \code{\link{physet-class}} object with updated information. 
#' 
#' @details The assign methods for updating whole slots actually generate a 
#' new \code{\link{physet}} object by passing the original object to the 
#' \code{\link{physet-constructor}}. While assign methods for adding variables 
#' or change names will keep the original one. \cr
#' New \code{sample_data} and \code{tax_table} are actually added to the same 
#' slot in the original \code{physet}object. Only those conflict attributes 
#' are replaced. Use \dQuote{._AddToS_} and \dQuote{._AddToT_} suffix to add 
#' extra attributes. \cr
#' Assign new values to \code{edge_com}, \code{edge_len} and \code{edge_mat} 
#' manually are not recommended as they are linked with \code{otu_table} and 
#' \code{py_tree}. Provide a new \code{otu_table} or \code{phy_tree} or both 
#' automatically trigger the \code{\link{physet-constructor}} to re-build 
#' these information securely.
#' 
#' @seealso \code{\link{access}}, \code{\link{extract}}, 
#' \code{\link{physet-constructor}}
#' @examples
#' data(oral)
#' data(oral_raw)
#' otu_table(oral) <- oral_raw$rawdata[, 1:4]
#' oral$newtaxaID._AddToT_ <- paste("ID", c(1:Ntaxa(oral)), sep="")
#' 
#' @name assign
#' @rdname assign-method
#' @usage NULL
NULL
#' @rdname assign-method
#' @usage NULL
#' @export
setGeneric("otu_table<-", function(x, value) standardGeneric("otu_table<-"))
#' @rdname assign-method
#' @usage NULL
setMethod("otu_table<-", signature(x="physet", value="ANY"), 
          function(x, value) physet(physet=x, otu_table=value))
#' @rdname assign-method
#' @usage NULL
#' @export
setGeneric("sample_data<-", function(x, value) standardGeneric("sample_data<-"))
#' @rdname assign-method
#' @usage NULL
setMethod("sample_data<-", signature(x="physet", value="ANY"), 
          function(x, value) physet(physet=x, sample_data=value))
#' @rdname assign-method
#' @usage NULL
#' @export
setGeneric("tax_table<-", function(x, value) standardGeneric("tax_table<-"))
#' @rdname assign-method
#' @usage NULL
setMethod("tax_table<-", signature(x="physet", value="ANY"), 
          function(x, value) physet(physet=x, tax_table=value))
#' @rdname assign-method
#' @usage NULL
#' @export
setGeneric("phy_tree<-", function(x, value) standardGeneric("phy_tree<-"))
#' @rdname assign-method
#' @usage NULL
setMethod("phy_tree<-", signature(x="physet", value="phyloOrNULL"), 
          function(x, value) physet(physet=x, phy_tree=value))
#' @rdname assign-method
#' @usage NULL
#' @export
setGeneric("seqdep<-", function(x, value) standardGeneric("seqdep<-"))
#' @rdname assign-method
#' @usage NULL
setMethod("seqdep<-", signature(x="physet", value="numeric"), 
          function(x, value) physet(physet=x, seqdep=value))
#' @rdname assign-method
#' @usage NULL
#' @export
setGeneric("Snames<-", function(x, value) standardGeneric("Snames<-"))
#' @rdname assign-method
#' @usage NULL
setMethod("Snames<-", signature(x="physet"), function(x, value) {
  rownames(x@otu_table) <- value
  rownames(x@sample_data) <- value
  if (!is.null(x@edge_mat))
    rownames(x@edge_mat) <- value
  return(x)
})
#' @rdname assign-method
#' @usage NULL
#' @export
setGeneric("Tnames<-", function(x, value) standardGeneric("Tnames<-"))
#' @rdname assign-method
#' @usage NULL
setMethod("Tnames<-", signature(x="physet"), function(x, value) {
  names(value) <- colnames(otu_table(x))
  colnames(x@otu_table) <- value
  rownames(x@tax_table) <- value
  if (!is.null(x@phy_tree)) {
    tree <- phy_tree(x)
    repname <- value[tree$tip.label]
    x@phy_tree$tip.label[is.na(repname)] <- repname[!is.na(repname)]
    edge_name <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]
    dimnames(x@edge_com) <- list(tree$tip.label, edge_name)
    repname <- value[colnames(x@edge_mat)]
    colnames(x@edge_mat) <- ifelse(is.na(repname), 
                                   colnames(x@edge_mat), repname)
    names(x@edge_len) <- colnames(x@edge_mat)
  }
  return(x)
})
#' @rdname assign-method
#' @usage NULL
#' @export
setMethod("$<-", signature(x="physet"), function(x, name, value) {
  if (name %in% names(sample_data(x)))
    x@sample_data[, name] <- as.vector(value)
  else if (name %in% names(tax_table(x)))
    x@tax_table[, name] <- as.vector(value)
  else if (length(grep("\\._AddToS_$", name)))
    x@sample_data[, sub("\\._AddToS_$", "", name)] <- as.vector(value)
  else if (length(grep("\\._AddToT_$", name)))
    x@tax_table[, sub("\\._AddToT_$", "", name)] <- as.vector(value)
  else
    stop(paste(name, " is undefined in sample_data and tax_table.\n", 
               "Use \'", paste(name, "._AddToS_", sep=""),  
               " to update sample_data.\n", 
               "Use \'", paste(name, "._AddToT_", sep=""),  
               " to update tax_table.\n", sep=""))
  return(x)
})
################################################################################
#' @title Extract Parts of an physet object
#' @description Extract from a physet object
#' 
#' @param x A physet object. 
#' @param i,j,... Indices specifying elements to extract: samples, taxa and 
#' edges. The function works as same as extractor in \code{\link{Extract}}, 
#' but it always set \code{drop = FALSE} to maintain slot's structure.
#' @param drop This variable is fixed as \code{FALSE}.
#' @return A subset of original \code{\link{physet-class}} object.
#' 
#' @seealso \code{\link{access}}, \code{\link{assign}}
#' @examples
#' library(Matrix)
#' data(oral)
#' oral[1:10, 2:5]
#' oral[1:5, grep("^[AP]", oral$Phylum), 
#'      colMeans(edge_mat(oral)/seqdep(oral)) > 0.3]
#' @aliases extract
#' @rdname extract-method
#' @export
setMethod("[", signature(x="physet", i="ANY", j="ANY"), function(x, i, j, ..., drop = FALSE) {
  value <- list(i, j, ...)[1:3]
  names(value) <- c("sample", "taxa", "edge")
  if (!is.null(value$sample)) {
    x@otu_table <- x@otu_table[value$sample, , drop = FALSE]
    x@sample_data <- x@sample_data[value$sample, , drop = FALSE]
    x@edge_mat <- x@edge_mat[value$sample, , drop = FALSE]
    x@seqdep <- x@seqdep[value$sample]
  }
  # eval(parse(text=expr))
  if (!is.null(value$taxa)) {
    x@otu_table <- x@otu_table[, value$taxa, drop = FALSE]
    x@tax_table <- x@tax_table[value$taxa, , drop = FALSE]
  }
  if (!is.null(value$edge)) {
    x@edge_mat <- x@edge_mat[, value$edge, drop = FALSE]
    x@edge_len <- x@edge_len[value$edge, drop = FALSE]
  }
  return(x)
})
################################################################################
