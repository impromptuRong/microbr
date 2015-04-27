################################################################################
#' @title Dong Lab Microbiome Analysis Pipeline.
#' 
#' @description
#' The microbr package is designed is a general purpose wrapper and pipeline 
#' for lots of routine microbiome and metagenomics analysis. The pipeline takes 
#' meta data, abundance table, taxonomy information and phylogenetic tree from 
#' upstream sequencing analysis as input and simplify downstream data exploring. 
#' \tabular{ll}{
#' Diversity analysis: \tab alpha/beta-diversity, richness, eveness \cr
#' Univariate analysis: \tab (non)-parametric, multi-group, p-value correction \cr
#' Multivariate analysis: \tab ordination analysis, distance analysis \cr
#' Enterotype analysis: \tab enterotype clustering \cr
#' Machine learning: \tab classification, regression, feature selection, clustering \cr
#' }
#' 
#' @section Important:
#' This package is used in Dong's Lab only or for personal learning purpose. 
#' Similar packages (like \code{\link[phyloseq]{phyloseq}}) have been developed and well 
#' maintained in Github and Bioconductor. Some codes are copied from or highly 
#' similar to other contributions. Please \emph{CONTACT the authors} for 
#' permission if codes will be published or be used for other purpose. 
#' 
#' @details
#' \tabular{ll}{
#' Package: \tab microbr\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2015-03-26\cr
#' License: \tab UNT\cr
#' }
#' 
# @import methods
#' @docType package
#' @name microbr-package
#' @keywords package
#' @author 
#' Developer: Ruichen Rong <\email{imprompturong@@gmail.com}> \cr
#' Maintainer: Kashi <\email{}>, Huaiying <\email{}>, Yuqing<\email{}>
#' @references \url{Too Many}
#' @seealso \code{\link[phyloseq]{phyloseq-package}}
#' @examples
#' ## ~~ simple examples of the most important functions ~~
NULL
################################################################################
#' @title Oral Microbiome Dataset - raw input
#'
#' @description A microbiome dataset containing the phylum count table, sample
#' information, phyla information and phylogenetic tree. 
#' The variables are as follows:
#'
#' \itemize{
#'   \item rawdata. The phylum count table for Oral data. 
#'   \item metainfo. The meta information for samples.
#'   \item taxonomy. The taxa information for phyla.
#'   \item rawtree. The phylogenetic tree on phylum level.
#' }
#'
#' @format A list with 4 slots for oral microbiome. 
#' @source \url{http://aaa/}
#' @name oral_raw
NULL
################################################################################
#' @title Oral Microbiome Dataset - physet object
#'
#' @description The physet object constructed from \code{data(oral_raw)}. It is
#' a sample physet object for testing purpose. See \code{\link{physet-class}} 
#' and \code{\link{physet}} for details.
#'
#' @format A physet object with otu_table, sample_data, tax_table, phy_tree. 
#' @source \url{http://aaa/}
#' @name oral
NULL
################################################################################
if(getRversion() >= "2.15.1")
  utils::globalVariables(names = c("MICROBR_GLOBAL_x", "MICROBR_GLOBAL_i", 
                                   "MICROBR_GLOBAL_j"), package = "microbr")

# if(getRversion() >= "3.1.0") {
#   utils::suppressForeignCheck("%dopar%")
#   utils::suppressForeignCheck("%foreach%")
# }