% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ordinate-methods.R
\docType{methods}
\name{ordinate}
\alias{dca}
\alias{ordinate}
\alias{ordinate,physet-method}
\alias{ordinate.default}
\alias{ordinate.physet}
\alias{phycap}
\alias{phycap,physet,formula-method}
\alias{phycap.default}
\alias{phycap.dist}
\alias{phycap.physet}
\alias{phycca}
\alias{phycca,physet,formula-method}
\alias{phycca.default}
\alias{phycca.physet}
\alias{phydca}
\alias{phydca,physet-method}
\alias{phydca.default}
\alias{phydca.physet}
\alias{phymca}
\alias{phymca,physet-method}
\alias{phymca.default}
\alias{phymca.physet}
\alias{phymds}
\alias{phymds,physet-method}
\alias{phymds.default}
\alias{phymds.dist}
\alias{phymds.physet}
\alias{phypcoa}
\alias{phypcoa,physet-method}
\alias{phypcoa.default}
\alias{phypcoa.dist}
\alias{phypcoa.physet}
\alias{phyrda}
\alias{phyrda,physet,formula-method}
\alias{phyrda.default}
\alias{phyrda.physet}
\title{Ordination Analysis for Communities}
\arguments{
\item{x}{A \code{matrix}, \code{data.frame} or any object inherits these
methods, like \code{\link{physet-class}}. \code{x} can also be a \code{dist}
matrix if ordination method use dissimilarity matrix.}

\item{...}{Additional parameters pass to analysis method or dfun.}

\item{apply}{The ordination analysis method. Current supported methods: CA,
CCA, DCA, MCA, RDA, PCA, CAP, (N)MDS, PCoA, DPCoA.}

\item{formula}{a model \code{formula} with the left hand side gives the
community variables and the right hand side gives the constrains. The LHS
part is matched to \code{x}, \code{comm} or \code{data} and The RHS part is
matched to \code{x} or \code{data}. \dQuote{otu_table, edge_mat} can be
used for \code{physet-class}. Conditioning variables can be given within a
function \code{Condition}. See \code{\link[vegan]{cca}} for details.}

\item{data}{A data frame contains variables on the RHS of the formula.
If \code{data} is missing, the function will seek RHS in \code{x}.}

\item{comm}{A (extra) community data used for calculate taxa scores based on
LHS of the formula. If \code{comm} is missing, the function will seek LHS in
\code{x}, then \code{data} if \code{x} is a \code{dist}. \code{comm} is not
used to calculate distance but only used for taxa scores.}

\item{dfun}{The function used to calculate disimilarity matrix for analysis.
Default method is \code{\link{phydist}} for \code{physet-class} and
\code{\link{vegdist}} for default functions. The function parses \code{x}
and LHS of \code{formula} to \code{dfun} if \code{x} is not a matrix.}

\item{method}{The distance/dissimilarity method pass to \code{dfun}.}
}
\value{
An ordination object. The class of the returned object depends upon
the ordination method, the function and package that perform it. Also the
object has a \code{data.frame} under \dQuote{plotdata} slot for ploting.
}
\description{
This is a wrapper for supported ordination analysis.
}
\details{
[Partial] [Constrained] Correspondence Analysis and Redundancy Analysis \cr
[Partial] Constrained Analysis of Principal Coordinates \cr
Nonmetric Multidimensional Scaling \cr
Detrended Correspondence Analysis \cr
Multiple Correspondence Analysis \cr
}
\section{Usage}{

\code{ordinate(x, apply, formula, dfun = vegdist, method = "bray", ...)} \cr
\code{phycca(x, formula, data, ...)} \cr
\code{phyrda(x, formula, data, ...)} \cr
\code{phycap(x, formula, data, comm, method, ...)} \cr
\code{phymds(x, formula, comm, method, ...)} \cr
\code{phydca(x, ...)} \cr
\code{phymca(x, ...)} \cr
}
\examples{
data(oral)
ordinate(oral, "CCA", otu_table ~ Group)
ordinate(oral, "PCoA", method = "unifrac.uw")
ordinate(oral, "CAP", edge_mat ~ Group, distance = "unifrac.w.un")
}
\seealso{
\code{\link[phyloseq]{ordinate}}, \code{\link[vegan]{decorana}},
\code{\link[vegan]{cca}}/\code{\link[vegan]{rda}}, \code{\link[MASS]{mca}},
\code{\link[vegan]{metaMDS}}, \code{\link[vegan]{capscale}}.
}

