% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/internal-methods.R
\docType{methods}
\name{.m2df}
\alias{.m2df}
\alias{.m2df,physet-method}
\alias{.m2df.default}
\alias{.m2df.matrix}
\alias{.m2df.physet}
\title{Convert Objects to Data.frame}
\usage{
.m2df(x, ...)

.m2df.default(x, ...)

.m2df.matrix(x, i, j, threshold = NULL, ...)

\S4method{.m2df}{physet}(x, i, j, threshold = NULL, ...)
}
\arguments{
\item{x}{A \code{list} or a object inherits this method: \code{physet-class}.}

\item{i,j,...}{The indices specifying elements in object. Indices can be
+/- \code{numeric}, \code{character} and \code{logical} \code{vectors} or
empty (missing).}

\item{threshold}{A \code{character} string indicates threshold like:
\code{"> 0"}, \code{"== 10"}, \code{"<= 0.4"}.}
}
\value{
A \code{data.frame}, each line indicates one record plus its row
and column information in x.
}
\description{
Internal function to convert matrix, data.frame, physet objects
to long-shape data.frame (sparse-Matrix format) for ploting and testing.
}
\examples{
data(oral)
.m2df(as.matrix(otu_table(oral)), 1:20, -c(1:10), "> 100")
.m2df(oral, 1:20, -c(1:10), "> 100")
}
\seealso{
\code{\link{reshape}}
}
\keyword{internal}

