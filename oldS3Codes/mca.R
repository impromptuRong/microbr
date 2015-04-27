#' @name mca
#' @rdname ordinate-method
#' @usage NULL
#' @import MASS
#' @export
mca <- function(x, ...) UseMethod("mca")
#' @rdname ordinate-method
#' @usage NULL
#' @export
mca.default <- function(x, ...) {
  return(MASS::mca(x, ...))
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
mca.physet <- function(x, ...) {
  return(MASS::mca(as.matrix(otu_table(x)/seqdep(x)), ...))
}
#' @rdname ordinate-method
#' @usage NULL
#' @export
setGeneric("mca")
#' @rdname ordinate-method
#' @usage NULL
setMethod("mca", signature(x = "physet"), mca.physet)
################################################################################