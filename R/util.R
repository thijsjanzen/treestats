
#' function to determine of object is ltable
#' this function is kind of crude, as it only checks some aspects that make an
#' ltable
#' @keywords internal
#' @rawNamespace import(Rcpp)
#' @rawNamespace import(nloptr)
#' @rawNamespace useDynLib(treestats)
is_ltable <- function(x) {
  if (!is.matrix(x)) return(FALSE)
  if (ncol(x) != 4) return(FALSE)
  if (length(which(x[, 4] == -1)) < 1) return(FALSE)
  if (x[1,2] != 0) return(FALSE)
  return(TRUE)
}
