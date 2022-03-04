#' @keywords internal
#' @rawNamespace import(Rcpp)
#' @rawNamespace import(nloptr)
#' @rawNamespace useDynLib(treestats)
apply_function_phy <- function(input_phy, stat_function, ...) {

  if (!inherits(input_phy, "phylo"))  {
    stop("object \"phy\" is not of class \"phylo\"")
  }

  return(stat_function(input_phy, ...))
}

#' @keywords internal
apply_function_phy_ltable <- function(input_obj,
                                      phy_function,
                                      ltable_function,
                                      only_extant = TRUE,
                                      ...) {

  if (inherits(input_obj, "phylo")) {
    if (only_extant && !ape::is.ultrametric(input_obj)) {
      stop("can only calculate statistic for ultrametric tree")
    }

    return(phy_function(input_obj, ...))
  }

  if (inherits(input_obj, "matrix")) {
    return(ltable_function(input_obj, ...))
  }

  stop("input object has to be phylo or ltable")
}
