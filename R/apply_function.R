#' @keywords internal
apply_function_phy <- function(input_phy, stat_function, ...) {

  if (!inherits(input_phy, "phylo"))  {
    stop("object \"phy\" is not of class \"phylo\"")
  }

  return(stat_function(input_phy, ...))
}
