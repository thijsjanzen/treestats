#' Double Cherry index
#' @description Calculate the number of double cherries, where a single cherry
#' is a node connected to two tips, and a double cherry is a node connected
#' to two cherry nodes.
#' @param input_obj phylo object or ltable
#' @return number of cherries
#' @references Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." PloS one 16.12 (2021): e0259877.
#' @export
double_cherries <- function(input_obj) {
  check_tree(input_obj,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(input_obj, "matrix")) {
    return(calc_double_cherries_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(calc_double_cherries_cpp(as.vector(t(input_obj$edge))))
  }
  stop("input object has to be phylo or ltable")
}
