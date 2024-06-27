#' Four prong index
#' @description Calculate the number of 4-tip caterpillars.
#' @param input_obj phylo object or ltable
#' @return number of 4-tip caterpillars
#' @references Chindelevitch, Leonid, et al. "Network science inspires novel
#' tree shape statistics." PloS one 16.12 (2021): e0259877.
#' Rosenberg, Noah A. "The mean and variance of the numbers of r-pronged nodes
#' and r-caterpillars in Yule-generated genealogical trees."
#' Annals of Combinatorics 10 (2006): 129-146.
#' @export
four_prong <- function(input_obj) {
  check_tree(input_obj,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(input_obj, "matrix")) {
    return(calc_four_prong_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(calc_four_prong_cpp(as.vector(t(input_obj$edge))))
  }
  stop("input object has to be phylo or ltable")
}
