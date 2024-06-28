#' J^1 index.
#' @description The J^1 index calculates the Shannon Entropy of a tree, where
#' at each node with two children, the Shannon Entropy is the sum of
#' p_i log_2(p_i) over the two children i, and p_i is L / (L + R), where L and
#' R represent the number of tips connected to the two daughter branches.
#' @param input_obj phylo object or ltable
#' @return j^1 index
#' @references Jeanne Lemant, Cécile Le Sueur, Veselin Manojlović, Robert Noble,
#' Robust, Universal Tree Balance Indices, Systematic Biology, Volume 71,
#' Issue 5, September 2022, Pages 1210–1224,
#' https://doi.org/10.1093/sysbio/syac027
#' @export
j_one <- function(input_obj) {

  check_tree(input_obj,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(input_obj, "matrix")) {
    return(calc_j_one_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(calc_j_one_cpp(as.vector(t(input_obj$edge))))
  }

  stop("input object has to be phylo or ltable")
}
