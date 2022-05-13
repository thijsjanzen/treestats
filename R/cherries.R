#' Calculate number of cherries, from the phyloTop package. A cherry is a pair
#' of sister tips.
#' @param input_obj phylo object or ltable
#' @return number of cherries
#' @export
cherries <- function(input_obj) {

  if (inherits(input_obj, "matrix")) {
    return(cherries_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(cherries_cpp(as.vector(t(input_obj$edge))))
  }
}
