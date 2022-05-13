#' Calculate the avgLadder index, from the phyloTop package
#' @param input_obj phylo object or ltable
#' @return average number of ladders
#' @export
avgLadder <- function(input_obj) { # nolint
  if (inherits(input_obj, "matrix")) {
    return(avgLadder_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(avgLadder_cpp(as.vector(t(input_obj$edge))))
  }
}
