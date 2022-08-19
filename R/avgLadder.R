#' Calculate the avgLadder index, from the phyloTop package.
#' Higher values indicate more unbalanced trees.
#' @param input_obj phylo object or ltable
#' @return average number of ladders
#' @export
avgLadder <- function(input_obj) { # nolint
  if (inherits(input_obj, "matrix")) {
    input_obj <- treestats::l_to_phylo(input_obj, drop_extinct = FALSE)
  }
  if (inherits(input_obj, "phylo")) {
    return(avgLadder_cpp(as.vector(t(input_obj$edge))))
  }
  stop("input object has to be phylo or ltable")
}
