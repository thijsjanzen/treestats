#' Average ladder index
#' @description Calculate the avgLadder index, from the phyloTop package.
#' Higher values indicate more unbalanced trees. To calculate the average ladder
#' index, first all potential ladders in the tree are calculated. A ladder is
#' defined as a sequence of nodes where one of the daughter branches is a
#' terminal branch, resulting in a 'ladder' like pattern. The average ladder
#' index then represents the average lenght across all observed ladders in the
#' tree.
#' @param input_obj phylo object or ltable
#' @return average number of ladders
#' @export
avg_ladder <- function(input_obj) {  #nolint
  check_tree(input_obj,
             require_binary = TRUE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)
  if (inherits(input_obj, "matrix")) {
    input_obj <- treestats::l_to_phylo(input_obj, drop_extinct = FALSE)
  }
  if (inherits(input_obj, "phylo")) {
    return(avgLadder_cpp(as.vector(t(input_obj$edge))))
  }
  stop("input object has to be phylo or ltable")
}
