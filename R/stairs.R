#' Stairs index
#' @description Calculates the staircase-ness measure, from the phyloTop
#' package. The staircase-ness reflects the number of subtrees that are
#' imbalanced, e.g. subtrees where the left child has more extant tips than the
#' right child, or vice versa.
#' @param input_obj phylo object or ltable
#' @return number of stairs
#' @references Norström, Melissa M., et al. "Phylotempo: a set of r scripts for
#' assessing and visualizing temporal clustering in genealogies inferred from
#' serially sampled viral sequences." Evolutionary Bioinformatics 8 (2012):
#' EBO-S9738.
#' @export
stairs <- function(input_obj) {
  check_tree(input_obj,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(input_obj, "matrix")) {
    return(stairs_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(stairs_cpp(as.vector(t(input_obj$edge))))
  }
  stop("input object has to be phylo or ltable")
}

#' Stairs2 index
#' @description Calculates the stairs2 measure, from the phyloTop package. The
#' stairs2 reflects the imbalance at each node, where it represents the average
#' across measure at each node, the measure being min(l, r) / max(l, r), where
#' l and r reflect the number of tips connected at the left (l) and right (r)
#' daughter.
#' @param input_obj phylo object or ltable
#' @return number of stairs
#' @references Norström, Melissa M., et al. "Phylotempo: a set of r scripts for
#' assessing and visualizing temporal clustering in genealogies inferred from
#' serially sampled viral sequences." Evolutionary Bioinformatics 8 (2012):
#' EBO-S9738.
#' @export
stairs2 <- function(input_obj) {
  check_tree(input_obj,
             require_binary = TRUE,
             require_ultrametric = FALSE)

  if (inherits(input_obj, "matrix")) {
    return(stairs2_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(stairs2_cpp(as.vector(t(input_obj$edge))))
  }
  stop("input object has to be phylo or ltable")
}
