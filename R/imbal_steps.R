#' imbalance steps index
#' @description Calculates the number of moves required to transform the focal
#' tree into a fully imbalanced (caterpillar) tree. Higher value indicates a
#' more balanced tree.
#' @param input_obj phylo object or ltable
#' @param normalize if true, the number of steps taken is divided by the tree
#' size
#' @return required number of moves
#' @export
imbalance_steps <- function(input_obj,
                            normalize = FALSE) {

  if (inherits(input_obj, "phylo")) {
    input_obj <- treestats::phylo_to_l(input_obj)
  }

  if (!inherits(input_obj, "matrix")) {
    stop("input object has to be phylo or ltable")
  }

  ltab <- input_obj

  parent1 <- length(which(ltab[, 2] == -1)) - 1
  parent2 <- length(which(ltab[, 2] == 2))

  target_parent <- -1
  donor_parent <- 2
  if (parent2 > parent1) {
    target_parent <- 2
    donor_parent <- -1
  }

  to_sample_from <- which(ltab[, 2] != target_parent &
                            ltab[, 3] != -1 & # not the root
                            ltab[, 3] != donor_parent)

  steps_taken <- length(to_sample_from)
  if (normalize == TRUE) {
    tree_size <- length(ltab[, 1])
    max_expected <- tree_size - log2(tree_size) - 1
    steps_taken <- steps_taken / max_expected
  }

  return(steps_taken)
}
