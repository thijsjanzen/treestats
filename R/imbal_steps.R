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

  ltab <- rebase_ltable(input_obj)

  num_parent1 <- length(which(ltab[, 3] < 0)) - 1
  num_parent2 <- length(which(ltab[, 3] > 0)) - 1

  target_parent <- ltab[1, 3]
  donor_parent  <- ltab[2, 3]
  if (num_parent2 > num_parent1) {
    target_parent <- ltab[2, 3]
    donor_parent <- ltab[1, 3]
  }


  to_sample_from <- which(ltab[, 2] != target_parent &
                          ltab[, 3] != donor_parent &
                          ltab[, 3] != target_parent)

  steps_taken <- length(to_sample_from)
  if (normalize == TRUE) {
    tree_size <- length(ltab[, 1])
    max_expected <- floor(tree_size - log2(tree_size) - 1)

    if (steps_taken > max_expected) {
      cat("error!")
      ax <- steps_taken / max_expected

      cat(num_parent1, "\n")
      cat(num_parent2, "\n")
      cat(ltab, "\n")
    }

    steps_taken <- steps_taken / max_expected
  }

  return(steps_taken)
}
