#' imbalance steps index
#' @description Calculates the number of moves required to transform the focal
#' tree into a fully imbalanced (caterpillar) tree. Higher value indicates a
#' more balanced tree.
#' @param input_obj phylo object or ltable
#' @param normalize if true, the number of steps taken is divided by the tree
#' size
#' @return required number of moves
#' @export
imbal_steps <- function(input_obj,
                        normalize = FALSE) {

  if (inherits(input_obj, "phylo")) {
    input_obj <- treestats::phylo_to_l(input_obj)
  }

  if (!inherits(input_obj, "matrix")) {
    stop("input object has to be phylo or ltable")
  }

  ltab <- input_obj

  to_sample_from <- which(ltab[, 2] != 2 &
                            ltab[, 3] != -1 &
                            ltab[, 3] != 2)

  tree_size <- length(ltab[, 1])

  expected_max_depth <- tree_size - 1

  steps_taken <- 0
  while (treestats::max_depth(ltab) < expected_max_depth) {
    ages <- ltab[to_sample_from, 1]
    focal_step <- to_sample_from[which.min(ages)]
    ltab[focal_step, 2] <- 2
    to_sample_from <- which(ltab[, 2] != 2 & ltab[, 3] != -1 & ltab[, 3] != 2)
    if (length(to_sample_from) < 1) break
    steps_taken <- steps_taken + 1
  }

  if (normalize == TRUE) {
    max_num_steps <- tree_size #- log(tree_size, 2) - 1
    steps_taken <- steps_taken / max_num_steps
  }

  return(steps_taken)
}
