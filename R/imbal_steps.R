#' imbalance steps index
#' @description Calculates the number of moves required to transform the focal
#' tree into a fully imbalanced (caterpillar) tree. Higher value indicates a
#' more balanced tree.
#' @param input_obj phylo object or ltable
#' @param normalization if true, the number of steps taken is normalized by tree
#' size, by dividing by the maximum number of moves required to move from a
#' fully balanced to a fully imbalanced tree, which is N - log2(N) - 1, where
#' N is the number of extant tips.
#' @return required number of moves
#' @export
imbalance_steps <- function(input_obj,
                            normalization = FALSE) {
  normalization <- check_normalization_key(normalization, TRUE)

  if (inherits(input_obj, "phylo")) {
    input_obj <- treestats::phylo_to_l(input_obj)
  }

  if (!inherits(input_obj, "matrix")) {
    stop("input object has to be phylo or ltable")
  }
  if (length(input_obj[, 1]) == 2) {
    message("Balance of a tree with two tips is undefined")
    return(NA)
  }

  ltab <- rebase_ltable(input_obj)

  attractor <- get_attractor(ltab)
  #cat(attractor, "\n")

  to_sample_from <- which(ltab[, 2] != attractor &
                          ltab[, 3] != -1 &
                          ltab[, 3] != 2)

  steps_taken <- length(to_sample_from)
  if (normalization == TRUE) {
    tree_size <- length(ltab[, 1])
    max_expected <- tree_size - ceiling(log2(tree_size)) - 1
    steps_taken <- steps_taken / max_expected
  }

  return(steps_taken)
}
