#' Imbalance steps index
#' @description Calculates the number of moves required to transform the focal
#' tree into a fully imbalanced (caterpillar) tree. Higher value indicates a
#' more balanced tree.
#' @param input_obj phylo object or ltable
#' @param normalization if true, the number of steps taken is normalized by tree
#' size, by dividing by the maximum number of moves required to move from a
#' fully balanced to a fully imbalanced tree, which is \eqn{N - log2(N) - 1},
#' where N is the number of extant tips.
#' @return required number of moves
#' @export
imbalance_steps <- function(input_obj,
                            normalization = FALSE) {

  check_tree(input_obj,
             require_binary = TRUE,
             require_ultrametric = TRUE,
             require_rooted = TRUE)

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

  steps_taken <- imbalance_steps_cpp(input_obj,
                                     normalization)
  return(steps_taken)
}
