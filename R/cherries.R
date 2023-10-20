#' Cherry index
#' @description Calculate the number of cherries, from the phyloTop package.
#' A cherry is a pair of sister tips.
#' @param input_obj phylo object or ltable
#' @param normalization "none", "yule", or "pda", the found number of
#' cherries is divided by the expected number, following
#' McKenzie & Steel 2000.
#' @return number of cherries
#' @references McKenzie, Andy, and Mike Steel. "Distributions of cherries for
#' two models of trees." Mathematical biosciences 164.1 (2000): 81-92.
#' @export
cherries <- function(input_obj, normalization = "none") {
  normalization <- check_normalization_key(normalization)

  num_cherries <- 0

  if (inherits(input_obj, "matrix")) {
    num_cherries <- cherries_ltable_cpp(input_obj)
  } else if (inherits(input_obj, "phylo")) {
    num_cherries <- cherries_cpp(as.vector(t(input_obj$edge)))
  } else {
    stop("input object has to be phylo or ltable")
  }

  if (normalization == "yule") {
    if (inherits(input_obj, "matrix")) {
      n <- length(input_obj[, 1])
    }
    if (inherits(input_obj, "phylo")) {
      n <- length(input_obj$tip.label)
    }
    expected_num_cherries <- n / 3
    num_cherries <- num_cherries / expected_num_cherries
  }

  if (normalization == "pda") {
    if (inherits(input_obj, "matrix")) {
      n <- length(input_obj[, 1])
    }
    if (inherits(input_obj, "phylo")) {
      n <- length(input_obj$tip.label)
    }
    expected_num_cherries <- n / 4
    num_cherries <- num_cherries / expected_num_cherries
  }

  return(num_cherries)
}
