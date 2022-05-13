#' Calculate ILnumber, from the phyloTop package. The ILnumber is the number
#' of internal nodes with a single tip child.
#' @param input_obj phylo object or ltable
#' @return ILnumber
#' @export
ILnumber <- function(input_obj) { # nolint
  if (inherits(input_obj, "matrix")) {
    return(ILnumber_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(ILnumber_cpp(as.vector(t(input_obj$edge))))
  }
}
