#' Calculate ILnumber, from the phyloTop package.
#' @description The ILnumber is the number of internal nodes with a
#' single tip child. Higher values typically indicate a tree that
#' is more unbalanced.
#' @param input_obj phylo object or ltable
#' @param normalization "none" or "tips", in which case the result is normalized
#' by dividing by N - 2, where N is the number of tips.
#' @return ILnumber
#' @export
ILnumber <- function(input_obj, normalization = "none") { # nolint
  if (inherits(input_obj, "matrix")) {
    il_stat <- ILnumber_ltable_cpp(input_obj)
    if (normalization == "tips") {
      n <- length(input_obj[, 1])
      il_stat <- il_stat / (n - 2)
    }
    return(il_stat)
  }
  if (inherits(input_obj, "phylo")) {
    il_stat <- ILnumber_cpp(as.vector(t(input_obj$edge)))
    if (normalization == "tips") {
      n <- length(input_obj$tip.label)
      il_stat <- il_stat / (n - 2)
    }
    return(il_stat)
  }
}
