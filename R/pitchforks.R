#' Calculate pitchforks, from the phyloTop package, a pitchfork is a clade
#' with three tips.
#' @param input_obj phylo object or ltable
#' @return number of pitchforks
#' @export
pitchforks <- function(input_obj) {
  if (inherits(input_obj, "matrix")) {
    return(pitchforks_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(pitchforks_cpp(as.vector(t(input_obj$edge))))
  }
}
