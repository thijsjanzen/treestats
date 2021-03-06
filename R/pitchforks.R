#' Calculate pitchforks, from the phyloTop package, a pitchfork is a clade
#' with three tips.
#' @param input_obj phylo object or ltable
#' @param normalization "none" or "tips", in which case the found number of
#' pitchforks is divided by the expected number, following
#' McKenzie & Steel 2000.
#' @return number of pitchforks
#' @export
pitchforks <- function(input_obj, normalization = "none") {
  if (inherits(input_obj, "matrix")) {
    pitch_stat <- pitchforks_ltable_cpp(input_obj)
    if (normalization == "tips") {
      pitch_stat <- pitch_stat * 3 / length(input_obj[, 1])
    }
    return(pitch_stat)
  }
  if (inherits(input_obj, "phylo")) {
    pitch_stat <- pitchforks_cpp(as.vector(t(input_obj$edge)))
    if (normalization == "tips") {
      pitch_stat <- pitch_stat * 3 / length(input_obj$tip.label)
    }
    return(pitch_stat)
  }
}
