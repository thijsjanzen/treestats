#' Calculate the avgLadder index, from the phyloTop package
#' @param input_obj phylo object or ltable
#' @return average number of ladders
#' @export
avgLadder <- function(input_obj) { # nolint
  if (inherits(input_obj, "matrix")) {
    input_obj <- treestats::l_to_phylo(input_obj)
    return(avgLadder_cpp(as.vector(t(input_obj$edge))))
  }
  if (inherits(input_obj, "phylo")) {
    return(avgLadder_cpp(as.vector(t(input_obj$edge))))
  }
}

#' Calculate number of cherries, from the phyloTop package. A cherry is a pair
#' of sister tips.
#' @param input_obj phylo object or ltable
#' @return number of cherries
#' @export
cherries <- function(input_obj) {

  if (inherits(input_obj, "matrix")) {
    return(cherries_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(cherries_cpp(as.vector(t(input_obj$edge))))
  }
}

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



#' Calculates the staircase-ness measure, from the phyloTop package. The
#' staircase-ness reflects the number of subtrees that are imbalanced, e.g.
#' subtrees where the left child has more extant tips than the right child, or
#' vice versa.
#' @param input_obj phylo object or ltable
#' @return number of stairs
#' @references Norström, Melissa M., et al. "Phylotempo: a set of r scripts for
#' assessing and visualizing temporal clustering in genealogies inferred from
#' serially sampled viral sequences." Evolutionary Bioinformatics 8 (2012):
#' EBO-S9738.
#' @export
stairs <- function(input_obj) {
  if (inherits(input_obj, "matrix")) {
    return(stairs_ltable_cpp(input_obj))
  }
  if (inherits(input_obj, "phylo")) {
    return(stairs_cpp(as.vector(t(input_obj$edge))))
  }
}
