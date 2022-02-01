#' Calculate the avgLadder index, from the phyloTop package
#' @param phy phylo object
#' @return average number of ladders
#' @export
avgLadder <- function(phy) { # nolint
  avgLadders <- apply_function_phy(phy, phyloTop::avgLadder) # nolint
  return(avgLadders)
}

#' Calculate number of cherries, from the phyloTop package. A cherry is a pair
#' of sister tips.
#' @param phy phylo object
#' @return number of cherries
#' @export
cherries <- function(phy) {
  num_cherries <- apply_function_phy(phy, phyloTop::cherries)
  return(num_cherries)
}

#' Calculate ILnumber, from the phyloTop package. The ILnumber is the number
#' of internal nodes with a single tip child.
#' @param phy phylo or multiPhylo object
#' @return ILnumber
#' @export
ILnumber <- function(phy) { # nolint
  ILnum <- apply_function_phy(phy, phyloTop::ILnumber) # nolint
  return(ILnum) # nolint
}

#' Calculate pitchforks, from the phyloTop package, a pitchfork is a clade
#' with three tips.
#' @param phy phylo object
#' @return number of pitchforks
#' @export
pitchforks <- function(phy) {
  num_pitchforks <- apply_function_phy(phy, phyloTop::pitchforks)
  return(num_pitchforks)
}

#' Calculates the staircase-ness measure, from the phyloTop package. The
#' staircase-ness reflects the number of subtrees that are imbalanced, e.g.
#' subtrees where the left child has more extant tips than the right child, or
#' vice versa.
#' @param phy phylo object
#' @return number of stairs
#' @references Norström, Melissa M., et al. "Phylotempo: a set of r scripts for
#' assessing and visualizing temporal clustering in genealogies inferred from
#' serially sampled viral sequences." Evolutionary Bioinformatics 8 (2012):
#' EBO-S9738.
#' @export
stairs <- function(phy) {
  num_stairs <- apply_function_phy(phy, phyloTop::stairs)
  return(num_stairs[[1]])
}
