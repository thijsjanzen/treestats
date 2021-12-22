#' calculate avgLadder index, from the phyloTop package
#' @param phy phylo or multiPhylo object
#' @return average number of ladders
#' @export
avgLadder <- function(phy) { # nolint
  avgLadders <- apply_function_phy(phy, phyloTop::avgLadder) # nolint
  return(avgLadders)
}

#' calculate number of cherries, from the phyloTop package
#' @param phy phylo or multiPhylo object
#' @return number of cherries
#' @export
cherries <- function(phy) {
  num_cherries <- apply_function_phy(phy, phyloTop::cherries)
  return(num_cherries)
}

#' calculate ILnumber, from the phyloTop package
#' @param phy phylo or multiPhylo object
#' @return ILnumber
#' @export
ILnumber <- function(phy) { # nolint
  ILnum <- apply_function_phy(phy, phyloTop::ILnumber) # nolint
  return(ILnum) # nolint
}

#' calculate pitchforks, from the phyloTop package
#' @param phy phylo or multiPhylo object
#' @return number of pitchforks
#' @export
pitchforks <- function(phy) {
  num_pitchforks <- apply_function_phy(phy, phyloTop::pitchforks)
  return(num_pitchforks)
}

#' calculate number of stairs, from the phyloTop package
#' @param phy phylo or multiPhylo object
#' @return number of stairs
#' @export
stairs <- function(phy) {
  num_stairs <- apply_function_phy(phy, phyloTop::stairs)
  return(num_stairs)
}
