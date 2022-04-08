#' performs all treebalance measures
#' @description This function can be used as a fast-track power user method to
#' quickly assess all balance measures provided in the treebalance package.
#' @param phylo phylo object
#' @return crown age
#' @export
calc_treebalance_stats <- function(phylo) {

  stats <- list()
  stats$B1           <- treebalance::B1I(phylo)
  stats$B2           <- treebalance::B2I(phylo)
  stats$aPP          <- treebalance::areaPerPairI(phylo)
  stats$aLD          <- treebalance::avgLeafDepI(phylo)
#  stats$cherry       <- treebalance::cherryI(phylo)  # this is cherries
 # stats$colless      <- treebalance::collessI(phylo) # this is Colless
#  stats$colPlaLab    <- treebalance::colPlaLab(phylo, method = "binary")
#  stats$furnas       <- treebalance::furnasI(phylo)  # only yields unique numbers, not really balancy
  stats$Ibased       <- treebalance::IbasedI(phylo)
  stats$ewColless    <- treebalance::ewCollessI(phylo)
  stats$maxDelW      <- treebalance::maxDelW(phylo)
  stats$maxDepth     <- treebalance::maxDepth(phylo)
  stats$maxWidth     <- treebalance::maxWidth(phylo)
  stats$rogers       <- treebalance::rogersI(phylo)
  stats$rquarted     <- treebalance::rQuartetI(phylo)
 # stats$sackin       <- treebalance::sackinI(phylo)  # this is Sackin
 # stats$sshape1      <- treebalance::sShapeI(phylo)  # this is Blum
#  stats$stairs1      <- treebalance::stairs1(phylo)  # this is stairs
 # stats$stairs2      <- treebalance::stairs2(phylo)  # this is stairs2
  stats$symNodes     <- treebalance::symNodesI(phylo)
  stats$totCoph      <- treebalance::totCophI(phylo)
  stats$varLeafDepth <- treebalance::varLeafDepI(phylo)

  return(stats)
}
