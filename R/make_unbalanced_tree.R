#' this function increasingly increases the imbalance of a tree
#' @description the goal of this function is to increasingly imbalance a tree,
#' by changing the topology, one move at a time. It does so by re-attaching
#' terminal branches to the root lineage, through the ltable. In effect, this
#' causes the tree to become increasingly caterpillarlike. When started with
#' a balanced tree, this allows for exploring the gradient between a fully
#' balanced tree, and a fully unbalanced tree.
#' Please note that the algorithm will try to increase imbalance, until a fully
#' caterpillar like tree is reached, which may occur before unbal_steps is
#' reached.
#' @param init_tree starting tree to work with
#' @param unbal_steps number of imbalance generating steps
#' @return phylo object
#' @export
#' @examples
#' simulated_tree <- ape::rphylo(n = 16, birth = 1, death = 0)
#' brts <- branching_times(simulated_tree)
#' balanced_tree <- nodeSub::create_balanced_tree(brts)
#' unbalanced_tree <- nodeSub::create_unbalanced_tree(brts)
#' intermediate_tree <- make_unbalanced_tree(balanced_tree, 8)
#' colless(balanced_tree)
#' colless(intermediate_tree) # should be intermediate value
#' colless(unbalanced_tree) # should be highest colless value
make_unbalanced_tree <- function(init_tree, unbal_steps) {
  ltab <- treestats::phylo_to_l(init_tree)

  # first, we move all daughter branches from -1
  to_sample_from <- which(ltab[, 2] == -1 & ltab[, 3] != 2)
  steps_taken <- 0
  while (TRUE) {
    focal_step <- sample(to_sample_from, 1)
    if (length(to_sample_from) == 1) focal_step <- to_sample_from
    ltab[focal_step, 2] <- 2
    to_sample_from <- which(ltab[, 2] == -1 & ltab[, 3] != 2)
    steps_taken <- steps_taken + 1
    if (steps_taken >= unbal_steps) break
    if (length(to_sample_from) < 1) break
  }

  if (steps_taken >= unbal_steps) {
    output_phy <- treestats::l_to_phylo(ltab)
    return(output_phy)
  }

  # now, we start attaching directly to 2
  to_sample_from <- which(ltab[, 2] != 2 & ltab[, 3] != -1 & ltab[, 3] != 2)

  while (steps_taken < unbal_steps) {
    ages <- ltab[to_sample_from, 1]
    focal_step <- to_sample_from[which.min(ages)]
    ltab[focal_step, 2] <- 2
    to_sample_from <- which(ltab[, 2] != 2 & ltab[, 3] != -1 & ltab[, 3] != 2)
    if (length(to_sample_from) < 1) break
    steps_taken <- steps_taken + 1
  }

  output_phy <- treestats::l_to_phylo(ltab)
  return(output_phy)
}
