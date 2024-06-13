#' Create an unbalanced tree (caterpillar tree)
#' @param phy phylo object
#' @return phylo phylo object
#' @description This function takes an input phylogeny, and returns a phylogeny
#' that is a perfectly imbalanced tree (e.g. a full caterpillar tree), that
#' has the same branching times as the original input tree.
#' @export
#' @examples
#' phy <- ape::rphylo(n = 16, birth = 1, death = 0)
#' bal_tree <- treestats::create_fully_unbalanced_tree(phy)
#' treestats::colless(phy)
#' treestats::colless(bal_tree) # much higher
create_fully_unbalanced_tree <- function(phy) {

  if (!inherits(phy, "phylo")) {
    stop("This function requires a phylogeny as input")
  }

  brts <- as.vector(sort(treestats::branching_times(phy), decreasing = TRUE))

  ltab <- rbind(c(brts[1], 0, -1, -1),
                c(brts[1], -1, 2, -1))

  leftcnt <- 2
  cnt <- 2
  for (i in 2:length(brts)) {
    focal_brt <- brts[i]
    cnt <- cnt + 1
    parent <- leftcnt
    daughter <- cnt * -1
    leftcnt <- daughter
    to_add <- c(focal_brt, parent, daughter, -1)
    ltab <- rbind(ltab, to_add)
  }

  return(treestats::l_to_phylo(ltab))
}


#' Create a fully balanced tree
#' @param phy phylo object
#' @return phylo phylo object
#' @description This function takes an input phylogeny, and returns a phylogeny
#' that is most ideally balanced tree, whilst having the same branching times as
#' the original input tree. Please note that if the number of tips is not even
#' or not a power of two, the tree may not have perfect balance, but the most
#' ideal balance possible.
#' @export
#' @examples phy <- ape::rphylo(n = 16, birth = 1, death = 0)
#' bal_tree <- treestats::create_fully_balanced_tree(phy)
#' treestats::colless(phy)
#' treestats::colless(bal_tree) # much lower
create_fully_balanced_tree <- function(phy) {
  if (!inherits(phy, "phylo")) {
    stop("This function requires a phylogeny as input")
  }

  brts <- as.vector(sort(treestats::branching_times(phy), decreasing = TRUE))

  ltab <- rbind(c(brts[1], 0, -1, -1),
                c(brts[1], -1, 2, -1))

  cnt <- 3
  current_tips_left <- c(2)
  current_tips_right <- c(-1)
  brts <- brts[-1]

  while (length(brts) > 0) {
    num_to_do <- length(current_tips_left)
    for (i in 1:num_to_do) {
      if (i > length(brts)) {
        break
      }
      brt_left <- brts[i]
      parent <- current_tips_left[i]
      daughter <- cnt
      to_add <- c(brt_left, parent, daughter, -1)
      current_tips_left <- c(current_tips_left, daughter)
      ltab <- rbind(ltab, to_add)
      cnt <- cnt + 1
    }
    brts <- brts[-c(1:num_to_do)]

    num_to_do <- length(current_tips_right)

    for (i in 1:num_to_do) {
      if (i > length(brts)) {
        break
      }
      brt_right <- brts[i]
      parent <- current_tips_right[i]
      daughter <- cnt * -1
      to_add <- c(brt_right, parent, daughter, -1)
      current_tips_right <- c(current_tips_right, daughter)
      ltab <- rbind(ltab, to_add)
      cnt <- cnt + 1
    }
    brts <- brts[-c(1:num_to_do)]

  }
  return(treestats::l_to_phylo(ltab))
}
