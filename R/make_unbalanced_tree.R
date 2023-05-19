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
#' Three methods are available: "youngest", reattaches branches in order of age,
#' starting with the branch originating from the most recent branching event
#' and working itself through the tree. "Random" picks a random branch to
#' reattach. "Terminal" also picks a random branch, but only from terminal
#' branches (e.g. branches that don't have any daughter lineages, which is
#' maximized in a fully imbalanced tree).
#' @param init_tree starting tree to work with
#' @param unbal_steps number of imbalance generating steps
#' @param group_method choice of "any" and "terminal"
#' @param selection_method choice of "random", "youngest" and "oldest"
#' @return phylo object
#' @export
#' @examples
#' simulated_tree <- ape::rphylo(n = 16, birth = 1, death = 0)
#' brts <- branching_times(simulated_tree)
#' if (requireNamespace("nodeSub")) {
#'   balanced_tree <- nodeSub::create_balanced_tree(brts)
#'   unbalanced_tree <- nodeSub::create_unbalanced_tree(brts)
#'   intermediate_tree <- make_unbalanced_tree(balanced_tree, 8)
#'   colless(balanced_tree)
#'   colless(intermediate_tree) # should be intermediate value
#'   colless(unbalanced_tree) # should be highest colless value
#' }
make_unbalanced_tree <- function(init_tree,
                                 unbal_steps,
                                 group_method = "any",
                                 selection_method = "random") {
  if (!group_method %in% c("any", "terminal")) {
    stop("group method unknown, pick from 'any' and 'terminal'")
  }
  if (!selection_method %in% c("youngest", "oldest", "random")) {
    stop("selection method unknown, pick from 'youngest', 'oldest' or 'random")
  }

  ltab_in <- treestats::phylo_to_l(init_tree)
  ltab_in <- rebase_ltable(ltab_in)
  ltab_out <- NULL

  if (group_method == "any") {
    ltab_out <- make_any(ltab_in, unbal_steps, selection_method)
  }

  if (group_method == "terminal") {
    ltab_out <- make_terminal(ltab_in, unbal_steps, selection_method)
  }

  phy_out <- treestats::l_to_phylo(ltab_out)
  return(phy_out)
}

#' @keywords internal
make_any <- function(ltab_in, unbal_steps, selection_method) {
  ltab_out <- c()
  if (selection_method == "youngest") {
    ltab_out <- make_unbalanced_tree_youngest(ltab_in, unbal_steps)
  } else if (selection_method == "oldest") {
    ltab_out <- make_unbalanced_tree_oldest(ltab_in, unbal_steps)
  } else if (selection_method == "random") {
    ltab_out <- make_unbalanced_tree_random(ltab_in, unbal_steps)
  }

  return(ltab_out)
}

#' @keywords internal
make_terminal <- function(ltab_in, unbal_steps, selection_method) {
  ltab_out <- NULL
  if (selection_method == "random") {
    ltab_out <- make_unbalanced_tree_terminal(ltab_in, unbal_steps)
  } else if (selection_method == "youngest") {
    ltab_out <- make_unbalanced_tree_terminal_youngest(ltab_in, unbal_steps)
  } else if (selection_method == "oldest") {
    ltab_out <- make_unbalanced_tree_terminal_oldest(ltab_in, unbal_steps)
  }
  return(ltab_out)
}


#' @keywords internal
get_attractor <- function(ltab) {
  attractor <- 2
  # determine attractor based on number of direct daughters,
  # more daughters = less required movements.
  num_one   <- length(which(ltab[, 2] == -1))
  num_two   <- length(which(ltab[, 2] == 2))

  if (num_one > num_two) {
    attractor <- -1
  }
  return(attractor)
}

#' @keywords internal
make_unbalanced_tree_youngest <- function(ltab,
                                          unbal_steps) {

  attractor <- get_attractor(ltab)

  to_sample_from <- which(ltab[, 2] != attractor &
                          ltab[, 3] != -1 &
                          ltab[, 3] != 2)
  steps_taken <- 0
  while (steps_taken < unbal_steps) {
    ages <- ltab[to_sample_from, 1]
    focal_step <- to_sample_from[which.min(ages)]
    ltab[focal_step, 2] <- attractor
    to_sample_from <- which(ltab[, 2] != attractor &
                            ltab[, 3] != -1 &
                            ltab[, 3] != 2)
    if (length(to_sample_from) < 1) break
    steps_taken <- steps_taken + 1
  }
  return(ltab)
}


#' @keywords internal
make_unbalanced_tree_oldest <- function(ltab,
                                          unbal_steps) {
  attractor <- get_attractor(ltab)

  to_sample_from <- which(ltab[, 2] != attractor &
                          ltab[, 3] != -1 &
                          ltab[, 3] != 2)
  steps_taken <- 0
  while (steps_taken < unbal_steps) {
    ages <- ltab[to_sample_from, 1]
    focal_step <- to_sample_from[which.max(ages)]
    ltab[focal_step, 2] <- attractor
    to_sample_from <- which(ltab[, 2] != attractor &
                            ltab[, 3] != -1 &
                            ltab[, 3] != 2)
    if (length(to_sample_from) < 1) break
    steps_taken <- steps_taken + 1
  }
  return(ltab)
}


#' @keywords internal
make_unbalanced_tree_random <- function(ltab,
                                        unbal_steps) {
  attractor <- get_attractor(ltab)

  to_sample_from <- which(ltab[, 2] != attractor &
                          ltab[, 3] != -1 &
                          ltab[, 3] != 2)

  if (length(to_sample_from) >= 1) {
    for (steps_taken in 1:unbal_steps) {
      focal_step <- DDD::sample2(to_sample_from, 1)
      if (length(to_sample_from) == 1) focal_step <- to_sample_from
      ltab[focal_step, 2] <- attractor
      to_sample_from <- which(ltab[, 2] != attractor &
                              ltab[, 3] != -1 &
                              ltab[, 3] != 2)

      if (length(to_sample_from) < 1) break
    }
  }
  return(ltab)
}


#' @keywords internal
make_unbalanced_tree_terminal <- function(ltab,
                                          unbal_steps) {
  attractor <- get_attractor(ltab)

  ltab <- cbind(ltab, 0)
  for (i in seq_along(ltab[, 1])) {
    ref <- ltab[i, 3]
    ltab[i, 5] <- length(which(ltab[, 2] == ref))
  }

  for (n in 1:unbal_steps) {
    to_sample <- which(ltab[, 5] ==  0 &
                       ltab[, 3] !=  2 &
                       ltab[, 3] != -1 &
                       ltab[, 2] != attractor)
    if (length(to_sample) < 1) {
      break
    }

    focal_spec <- DDD::sample2(to_sample, 1)
    if (length(to_sample) == 1) focal_spec <- to_sample

    parent_spec <- abs(ltab[focal_spec, 2])
    ltab[parent_spec, 5] <- ltab[parent_spec, 5] - 1
    ltab[focal_spec, 2] <- attractor
    ltab[abs(attractor), 5] <- ltab[abs(attractor), 5] + 1
  }

  return(ltab)
}

#' @keywords internal
make_unbalanced_tree_terminal_youngest <- function(ltab,   #nolint
                                          unbal_steps) {

  attractor <- get_attractor(ltab)

  ltab <- cbind(ltab, 0)
  for (i in seq_along(ltab[, 1])) {
    ref <- ltab[i, 3]
    ltab[i, 5] <- length(which(ltab[, 2] == ref))
  }

  for (n in 1:unbal_steps) {
    to_sample <- which(ltab[, 5] ==  0 &
                       ltab[, 3] !=  2 &
                       ltab[, 3] != -1 &
                       ltab[, 2] != attractor)

    ages <- ltab[to_sample, 1]
    focal_spec <- to_sample[which.min(ages)]

    parent_spec <- abs(ltab[focal_spec, 2])
    ltab[parent_spec, 5] <- ltab[parent_spec, 5] - 1
    ltab[focal_spec, 2] <- attractor
    ltab[abs(attractor), 5] <- ltab[abs(attractor), 5] + 1
    if (length(to_sample) < 1) {
      break
    }
  }

  return(ltab)
}

#' @keywords internal
make_unbalanced_tree_terminal_oldest <- function(ltab,   #nolint
                                                   unbal_steps) {

  attractor <- get_attractor(ltab)

  ltab <- cbind(ltab, 0)
  for (i in seq_along(ltab[, 1])) {
    ref <- ltab[i, 3]
    ltab[i, 5] <- length(which(ltab[, 2] == ref))
  }

  for (n in 1:unbal_steps) {
    to_sample <- which(ltab[, 5] ==  0 &
                       ltab[, 3] !=  2 &
                       ltab[, 3] != -1 &
                       ltab[, 2] != attractor)

    ages <- ltab[to_sample, 1]
    focal_spec <- to_sample[which.max(ages)]

    parent_spec <- abs(ltab[focal_spec, 2])
    ltab[parent_spec, 5] <- ltab[parent_spec, 5] - 1
    ltab[focal_spec, 2] <- attractor
    ltab[abs(attractor), 5] <- ltab[abs(attractor), 5] + 1
    if (length(to_sample) < 1) {
      break
    }
  }

  return(ltab)
}
