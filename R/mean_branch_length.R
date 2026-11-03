#' Mean branch length of a tree
#' @description Mean branch length of a tree, including extinct branches.
#' @param phy phylo object or Ltable
#' @return mean branch length
#' @export
mean_branch_length <- function(phy) {

  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    return(calc_mean_branch_length_ltable_cpp(phy))
  }
  if (inherits(phy, "phylo")) {
    return(calc_mean_branch_length_cpp(as.vector(phy$edge.length)))
  }

  stop("input object has to be phylo or ltable")
}

#' Variance of branch lengths of a tree
#' @description Variance of branch lengths of a tree,
#' including extinct branches.
#' @param phy phylo object or Ltable
#' @return variance of branch lengths
#' @export
var_branch_length <- function(phy) {
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    return(stats::var(phy$edge.length))
  }

  stop("input object has to be phylo or ltable")
}

#' Variance of internal branch lengths of a tree
#' @description Variance of internal branch lengths of a tree,
#' e.g. of branches not leading to a tip.
#' @param phy phylo object or Ltable
#' @return variance of internal branch lengths
#' @export
var_branch_length_int <- function(phy) {
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    int_branches <- get_internal_branches(phy)

    return(stats::var(int_branches))
  }

  stop("input object has to be phylo or ltable")
}

#' Mean length of internal branches of a tree
#' @description Mean length of internal branches of a tree, e.g. of
#' branches not leading to a tip.
#' @param phy phylo object or Ltable
#' @return mean of internal branch lengths
#' @export
mean_branch_length_int <- function(phy) {
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    int_branches <- get_internal_branches(phy)

    return(mean(int_branches))
  }

  stop("input object has to be phylo or ltable")
}

#' Variance of external branch lengths of a tree
#' @description Variance of external branch lengths of a tree, e.g. of
#' branches leading to a tip.
#' @param phy phylo object or Ltable
#' @return variance of external branch lengths
#' @export
var_branch_length_ext <- function(phy) {
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    ext_branches <- get_external_branches(phy)
    return(stats::var(ext_branches))
  }

  stop("input object has to be phylo or ltable")
}

#' Mean length of external branch lengths of a tree
#' @description Mean length of external branch lengths of a tree, e.g. of
#' branches leading to a tip.
#' @param phy phylo object or Ltable
#' @return mean of external branch lengths
#' @export
mean_branch_length_ext <- function(phy) {
  check_tree(phy,
             require_binary = FALSE,
             require_ultrametric = FALSE,
             require_rooted = TRUE)

  if (inherits(phy, "matrix")) {
    phy <- treestats::l_to_phylo(phy, drop_extinct = FALSE)
  }
  if (inherits(phy, "phylo")) {
    ext_branches <- get_external_branches(phy)
    return(mean(ext_branches))
  }

  stop("input object has to be phylo or ltable")
}


#' @keywords internal
get_internal_branches <- function(phy) {
  root_no <- min(phy$edge[, 1])
  int_branches <- phy$edge.length[phy$edge[, 2] >= root_no]
  return(int_branches)
}

#' @keywords internal
get_external_branches <- function(phy) {
  root_no <- min(phy$edge[, 1])
  ext_branches <- phy$edge.length[phy$edge[, 2] < root_no]
  return(ext_branches)
}
