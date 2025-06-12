#' @keywords internal
check_normalization_key <- function(normalization) {
  output <- normalization
  if (output == "Yule") output <- "yule"
  if (output == "PDA")  output <- "pda"
  if (output == "Tips") output <- "tips"

  return(output)
}

#' @keywords internal
check_binary <- function(phy) {
  # in some weird testing cases, ape::is.binary returned a vector of integers,
  # somehow this local (identical) version does not.
  n <- length(phy$tip.label)
  m <- phy$Nnode
  dgr <- tabulate(phy$edge, n + m)
  ref <- c(rep.int(1L, n), rep.int(3L, m))
  ## can use identical() as long as tabulate() returns integers
  if (ape::is.rooted(phy)) ref[n + 1L] <- 2L
  a1 <- identical(dgr, ref)
  # check that the root node is indeed binary
  a2 <- dgr[n + 1L] == 2
  return(a1 && a2)
}

#' @keywords internal
check_tree <- function(phy,
                       require_binary = FALSE,
                       require_ultrametric = FALSE,
                       require_rooted = FALSE) {

  # early exit
  if (!require_binary &&
      !require_ultrametric &&
      !require_rooted) return()


  if (inherits(phy, "phylo")) {
    if (require_binary) {
      valid <- check_binary(phy)
      if (!valid) {
        stop("Tree is non-binary, statistic not applicable")
      }
    }
    if (require_ultrametric) {
      valid <- ape::is.ultrametric(phy, tol = 1e-7, option = 1)

      if (!valid) {
        stop("Tree is not ultrametric, statistic not applicable")
      }
    }

    if (require_rooted) {
      valid <- ape::is.rooted(phy)

      if (!valid) {
        stop("Tree is not rooted, statistic not applicable")
      }
    }

  }
  if (inherits(phy, "matrix")) {
    if (require_ultrametric) {
      valid <- sum(phy[, 4] != -1)
      if (valid > 0) {
        stop("Tree is not ultrametric, statistic not applicable")
      }
    }

    if (require_binary) {
        max_num_branch_events <- max(table(phy[, 1]))
        if (max_num_branch_events > 2) {
          stop("Tree is non-binary, statistic not applicable")
        }
    }
  }
}
