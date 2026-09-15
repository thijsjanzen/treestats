#' @keywords internal
check_normalization_key <- function(normalization) {
  output <- normalization
  if (output == "Yule") output <- "yule"
  if (output == "PDA")  output <- "pda"
  if (output == "Tips") output <- "tips"

  return(output)
}

#' Faster than the ape version
#' @description Fast Rcpp function to check if a tree is binary
#' @param phy phylo object
#' @return boolean, true if binary
#' @export
check_binary <- function(phy) {
  return(check_is_binary_rcpp(as.vector(t(phy$edge))))
}

#' @keywords internal
binary_check <- function(phy, require_binary = TRUE) {
  if (require_binary) {
    valid <- check_binary(phy)
    if (!valid) {
      stop("Tree is non-binary, statistic not applicable")
    }
  }
}

#' @keywords internal
ultrametric_check <- function(phy, require_ultrametric = TRUE) {
  if (require_ultrametric) {
    valid <- ape::is.ultrametric(phy, tol = 1e-7, option = 1)

    if (!valid) {
      stop("Tree is not ultrametric, statistic not applicable")
    }
  }
}

#' @keywords internal
rooted_check <- function(phy, require_rooted = TRUE) {
  if (require_rooted) {
    valid <- ape::is.rooted(phy)

    if (!valid) {
      stop("Tree is not rooted, statistic not applicable")
    }
  }
}

binary_check_ltable <- function(phy, require_binary) {
  if (require_binary) {
    max_num_branch_events <- max(table(phy[, 1]))
    if (max_num_branch_events > 2) {
      stop("Tree is non-binary, statistic not applicable")
    }
  }
}

ultrametric_check_ltable <- function(phy, require_ultrametric) {
  if (require_ultrametric) {
    valid <- sum(phy[, 4] != -1)
    if (valid > 0) {
      stop("Tree is not ultrametric, statistic not applicable")
    }
  }
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
    binary_check(phy, require_binary)
    ultrametric_check(phy, require_ultrametric)
    rooted_check(phy, require_rooted)
  }
  if (inherits(phy, "matrix")) {
    ultrametric_check_ltable(phy, require_ultrametric)
    binary_check_ltable(phy, require_binary)
  }
}
