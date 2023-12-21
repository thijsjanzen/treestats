#' a function to modify an ltable, such that the longest path in the phylogeny
#' is a crown lineage.
#' @param ltable ltable
#' @return modified ltable
#' @export
rebase_ltable <- function(ltable) {
  if (length(ltable[, 1]) == 2) return(ltable)
  prev_main_attractor <- c()
  while (TRUE) {
    res <- swap_deepest(ltable)
    ltable <- res$ltab
    prev_main_attractor <- c(prev_main_attractor,
                             res$main_attractor)
    if (length(prev_main_attractor) > 5) {
      end <- length(prev_main_attractor)

      prev_main_attractor <- prev_main_attractor[(end - 5):end]
    }

    if (stats::sd(prev_main_attractor) == 0 &&
        length(prev_main_attractor) > 3) {
      stop("Stuck in endless loop, possibly due to polytomies")
    }

    if (res$stop) break
  }

  new_ltable <- renumber_ltable(ltable)

  return(new_ltable)
}

#' @keywords internal
renumber_ltable <- function(ltab) {
  temp_new_ltab <- ltab

  for (i in  seq_along(temp_new_ltab[, 1])) {
    current_label <- ltab[i, 3]
    if (abs(current_label) != i) {
      new_label <- i * sign(ltab[i, 3])

      temp_new_ltab[i, 3] <- new_label
      daughters <- which(ltab[, 2] == current_label &
                           ltab[, 1] <= ltab[i, 1])
      if (length(daughters) > 0) {
        temp_new_ltab[daughters, 2] <- i * sign(ltab[i, 3])
      }

      other_instances <- which(ltab[, 3] == i &
                                 ltab[, 1] < ltab[i, 1])
      if (length(other_instances) > 0) {
        temp_new_ltab[other_instances, 3] <- current_label
      }
    }
  }
  return(temp_new_ltab)
}

#' @keywords internal
swap_deepest <- function(ltab) {
  depths <- rep(0, length(ltab[, 1]))
  depths[c(1, 2)] <- 1
  for (i in 3:length(ltab[, 1])) {
    parent_index <- abs(ltab[i, 2])
    depths[parent_index] <- depths[parent_index] + 1
    depths[i] <- depths[parent_index]
  }
  main_attractor <- which.max(depths)
  focal_index <- which(abs(ltab[, 3]) == main_attractor)
  main_attractor <- ltab[focal_index, 3]

  new_ltab <- ltab
  finalized <- FALSE

  if (!focal_index %in% c(1, 2)) {
    parent <- ltab[focal_index, 2]
    index_parent <- which(ltab[, 3] == parent)
    # switch them around, label wise:
    new_ltab[focal_index, 3] <- parent
    new_ltab[focal_index, 2] <- main_attractor
    new_ltab[index_parent, 3] <- main_attractor
    daughters <- which(new_ltab[, 2] == parent &
                         new_ltab[, 1] > new_ltab[focal_index, 1])
    if (length(daughters) > 0) {
      new_ltab[daughters, 2] <- main_attractor
    }

    new_ltab <- renumber_ltable(new_ltab)

  } else {
    finalized <- TRUE
  }

  return(list("ltab" = new_ltab,
              "stop" = finalized,
              "main_attractor" = main_attractor))
}
