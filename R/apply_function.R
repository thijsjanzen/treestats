#' @keywords internal
apply_function_ltab <- function(input_data, stat_function, ...) {

  if (is.matrix(input_data)) {
    if (is_ltable(input_data)) {
      return(stat_function(input_data, ...))
    } else {
      stop("input data has to be an ltable or phylo object")
    }
  }

  if (is.list(input_data)) {
    if (class(input_data) == "multiPhylo") {
      warning("removing extinct branches to be sure!")
      input_data <- lapply(input_data, geiger::drop.extinct)
      input_data <- lapply(input_data, DDD::phylo2L)
      output_data <-  lapply(input_data, stat_function, ...)
      return(unlist(output_data))
    }
  }

  if (class(input_data) == "phylo") {
    if (length(geiger::is.extinct(input_data)) > 0) {
      warning("removed extinct branches")
      input_data <- geiger::drop.extinct(input_data)
    }

    input_data <- DDD::phylo2L(input_data)
    return(stat_function(input_data, ...))
  }

  stop("input data has to be an ltable or phylo object")
}

#' @keywords internal
apply_function <- function(input_data, stat_function, ...) {

  if (is.matrix(input_data)) {
    if (is_ltable(input_data)) {
      warning("received ltable, converting to phylo")
      input_data <- DDD::L2phylo(input_data, dropextinct = TRUE)
    } else {
      stop("input data has to be an ltable or phylo object")
    }
  }

  if (is.list(input_data)) {
    if (class(input_data) == "multiPhylo") {
      input_data <- lapply(input_data, geiger::drop.extinct)
      output_data <-  lapply(input_data, stat_function, ...)
      return(unlist(output_data))
    }
  }

  if (class(input_data) == "phylo") {
    inptu_data <- geiger::drop.extinct(input_data)
    return(stat_function(input_data, ...))
  }

  stop("input data has to be an ltable or phylo object")
}
