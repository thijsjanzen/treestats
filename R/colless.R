#' fast function using C++ to calculate the colless index of (im)balance
#' @param phy phylo object
#' @param normalization A character string equals to NULL (default) for no
#' normalization or one of "pda" or "yule".
#' @return colless index
#' @export
colless <- function(phy,
                    normalization = "none") {
  colless_index <- apply_function_phy(phy, calc_colless_cpp, normalization)
  return(colless_index)
}
