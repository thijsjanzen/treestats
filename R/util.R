#' @keywords internal
check_normalization_key <- function(normalization) {
  output <- normalization
  if (output == "Yule") output <- "yule"
  if (output == "PDA")  output <- "pda"
  if (output == "Tips") output <- "tips"

  return(output)
}
