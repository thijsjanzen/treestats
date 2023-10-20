#' @keywords internal
check_normalization_key <- function(normalization, expect_boolean = FALSE) {
  output <- normalization
  if (output == "Yule") output <- "yule"
  if (output == "PDA") output <- "pda"
  if (output == "Tips") output <- "tips"

  if (expect_boolean) {
    if (output %in% c("Yule", "yule", "pda", "PDA")) output <- TRUE
    if (output == "none") output <- FALSE
  }
  return(output)
}
