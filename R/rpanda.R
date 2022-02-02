#' @keywords internal
integr <- function(x, f)
{
  # get lengths of var and integrand
  n = length(x)

  # trapezoidal integration
  integral = 0.5 * sum((x[2:n] - x[1:(n - 1)]) * (f[2:n] + f[1:(n - 1)]))

  # print definite integral
  return(integral)
}

#' @keywords internal
skewness <- function(x, na.rm = FALSE) {
  if (is.matrix(x))
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x))
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

#' function to calculate the laplacian spectrum, from RPANDA
#' @param phy phy
#' @return list
#' @export
calc_lapl_spectrum <- function(phy) {

  lapl_mat <- prep_lapl_spec(phy)

  e = eigen(lapl_mat, symmetric = TRUE, only.values = TRUE)

  x = subset(e$values, e$values >= 1)
  d = density(log(x))
  dsc = d$y / (integr(d$x, d$y))
  principal_eigenvalue <- max(x)
  skewness <- skewness(x)
  peak_height <- max(dsc)
  gaps <- abs(diff(x))
  gapMat <- as.matrix(gaps)
  modalities <- c(1:length(gapMat))
  gapMatCol <- cbind(modalities, gapMat)
  eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[,2]))
  res <- list(eigenvalues = x,
              principal_eigenvalue = principal_eigenvalue,
              asymmetry = skewness,
              peakedness = peak_height,
              eigengap = eigenGap[,1])
  return(res)
}
