#' calculate the spectra of eigenvalues for the modified graph Laplacian
#' of a phylogenetic tree. These are taken from the RPANDA package.
#' @param phy phylo object
#' @param meth method used to calculate the spectral density, either 'standard'
#' or 'normal'. Standard calculates the unnormalized spectral density, 'normal'
#' normalizes the spectral density to the degree matrix.
#' @return a list with the following values:
#' \itemize{
#' \item{eigenvalues} the vector of eigenvalues
#' \item{principal_eigenvalue} the largest eigenvalue
#' \item{asymmetry} the skewness of the psectral density profile
#' \item{peak_height} the larges y-axis value of the spectral density profile
#' \item{eigengap} the position of the largest difference between eigenvalues,
#' given the number of modalities in the tree.
#' }
#' @references Lewitus, E., Morlon, H., Characterizing and comparing phylogenies
#' from their Laplacian spectrum, Systematic Biology, Volume 65, Issue 3,
#' 2016, Pages 495–507, https://doi.org/10.1093/sysbio/syv116
#' @export
spectral_density <- function(phy,
                    meth = "standard") {
  if (class(phy) != "phylo") {
    stop("input has to be phylo object")
  }
  return(RPANDA::spectR(phy, meth = meth))
}
