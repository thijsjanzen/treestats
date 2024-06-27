context("rpanda")

test_that("usage", {
  testthat::skip_on_cran() # RPANDA uses iGraph for the Laplacian spectrum
                 # this can cause BLAS issues

  if (requireNamespace("RPANDA")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 5, birth = 1, death = 0)

    ref <- RPANDA::spectR(focal_tree)
    stat <- treestats::laplacian_spectrum(focal_tree)

    diff_eig <- sum(ref$eigenvalues - stat$eigenvalues)
    testthat::expect_equal(0, diff_eig)
    testthat::expect_equal(ref$asymmetry, stat$asymmetry)
    testthat::expect_equal(ref$peakedness, stat$peakedness, tolerance = 0.01)
    testthat::expect_equal(ref$principal_eigenvalue, stat$principal_eigenvalue)
    testthat::expect_equal(ref$eigengap, stat$eigengap)


    stat2 <- treestats::laplacian_spectrum(
      phy = treestats::phylo_to_l(focal_tree)
    )
    testthat::expect_true(all.equal(stat, stat2))
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::laplacian_spectrum(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::laplacian_spectrum(list()),
    "input object has to be phylo or ltable"
  )
})
