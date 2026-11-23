context("max_clade_basis")

test_that("usage", {
  ref <- function(phy) {
    n <- ape::Ntip(phy)
    k <- ape::node.depth(phy, method = 1)[-(1:(n + 1))]
    sum(1 / (k * (k + 1)))
  }

  standard_test(treestats::max_clade_basis, ref, tol = 1e-4)
})

test_that("yule", {
  test_yule(treestats::max_clade_basis)
})

test_that("wrong_object", {
  test_abuse(treestats::max_clade_basis)
})
