context("blum")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  blum1 <- treestats::blum(focal_tree)
  blum_check <- treebalance::sShapeI(focal_tree, logbase = exp(1))
  testthat::expect_equal(blum1, blum_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  blum2 <- treestats::blum(ltab)
  testthat::expect_equal(blum1, blum2)
})

test_that("usage", {
  set.seed(5)
  focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)

  blum1 <- treestats::blum(focal_tree)
  blum_check <- log(2 - 1) + log(3 - 1) + log(4 - 1)
  testthat::expect_equal(blum1, blum_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  blum2 <- treestats::blum(ltab)
  testthat::expect_equal(blum1, blum2)
})
