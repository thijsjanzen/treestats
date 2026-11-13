context("blum")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    standard_test(treestats::blum, treebalance::sShapeI)
  }
})

test_that("usage", {
  set.seed(5)
  focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)
  blum1 <- treestats::blum(focal_tree)
  blum_check <- log2(2 - 1) + log2(3 - 1) + log2(4 - 1)
  testthat::expect_equal(blum1, blum_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  blum2 <- treestats::blum(ltab)
  testthat::expect_equal(blum1, blum2)

  sshape1 <- treestats::sshape(focal_tree)
  testthat::expect_equal(sshape1, blum1)
  sshape2 <- treestats::sshape(ltab)
  testthat::expect_equal(sshape2, blum1)
})

test_that("polytomies", {
  if (requireNamespace("treebalance")) {
    test_polytomies(treestats::blum, treebalance::sShapeI)
  }
})


test_that("wrong_object", {
  testthat::expect_error(
    treestats::blum(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::blum(list()),
    "input object has to be phylo or ltable"
  )
})
