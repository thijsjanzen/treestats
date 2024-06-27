context("four prong")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)

  bal_tree <- treestats::create_fully_balanced_tree(focal_tree)

  c1 <- treestats::four_prong(bal_tree)
  testthat::expect_equal(c1, 0)

  ltab <- treestats::phylo_to_l(bal_tree)
  c2 <- treestats::four_prong(ltab)
  testthat::expect_equal(c2, 0)


  unbal_tree <- treestats::create_fully_unbalanced_tree(focal_tree)

  c1 <- treestats::four_prong(unbal_tree)
  testthat::expect_equal(c1, 1)

  ltab <- treestats::phylo_to_l(unbal_tree)
  c2 <- treestats::four_prong(ltab)
  testthat::expect_equal(c2, 1)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::four_prong(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::four_prong(list()),
    "input object has to be phylo or ltable"
  )
})
