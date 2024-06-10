context("double cherries")

test_that("usage", {
  if (requireNamespace("nodeSub")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)
    brts <- treestats::branching_times(focal_tree)
    bal_tree <- nodeSub::create_balanced_tree(brts)

    c1 <- treestats::double_cherries(bal_tree)
    testthat::expect_equal(c1, 1)

    ltab <- treestats::phylo_to_l(bal_tree)
    c2 <- treestats::double_cherries(ltab)
    testthat::expect_equal(c2, 1)

    focal_tree <- ape::rphylo(n = 8, birth = 1, death = 0)
    brts <- treestats::branching_times(focal_tree)
    bal_tree <- nodeSub::create_balanced_tree(brts)

    c1 <- treestats::double_cherries(bal_tree)
    testthat::expect_equal(c1, 2)

    ltab <- treestats::phylo_to_l(bal_tree)
    c2 <- treestats::double_cherries(ltab)
    testthat::expect_equal(c2, 2)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::double_cherries(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::double_cherries(list()),
    "input object has to be phylo or ltable"
  )
})
