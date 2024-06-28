context("j_one")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 128, birth = 1, death = 0)

  bal_tree <- treestats::create_fully_balanced_tree(focal_tree)
  unbal_tree <- treestats::create_fully_unbalanced_tree(focal_tree)

  j_one_1 <- treestats::j_one(bal_tree)
  j_one_2 <- treestats::j_one(unbal_tree)

  testthat::expect_equal(j_one_1, 1)
  testthat::expect_equal(j_one_2, 0.1085403, tolerance = 0.001)

  j_one_1_l <- treestats::j_one(treestats::phylo_to_l(bal_tree))
  j_one_2_l <- treestats::j_one(treestats::phylo_to_l(unbal_tree))

  testthat::expect_equal(j_one_1_l, 1)
  testthat::expect_equal(j_one_2_l, 0.1085403, tolerance = 0.001)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::j_one(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::j_one(list()),
    "input object has to be phylo or ltable"
  )
})
