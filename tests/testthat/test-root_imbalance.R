context("root_imbalance")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 64, birth = 1, death = 0)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::root_imbalance(focal_tree),
                         treestats::root_imbalance(ltab))
  bal_tree <- treestats::create_fully_balanced_tree(focal_tree)
  a1 <- treestats::root_imbalance(bal_tree)
  testthat::expect_equal(a1, 0.5)
  a2 <- treestats::root_imbalance(treestats::phylo_to_l(bal_tree))
  testthat::expect_equal(a2, 0.5)

  unbal_tree <- treestats::create_fully_unbalanced_tree(focal_tree)

  a1 <- treestats::root_imbalance(unbal_tree)
  testthat::expect_equal(a1, 63 / 64)
  a2 <- treestats::root_imbalance(treestats::phylo_to_l(unbal_tree))
  testthat::expect_equal(a2, 63 / 64)


  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::root_imbalance(focal_tree),
                         treestats::root_imbalance(ltab))

})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::root_imbalance(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::root_imbalance(list()),
    "input object has to be phylo or ltable"
  )
})
