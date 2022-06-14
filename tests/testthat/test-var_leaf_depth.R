context("var leaf depth")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1 <- treestats::var_leaf_depth(focal_tree)
  a2 <- treebalance::varLeafDepI(focal_tree)
  a2 <- a2 * (100) / (100 - 1)  # treebalance uses n, we use n - 1


  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::var_leaf_depth(focal_tree),
                         treestats::var_leaf_depth(ltab))

  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.3, fossils = TRUE)

  a1 <- treestats::var_leaf_depth(focal_tree)
  a2 <- treebalance::varLeafDepI(focal_tree)
  a2 <- a2 * (100) / (100 - 1)  # treebalance uses n, we use n - 1

  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::var_leaf_depth(focal_tree),
                         treestats::var_leaf_depth(ltab))

})
