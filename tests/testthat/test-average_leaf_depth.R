context("average leaf depth")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  ald <- treestats::average_leaf_depth(focal_tree)
  ald_check <- treebalance::avgLeafDepI(focal_tree)
  testthat::expect_equal(ald, ald_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::average_leaf_depth(focal_tree),
                         treestats::average_leaf_depth(ltab))


  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  ald <- treestats::average_leaf_depth(focal_tree)
  ald_check <- treebalance::avgLeafDepI(focal_tree)
  testthat::expect_equal(ald, ald_check)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::average_leaf_depth(focal_tree),
                         treestats::average_leaf_depth(ltab))
})
