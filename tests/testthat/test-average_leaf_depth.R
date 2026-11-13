context("average leaf depth")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    ald <- treestats::average_leaf_depth(extant_tree)
    ald_check <- treebalance::avgLeafDepI(extant_tree)
    testthat::expect_equal(ald, ald_check)

    ltab <- treestats::phylo_to_l(extant_tree)
    testthat::expect_equal(treestats::average_leaf_depth(extant_tree),
                           treestats::average_leaf_depth(ltab))


    # with extinct species:
    ald <- treestats::average_leaf_depth(complete_tree)
    ald_check <- treebalance::avgLeafDepI(complete_tree)
    testthat::expect_equal(ald, ald_check)

    ltab <- treestats::phylo_to_l(complete_tree)
    testthat::expect_equal(treestats::average_leaf_depth(complete_tree),
                           treestats::average_leaf_depth(ltab))
  }
})

test_that("polytomies", {
  if (requireNamespace("treebalance")) {
    test_polytomies(treestats::average_leaf_depth, treebalance::avgLeafDepI)
  }
})

test_that("normalization", {
  test_yule(treestats::average_leaf_depth)
})



test_that("abuse", {
  test_abuse(treestats::average_leaf_depth)
})
