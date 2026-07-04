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
  c1 <- treestats::average_leaf_depth(extant_tree)
  c2 <- treestats::average_leaf_depth(extant_tree, normalization = "yule")
  testthat::expect_lt(c2, c1)
  c3 <- treestats::average_leaf_depth(treestats::phylo_to_l(extant_tree),
                            normalization = "yule")
  testthat::expect_equal(c2, c3)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::average_leaf_depth(focal_tree))
    stats2 <- c(stats2, treestats::average_leaf_depth(focal_tree,
                                                      normalization = "yule"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.2)
})



test_that("wrong_object", {
  testthat::expect_error(
    treestats::average_leaf_depth(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::average_leaf_depth(list()),
    "input object has to be phylo or ltable"
  )
})
