context("var leaf depth")

test_that("usage", {
  if (requireNamespace("treebalance")) {
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
    n <- length(focal_tree$tip.label)
    a2 <- a2 * (n) / (n - 1)  # treebalance uses n, we use n - 1

    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::var_leaf_depth(focal_tree),
                           treestats::var_leaf_depth(ltab))
  }

})


test_that("normalization", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::var_leaf_depth(focal_tree)
  c2 <- treestats::var_leaf_depth(focal_tree, normalization = "yule")
  testthat::expect_lt(c2, c1)
  focal_ltab <- treestats::phylo_to_l(focal_tree)
  c3 <- treestats::var_leaf_depth(focal_ltab, normalization = "yule")
  testthat::expect_equal(c2, c3)


  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::var_leaf_depth(focal_tree))
    stats2 <- c(stats2, treestats::var_leaf_depth(focal_tree,
                                                  normalization = "yule"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.2)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::var_leaf_depth(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::var_leaf_depth(list()),
    "input object has to be phylo or ltable"
  )
})
