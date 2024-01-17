context("sym_nodes")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::sym_nodes(focal_tree)
    a1_check <- treebalance::symNodesI(focal_tree)
    testthat::expect_equal(a1, a1_check)
    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::sym_nodes(focal_tree),
                           treestats::sym_nodes(ltab))


    # with extinct species:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    a1 <- treestats::sym_nodes(focal_tree)
    a1_check <- treebalance::symNodesI(focal_tree)
    testthat::expect_equal(a1, a1_check)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::sym_nodes(focal_tree),
                           treestats::sym_nodes(ltab))
  }
})


test_that("normalization", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::sym_nodes(focal_tree)
  c2 <- treestats::sym_nodes(focal_tree, normalization = "tips")
  testthat::expect_lt(c2, c1)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::sym_nodes(focal_tree))
    stats2 <- c(stats2, treestats::sym_nodes(focal_tree,
                                             normalization = "tips"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.2)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::sym_nodes(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::sym_nodes(list()),
    "input object has to be phylo or ltable"
  )
})
