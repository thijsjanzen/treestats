context("tot_coph")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::tot_coph(focal_tree)
    a2 <- treebalance::totCophI(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::tot_coph(focal_tree),
                           treestats::tot_coph(ltab))

    # with extinct species:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.9, fossils = TRUE)

    a1 <- treestats::tot_coph(focal_tree)
    a2 <- treebalance::totCophI(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::tot_coph(focal_tree),
                           treestats::tot_coph(ltab))

    testthat::expect_equal(treestats::tot_coph(focal_tree,
                                               normalization = "yule"),
                           treestats::tot_coph(ltab,
                                               normalization = "yule"))
  }
})


test_that("normalization", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::tot_coph(focal_tree)
  c2 <- treestats::tot_coph(focal_tree, normalization = "yule")
  testthat::expect_lt(c2, c1)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::tot_coph(focal_tree))
    stats2 <- c(stats2, treestats::tot_coph(focal_tree, normalization = "yule"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.2)
})

test_that("abuse", {
  tree1 <- ape::rphylo(n = 2, birth = 1, death = 0)
  testthat::expect_warning(
    treestats::tot_coph(tree1,  normalization = "yule"),
    "normalization not valid for trees of size 2"
  )

  testthat::expect_warning(
    treestats::tot_coph(treestats::phylo_to_l(tree1),  normalization = "yule"),
    "normalization not valid for trees of size 2"
  )
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::tot_coph(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::tot_coph(list()),
    "input object has to be phylo or ltable"
  )
})
