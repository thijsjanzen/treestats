context("area per pair")

test_that("usage", {

  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::area_per_pair(focal_tree)
    a2 <- treebalance::areaPerPairI(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::area_per_pair(focal_tree),
                           treestats::area_per_pair(ltab))

    # with extinct species:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    a1 <- treestats::area_per_pair(focal_tree)
    a2 <- treebalance::areaPerPairI(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::area_per_pair(focal_tree),
                           treestats::area_per_pair(ltab))
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::area_per_pair(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::area_per_pair(list()),
    "input object has to be phylo or ltable"
  )
})
