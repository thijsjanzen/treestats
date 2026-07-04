context("area per pair")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    a1 <- treestats::area_per_pair(extant_tree)
    a2 <- treebalance::areaPerPairI(extant_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(extant_tree)
    testthat::expect_equal(treestats::area_per_pair(extant_tree),
                           treestats::area_per_pair(ltab))

    # with extinct species:
    a1 <- treestats::area_per_pair(complete_tree)
    a2 <- treebalance::areaPerPairI(complete_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(complete_tree)
    testthat::expect_equal(treestats::area_per_pair(complete_tree),
                           treestats::area_per_pair(ltab))

    # on polytomies
    test_polytomies(treestats::area_per_pair, treebalance::areaPerPairI)
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
