context("avg vertex depth")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    a1 <- treestats::avg_vert_depth(extant_tree)
    a2 <- treebalance::avgVertDep(extant_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(extant_tree)
    testthat::expect_equal(treestats::avg_vert_depth(extant_tree),
                           treestats::avg_vert_depth(ltab))

    a1 <- treestats::avg_vert_depth(complete_tree)
    a2 <- treebalance::avgVertDep(complete_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(complete_tree)
    testthat::expect_equal(treestats::avg_vert_depth(complete_tree),
                           treestats::avg_vert_depth(ltab))
  }
})

test_that("polytomies", {
  if (requireNamespace("treebalance")) {
    test_polytomies(treestats::avg_vert_depth, treebalance::avgVertDep)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::avg_vert_depth(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::avg_vert_depth(list()),
    "input object has to be phylo or ltable"
  )
})
