context("tot_coph")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    standard_test(treestats::tot_coph, treebalance::totCophI)

    test_polytomies(treestats::tot_coph, treebalance::totCophI)
  }
})


test_that("normalization", {
  testthat::skip_on_cran()
  test_yule(treestats::tot_coph)
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
  test_abuse(treestats::tot_coph)
})
