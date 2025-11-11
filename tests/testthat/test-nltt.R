context("nltt")

test_that("usage", {
  if (requireNamespace("nLTT")) {
    set.seed(42)
    tree1 <- ape::rphylo(n = 100, birth = 1, death = 0)
    tree2 <- ape::rphylo(n = 100, birth = 0.5, death = 0)

    nltt <- treestats::nLTT(tree1, tree2)
    nltt_check <- nLTT::nLTTstat(tree1, tree2)
    testthat::expect_equal(nltt, nltt_check)


    empty_tree <- ape::rphylo(n = 2, birth = 1, death = 0)
    nltt_base <- treestats::nLTT_base(tree1)
    nltt_base2 <- treestats::nLTT(tree1, empty_tree)
    nltt_base3 <-  nLTT::nLTTstat(tree1, empty_tree)

    testthat::expect_equal(nltt_base, nltt_base2, tolerance = 1e-3)
    testthat::expect_equal(nltt_base, nltt_base3, tolerance = 1e-3)

    empty_tree <- ape::read.tree(text = "(1:4,2:4):0;")
    nltt_stat <- nLTT::nltt_diff(tree1, empty_tree)
    testthat::expect_equal(nltt_stat, nltt_base, tolerance = 1e-3)
    testthat::expect_equal(nltt_stat, nltt_base2, tolerance = 1e-3)

    # test ltable functionality
    ltab1 <- treestats::phylo_to_l(tree1)
    ltab2 <- treestats::phylo_to_l(tree2)
    nltt3 <- treestats::nLTT(ltab1, ltab2)
    testthat::expect_equal(nltt3, nltt)

    # mixed ltable:
    nltt4 <- treestats::nLTT(ltab1, tree2)
    testthat::expect_equal(nltt4, nltt)

    nltt5 <- treestats::nLTT(tree1, ltab2)
    testthat::expect_equal(nltt5, nltt)
  }
})

test_that("wrong_object", {
  set.seed(42)
  tree1 <- ape::rphylo(n = 10, birth = 1, death = 0)
  testthat::expect_error(
    treestats::nLTT(10, tree1),
    "input needs to be phylo or ltable object"
  )

  testthat::expect_error(
    treestats::nLTT(list(), tree1),
    "input needs to be phylo or ltable object"
  )
})
