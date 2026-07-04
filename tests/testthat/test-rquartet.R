context("rquartet")

test_that("usage", {
  if (requireNamespace("treebalance")) {

    standard_test(treestats::rquartet, treebalance::rQuartetI)

    # normalization
    ltab <- treestats::phylo_to_l(extant_tree)
    a3 <- treestats::rquartet(extant_tree, normalization = "yule")
    a4 <- treestats::rquartet(ltab, normalization = "yule")
    testthat::expect_equal(a3, a4)
    testthat::expect_lt(a3, 2)

    a3 <- treestats::rquartet(extant_tree, normalization = "pda")
    a4 <- treestats::rquartet(ltab, normalization = "pda")
    testthat::expect_equal(a3, a4)
    testthat::expect_lt(a3, 2)
  }
})

test_that("polytomy", {
  if (requireNamespace("treebalance")) {
    focal_tree <- ape::read.tree(text = "(t1:1.5,(t2:1,t3:1,t4:1):0.5);")
    res <- treestats::rquartet(focal_tree)
    ref <- treebalance::rQuartetI(focal_tree)
    testthat::expect_equal(res, ref)

    test_polytomies(treestats::rquartet, treebalance::rQuartetI)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::rquartet(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::rquartet(list()),
    "input object has to be phylo or ltable"
  )
})
