context("rquartet")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::rquartet(focal_tree)
    a2 <- treebalance::rQuartetI(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::rquartet(focal_tree),
                           treestats::rquartet(ltab))

    # normalization
    a3 <- treestats::rquartet(focal_tree, normalization = "yule")
    a4 <- treestats::rquartet(ltab, normalization = "yule")
    testthat::expect_equal(a3, a4)
    testthat::expect_lt(a3, a1)
    testthat::expect_lt(a3, 2)

    a3 <- treestats::rquartet(focal_tree, normalization = "pda")
    a4 <- treestats::rquartet(ltab, normalization = "pda")
    testthat::expect_equal(a3, a4)
    testthat::expect_lt(a3, a1)
    testthat::expect_lt(a3, 2)

    # with extinct species:
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

    a1 <- treestats::rquartet(focal_tree)
    a2 <- treebalance::rQuartetI(focal_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::rquartet(focal_tree),
                           treestats::rquartet(ltab))
  }
})

test_that("polytomy", {
  testthat::skip_on_ci()
  testthat::skip_on_cran()
  focal_tree <- ape::read.tree(text = "(t1:1.5,(t2:1,t3:1,t4:1):0.5);")
  testthat::expect_error(
    treestats::rquartet(focal_tree),
    "Tree must be binary, for non binary trees use treebalance::rQuartetI"
  )
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
