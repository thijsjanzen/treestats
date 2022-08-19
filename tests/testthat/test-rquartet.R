context("rquartet")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1 <- treestats::rquartet(focal_tree)
  a2 <- treebalance::rQuartetI(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::rquartet(focal_tree),
                         treestats::rquartet(ltab))

  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1 <- treestats::rquartet(focal_tree)
  a2 <- treebalance::rQuartetI(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::rquartet(focal_tree),
                         treestats::rquartet(ltab))

})

test_that("abuse", {
  ctree <- ape::read.tree(text = "c(A:1,B:1,C:1,D:1);")
  testthat::expect_error(treestats::rquartet(ctree),
      "Tree must be binary, for non binary trees use treebalance::rQuartetI")
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

