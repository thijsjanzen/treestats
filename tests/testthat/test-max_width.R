context("max width")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1 <- treestats::max_width(focal_tree)
  a2 <- treebalance::maxWidth(focal_tree)

  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::max_width(focal_tree),
                         treestats::max_width(ltab))


  # with extinct species:
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1 <- treestats::max_width(focal_tree)
  a2 <- treebalance::maxWidth(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::max_width(focal_tree),
                         treestats::max_width(ltab))
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::max_width(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::max_width(list()),
    "input object has to be phylo or ltable"
  )
})
