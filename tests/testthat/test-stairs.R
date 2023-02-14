context("stairs")

test_that("usage stairs1", {
  if (requireNamespace("phyloTop")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)
    c1 <- treestats::stairs(focal_tree)
    c2 <- phyloTop::stairs(focal_tree)[[1]]
    testthat::expect_equal(c1, c2)

    c3 <- treestats::stairs(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)

    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0.5,
                              fossils = TRUE)

    c1 <- treestats::stairs(focal_tree)
    c2 <- phyloTop::stairs(focal_tree)[[1]]
    testthat::expect_equal(c1, c2)

    c3 <- treestats::stairs(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)
  }
})

test_that("usage stairs2", {
  if (requireNamespace("phyloTop")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

    c1 <- treestats::stairs2(focal_tree)
    c2 <- phyloTop::stairs(focal_tree)[[2]]
    testthat::expect_equal(c1, c2)

    c3 <- treestats::stairs2(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)


    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0.5, fossils = TRUE)

    c1 <- treestats::stairs2(focal_tree)
    c2 <- phyloTop::stairs(focal_tree)[[2]]
    testthat::expect_equal(c1, c2)

    c3 <- treestats::stairs2(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::stairs(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::stairs(list()),
    "input object has to be phylo or ltable"
  )
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::stairs2(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::stairs2(list()),
    "input object has to be phylo or ltable"
  )
})
