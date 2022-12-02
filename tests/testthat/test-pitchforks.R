context("pitchforks")

test_that("usage", {
  if (requireNamespace("phyloTop")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

    c1 <- treestats::pitchforks(focal_tree)
    c2 <- phyloTop::pitchforks(focal_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::pitchforks(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)


    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0.5, fossils = TRUE)

    c1 <- treestats::pitchforks(focal_tree)
    c2 <- phyloTop::pitchforks(focal_tree)
    testthat::expect_equal(c1, c2)

    c3 <- treestats::pitchforks(treestats::phylo_to_l(focal_tree))
    testthat::expect_equal(c1, c3)
  }
})


test_that("normalization", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::pitchforks(focal_tree)
  c2 <- treestats::pitchforks(focal_tree, normalization = "tips")
  testthat::expect_lt(c2, c1)
  c3 <- treestats::pitchforks(treestats::phylo_to_l(focal_tree),
                            normalization = "tips")
  testthat::expect_equal(c2, c3)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::pitchforks(focal_tree))
    stats2 <- c(stats2, treestats::pitchforks(focal_tree,
                                              normalization = "tips"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.2)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::pitchforks(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::pitchforks(list()),
    "input object has to be phylo or ltable"
  )
})
