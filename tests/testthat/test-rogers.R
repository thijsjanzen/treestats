context("rogers")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    standard_test(treestats::rogers, treebalance::rogersI)

    test_polytomies(treestats::rogers, treebalance::rogersI)
  }
})


test_that("normalization", {
  testthat::skip_on_cran()
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::rogers(focal_tree)
  c2 <- treestats::rogers(focal_tree, normalization = "tips")
  testthat::expect_lt(c2, c1)
  c3 <- treestats::rogers(treestats::phylo_to_l(focal_tree),
                          normalization = "tips")
  testthat::expect_equal(c2, c3)

  c2 <- treestats::rogers(focal_tree, normalization = "yule")
  testthat::expect_lt(c2, c1)
  c3 <- treestats::rogers(treestats::phylo_to_l(focal_tree),
                          normalization = "yule")
  testthat::expect_equal(c2, c3)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::rogers(focal_tree))
    stats2 <- c(stats2, treestats::rogers(focal_tree, normalization = "tips"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.2)



  found <- c()
  found2 <- c()
  for (r in 1:100) {
    focal_tree <- ape::rphylo(n = 31, birth = 1, death = 0)
    found[r] <- treestats::rogers(focal_tree, normalization = "yule")
    found2[r] <- treestats::rogers(treestats::phylo_to_l(focal_tree),
                                       normalization = "yule")
  }
  testthat::expect_equal(mean(found), 1.0, tolerance = 0.1)
  testthat::expect_equal(mean(found), mean(found2), tolerance = 0.01)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::rogers(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::rogers(list()),
    "input object has to be phylo or ltable"
  )
})
