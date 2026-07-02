context("psv")

test_that("usage", {
  if (requireNamespace("picante")) {

    test_func <- function(phy) {
      n <- length(phy$tip.label)
      sample_mat <- matrix(data = 1, nrow = n, ncol = n)
      colnames(sample_mat) <- phy$tip.label

      return(picante::psv(sample_mat,
                          phy,
                          compute.var = FALSE,
                          scale.vcv = FALSE)[1, 1])
    }

    standard_test(treestats::psv, test_func, drop_complete = TRUE)

    test_polytomies(treestats::psv, test_func)
  }
})

test_that("unrooted", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
  a1 <- treestats::psv(focal_tree)
  a2 <- treestats::psv(ape::unroot(focal_tree))
  testthat::expect_equal(a1, a2)
})

test_that("normalisation", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::psv(focal_tree)
  c2 <- treestats::psv(focal_tree, normalization = "tips")
  testthat::expect_lt(c2, c1)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::psv(focal_tree))
    stats2 <- c(stats2, treestats::psv(focal_tree, normalization = "tips"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.4)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::psv(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::psv(list()),
    "input object has to be phylo or ltable"
  )
})
