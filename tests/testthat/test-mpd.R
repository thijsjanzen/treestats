context("mean pair distance")

test_that("usage", {
  if (requireNamespace("picante")) {
    a1 <- treestats::mean_pair_dist(extant_tree)

    n <- length(extant_tree$tip.label)
    sample_mat <- matrix(data = 1, nrow = n, ncol = n)
    colnames(sample_mat) <- extant_tree$tip.label

    a2 <- picante::mpd(sample_mat, cophenetic(extant_tree),
                       abundance.weighted = FALSE)[[1]]
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(extant_tree)
    testthat::expect_equal(treestats::mean_pair_dist(extant_tree),
                           treestats::mean_pair_dist(ltab))
  }
})

test_that("polytomies", {
  for (focal_tree in poly_trees) {
    res <- treestats::mean_pair_dist(focal_tree)
    testthat::expect_true(!is.na(res))
  }
})

test_that("unrooted", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
  a1 <- treestats::mean_pair_dist(focal_tree)
  a2 <- treestats::mean_pair_dist(ape::unroot(focal_tree))
  testthat::expect_equal(a1, a2)
})

test_that("normalization", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::mean_pair_dist(focal_tree)
  c2 <- treestats::mean_pair_dist(focal_tree, normalization = "tips")
  testthat::expect_lt(c2, c1)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::mean_pair_dist(focal_tree))
    stats2 <- c(stats2, treestats::mean_pair_dist(focal_tree,
                                                  normalization = "tips"))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.3)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::mean_pair_dist(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mean_pair_dist(list()),
    "input object has to be phylo or ltable"
  )
})
