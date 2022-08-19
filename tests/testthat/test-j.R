context("mean pair distance")

test_that("usage", {
  if (requireNamespace("picante")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::entropy_j(focal_tree)

    n <- length(focal_tree$tip.label)
    sample_mat <- matrix(data = 1, nrow = n, ncol = n)
    colnames(sample_mat) <- focal_tree$tip.label

    a2 <- picante::mpd(sample_mat, cophenetic(focal_tree),
                       abundance.weighted = FALSE)[[1]]
    testthat::expect_equal(a1, a2 / n)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::entropy_j(focal_tree),
                           treestats::entropy_j(ltab))
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::entropy_j(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::entropy_j(list()),
    "input object has to be phylo or ltable"
  )
})
