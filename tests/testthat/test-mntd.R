context("mean nearest taxon distance")

test_that("usage", {
  if (requireNamespace("picante")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    a1 <- treestats::mntd(focal_tree)

    n <- length(focal_tree$tip.label)
    sample_mat <- matrix(data = 1, nrow = n, ncol = n)
    colnames(sample_mat) <- focal_tree$tip.label

    a2 <- picante::mntd(sample_mat, cophenetic(focal_tree),
                        abundance.weighted = FALSE)[[1]]
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(focal_tree)
    testthat::expect_equal(treestats::mntd(focal_tree),
                           treestats::mntd(ltab))
  }
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  testthat::expect_error(
    treestats::mntd(focal_tree),
    "can only calculate mntd statistic for ultrametric tree")

  focal_ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_error(
    treestats::mntd(focal_ltab),
    "can only calculate mntd statistic for ultrametric tree")

})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::mntd(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mntd(list()),
    "input object has to be phylo or ltable"
  )
})
