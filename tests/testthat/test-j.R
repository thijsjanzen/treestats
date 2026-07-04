context("entropy j")

test_that("usage", {
  if (requireNamespace("picante")) {
    a1 <- treestats::entropy_j(extant_tree)

    n <- length(extant_tree$tip.label)
    sample_mat <- matrix(data = 1, nrow = n, ncol = n)
    colnames(sample_mat) <- extant_tree$tip.label

    a2 <- picante::mpd(sample_mat, stats::cophenetic(extant_tree),
                       abundance.weighted = FALSE)[[1]]
    testthat::expect_equal(a1, a2 / n)

    ltab <- treestats::phylo_to_l(extant_tree)
    testthat::expect_equal(treestats::entropy_j(extant_tree),
                           treestats::entropy_j(ltab))
  }
})

test_that("polytomies", {
  if (requireNamespace("picante")) {

    test_func <- function(phy) {
      n <- length(phy$tip.label)
      sample_mat <- matrix(data = 1, nrow = n, ncol = n)
      colnames(sample_mat) <- phy$tip.label

      a2 <- picante::mpd(sample_mat, stats::cophenetic(phy),
                         abundance.weighted = FALSE)[[1]]
      return(a2 / n)
    }

    test_polytomies(treestats::entropy_j, test_func)
  }
})

test_that("unrooted", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
  a1 <- treestats::entropy_j(focal_tree)
  a2 <- treestats::entropy_j(ape::unroot(focal_tree))
  testthat::expect_equal(a1, a2)
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
