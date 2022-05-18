context("mean_branch_length")

test_that("usage", {
  calc_mean_br_r <- function(focal_tree) {
    return(mean(focal_tree$edge.length, na.rm = TRUE))
  }

  focal_tree <- ape::read.tree(text = "(1:4,2:4):0;")

  a1 <- treestats::mean_branch_length(focal_tree)
  a2 <- calc_mean_br_r(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::mean_branch_length(focal_tree),
                         treestats::mean_branch_length(ltab))



  focal_tree <- TreeSim::sim.bd.taxa(n = 140,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  a1 <- treestats::mean_branch_length(focal_tree)
  a2 <- calc_mean_br_r(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(a2,
                         treestats::mean_branch_length(ltab))


  # with extinct species:
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0.2)[[1]]

  a1 <- treestats::mean_branch_length(focal_tree)
  a2 <- calc_mean_br_r(focal_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(a2,
                         treestats::mean_branch_length(ltab))
})
