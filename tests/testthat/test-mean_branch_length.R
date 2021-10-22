context("mean_branch_length")

test_that("usage", {
  focal_tree <- phytools::read.newick(text = "(1:4,2:4):0;")

  mbr <- treestats::mean_branch_length(focal_tree)
  testthat::expect_equal(mbr, 4)
})
