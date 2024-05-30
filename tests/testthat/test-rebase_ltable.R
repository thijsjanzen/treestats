context("rebase_ltable")

test_that("usage", {
  set.seed(42)
  # this is a random tree that happens to trigger a rebase:
  focal_tree <- ape::rphylo(n = 20, birth = 1, death = 0)

  ltab <- treestats::phylo_to_l(focal_tree)

  ltab2 <- treestats::rebase_ltable(ltab)
  testthat::expect_true(all.equal(ltab, ltab2) != TRUE)

  # test small tree:
  focal_tree <- ape::rphylo(n = 2, birth = 1, death = 0)

  ltab <- treestats::phylo_to_l(focal_tree)

  focal_tree2 <- treestats::rebase_ltable(ltab)
  testthat::expect_true(all.equal(ltab, focal_tree2))
})
