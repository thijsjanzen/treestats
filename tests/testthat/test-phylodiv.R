context("phylodiv")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 3, birth = 1, death = 0)

  div1 <- treestats::phylogenetic_diversity(focal_tree)
  div2 <- sum(focal_tree$edge.length)
  testthat::expect_equal(div1, div2, tolerance = 1e-4)

  ca <- max(treestats::branching_times(focal_tree))
  pds <- treestats::phylogenetic_diversity(focal_tree,
                                           t = seq(ca, 0, length.out =  100))
  testthat::expect_equal(length(pds), 100)
  for (i in 2:length(pds)) {
    testthat::expect_gt(pds[i], pds[i - 1])
  }

  # now check sub time
  brts <- ape::branching.times(focal_tree)
  tt <- (brts[1] + brts[2]) / 2
  div1 <- treestats::phylogenetic_diversity(focal_tree, tt)
  div2 <- (brts[1] - tt) * 2
  testthat::expect_equal(div1[[1]], div2[[1]])

  # now with extinct lineages
  focal_tree <- ape::read.tree(text =
                          "((t1:2.0, t2:2.0):1.0, (t3:1.0, t4:2.0):1.0):1.0;")

  div1 <- treestats::phylogenetic_diversity(focal_tree)
  testthat::expect_equal(div1, 8)

  div1 <- treestats::phylogenetic_diversity(focal_tree, 2.0)
  div2 <- treestats::phylogenetic_diversity(focal_tree, 2.02)
  testthat::expect_gt(div1, div2)

  # now with double extinct lineages
  focal_tree <- ape::read.tree(text =
                      "((:2.0, :2.0):1.0, ((:0.5, :1.0):0.5, :2.0):1.0):1.0;")
  div1 <- treestats::phylogenetic_diversity(focal_tree)
  testthat::expect_equal(div1, 8)
  div1 <- treestats::phylogenetic_diversity(focal_tree, 2)
  div2 <- treestats::phylogenetic_diversity(focal_tree, 2.01)
  testthat::expect_gt(div1, div2)
})

test_that("usage 2", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 50, birth = 1, death = 0)

  div1 <- treestats::phylogenetic_diversity(focal_tree)
  div2 <- sum(focal_tree$edge.length)
  testthat::expect_equal(div1, div2, tolerance = 1e-4)

  ltab <- treestats::phylo_to_l(focal_tree)
  div3 <- treestats::phylogenetic_diversity(ltab)
  testthat::expect_equal(div3, div2)
})

test_that("ltab", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 3, birth = 1, death = 0)

  focal_ltab <- treestats::phylo_to_l(focal_tree)

  div1 <- treestats::phylogenetic_diversity(focal_ltab)

  testthat::expect_error(treestats::phylogenetic_diversity(focal_ltab,
                                                           t = 1))

  testthat::expect_error(treestats::phylogenetic_diversity(focal_ltab,
                                                           t = c(0.0, 1.0)))

  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0.3, fossils = TRUE)
  focal_ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_error(treestats::phylogenetic_diversity(focal_ltab))
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::phylogenetic_diversity(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::phylogenetic_diversity(list()),
    "input object has to be phylo or ltable"
  )
})
