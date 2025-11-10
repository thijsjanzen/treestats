context("blum")

test_that("usage", {
  if (requireNamespace("treebalance")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

    blum1 <- treestats::blum(focal_tree)
    blum_check <- treebalance::sShapeI(focal_tree, logbase = exp(1))
    testthat::expect_equal(blum1, blum_check)
    ltab <- treestats::phylo_to_l(focal_tree)
    blum2 <- treestats::blum(ltab)
    testthat::expect_equal(blum1, blum2)
  }
})

test_that("usage", {
  set.seed(5)
  focal_tree <- ape::rphylo(n = 4, birth = 1, death = 0)

  blum1 <- treestats::blum(focal_tree)
  blum_check <- log(2 - 1) + log(3 - 1) + log(4 - 1)
  testthat::expect_equal(blum1, blum_check)
  ltab <- treestats::phylo_to_l(focal_tree)
  blum2 <- treestats::blum(ltab)
  testthat::expect_equal(blum1, blum2)

  sshape1 <- treestats::sshape(focal_tree)
  testthat::expect_equal(sshape1, blum1)
  sshape2 <- treestats::sshape(ltab)
  testthat::expect_equal(sshape2, blum1)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::blum(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::blum(list()),
    "input object has to be phylo or ltable"
  )
})

test_that("normalisation", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)

  c1 <- treestats::blum(focal_tree, normalization = FALSE)
  c2 <- treestats::blum(focal_tree, normalization = TRUE)
  testthat::expect_lt(c2, c1)
  c3 <- treestats::blum(treestats::phylo_to_l(focal_tree),
                        normalization = TRUE)
  testthat::expect_equal(c2, c3)

  stats1 <- c()
  stats2 <- c()
  for (n in seq(100, 200, by = 10)) {
    focal_tree <- ape::rphylo(n = n, birth = 1, death = 0)
    stats1 <- c(stats1, treestats::blum(focal_tree))
    stats2 <- c(stats2, treestats::blum(focal_tree, normalization = TRUE))
  }

  a1 <- cor(stats1, seq(100, 200, by = 10))
  a2 <- cor(stats2, seq(100, 200, by = 10))

  testthat::expect_lt(a2, a1)
  testthat::expect_lt(a2, 0.5)
})
