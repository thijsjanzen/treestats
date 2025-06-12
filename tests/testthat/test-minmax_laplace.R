context("minmax_laplace")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 5, birth = 1, death = 0)

  a1_1 <- treestats::minmax_laplace(focal_tree)

  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  # ref <-
  #    treeCentrality::computeSpectrum(focal_tree, weight = TRUE, dist = FALSE,
  #                                    full = TRUE, lap = TRUE, norm = FALSE)


  ref <- c(5.684185e-01, 4.343028e-01, 1.754782e-01, 1.248247e-01, 7.433148e-02,
           3.966736e-02, 3.384060e-02, 4.262393e-03, 3.699463e-17)
  ref <- round(ref, digits = 10)
  testthat::expect_equal(a1_1$min, min(ref[ref > 0]), tolerance = 0.01)

  testthat::expect_equal(a1_1$max, max(ref), tolerance = 0.01)

  ltab <- treestats::phylo_to_l(focal_tree)

  testthat::expect_equal(
    treestats::minmax_laplace(focal_tree),
    treestats::minmax_laplace(ltab))

  focal_tree <- ape::rphylo(n = 5, birth = 1, death = 0.2, fossils = TRUE)

  ref <- c(3.609466e+00,  2.997761e+00,  2.771082e+00,  1.951827e+00,
           1.704470e+00,  1.384831e+00,  9.340857e-01,  7.710821e-01,
           7.093589e-01,  3.099265e-01,  2.548943e-01,  1.853926e-01,
           1.459424e-01,  3.611505e-02,  1.862803e-02,  6.809262e-03,
           -1.202488e-16)
  ref <- round(ref, digits = 10)
  # ref <-
  #     treeCentrality::computeSpectrum(focal_tree, weight = TRUE, dist = FALSE,
  #                                      full = TRUE, lap = TRUE, norm = FALSE)
  a1_1 <- treestats::minmax_laplace(focal_tree)

  testthat::expect_equal(a1_1$min, min(ref[ref > 0]), tolerance = 0.001)

  testthat::expect_equal(a1_1$max, max(ref), tolerance = 0.001)

  ltab <- treestats::phylo_to_l(focal_tree)

  testthat::expect_warning(
    a2_1 <- treestats::minmax_laplace(ltab)
  )

  testthat::expect_equal(a2_1$min, min(ref[ref > 0]))
  testthat::expect_equal(a2_1$max, max(ref), tolerance = 0.001)

  # test rspectra use
  ref <- treestats::minmax_laplace(focal_tree,
                                   use_rspectra = FALSE)

  a2 <- treestats::minmax_laplace(focal_tree,
                                  use_rspectra = TRUE)
  testthat::expect_equal(ref$min, a2$min)
  testthat::expect_equal(ref$max, a2$max)


})

test_that("igraph", {
  if (requireNamespace("igraph")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 10, birth = 1, death = 0)

    a1_1 <- treestats::minmax_laplace(focal_tree)

    df <- as.data.frame(cbind(focal_tree$edge,
                              weight = focal_tree$edge.length))
    g <- igraph::graph_from_data_frame(df, directed = FALSE)

    lapl_mat <- igraph::laplacian_matrix(g,
                                         normalization = "unnormalized",
                                         sparse = FALSE)
    ref <- eigen(lapl_mat)$values
    ref <- round(ref, digits = 10)

    testthat::expect_equal(a1_1$min, min(ref[ref > 1e-5]))
    testthat::expect_equal(a1_1$max, max(ref))
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::minmax_laplace(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::minmax_laplace(list()),
    "input object has to be phylo or ltable"
  )
})
