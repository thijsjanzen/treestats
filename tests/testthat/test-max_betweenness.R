context("betweenness")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1_1 <- treestats::max_betweenness(focal_tree)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:

  a2_1 <- 12795  # max(treeCentrality::computeBetweenness(focal_tree)) #nolint

  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)

  if (requireNamespace("igraph")) {
    df <- as.data.frame(cbind(focal_tree$edge,
                              weight = focal_tree$edge.length))
    g <- igraph::graph_from_data_frame(df, directed = FALSE)
    ref_betweenness <- igraph::betweenness(g)
    testthat::expect_equal(a1_1, max(ref_betweenness), tolerance = 0.01)
  }



  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::max_betweenness(focal_tree, ),
                         treestats::max_betweenness(ltab))

  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1_1 <- treestats::max_betweenness(focal_tree)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  a2_1 <- 20315 # max(treeCentrality::computeBetweenness(focal_tree)) #nolint
  if (requireNamespace("igraph")) {
    df <- as.data.frame(cbind(focal_tree$edge,
                              weight = focal_tree$edge.length))
    g <- igraph::graph_from_data_frame(df, directed = FALSE)
    ref_betweenness <- igraph::betweenness(g)
    testthat::expect_equal(a1_1, max(ref_betweenness), tolerance = 0.01)
  }

  testthat::expect_equal(a1_1, a2_1, tolerance = 0.01)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(treestats::max_betweenness(focal_tree),
                         treestats::max_betweenness(ltab))
})


test_that("normalization", {
  if (requireNamespace("igraph")) {
    set.seed(42)
    focal_tree <- ape::rphylo(n = 30, birth = 1, death = 0)
    a1_1 <- treestats::max_betweenness(focal_tree, normalization = "tips")

    df <- as.data.frame(cbind(focal_tree$edge,
                             weight = focal_tree$edge.length))
    g <- igraph::graph_from_data_frame(df, directed = FALSE)

    ref_betweenness <- igraph::betweenness(g, normalized = TRUE)
    testthat::expect_equal(a1_1, max(ref_betweenness), tolerance = 0.01)

    ltab <- treestats::phylo_to_l(focal_tree)
    a1_2 <- treestats::max_betweenness(ltab, normalization = "tips")
    testthat::expect_equal(a1_2, max(ref_betweenness), tolerance = 0.01)
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::max_betweenness(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::max_betweenness(list()),
    "input object has to be phylo or ltable"
  )
})
