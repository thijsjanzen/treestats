context("eigenvector")

test_that("usage", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)

  a1_1 <- treestats::eigen_centrality(focal_tree, weight = TRUE, scale = FALSE)
  a1_2 <- treestats::eigen_centrality(focal_tree, weight = FALSE, scale = FALSE)



  # because treeCentrality is not available on CRAN, we precompute reference
  # values:
  a2_1 <- 3.585569  # treeCentrality::computeEigenvector(focal_tree
  #                                      weight = TRUE))
  a2_2 <- 2.609486  # treeCentrality::computeEigenvector(focal_tree
  #                                      weight = FALSE))

  testthat::expect_equal(a1_1$eigenvalue, a2_1, tolerance = 1e-4)
  testthat::expect_equal(a1_2$eigenvalue, a2_2, tolerance = 1e-4)
  a1_3 <- treestats::eigen_centrality(focal_tree, weight = TRUE, scale = TRUE)
  testthat::expect_equal(a1_1$eigenvalue,
                         a1_3$eigenvalue)
  testthat::expect_equal(max(a1_3$eigenvector), 1)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(
              treestats::eigen_centrality(focal_tree, weight = TRUE)$eigenvalue,
              treestats::eigen_centrality(ltab, weight = TRUE)$eigenvalue)

  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

  a1_1 <- treestats::eigen_centrality(focal_tree, weight = TRUE)
  a1_2 <- treestats::eigen_centrality(focal_tree, weight = FALSE)
  # because treeCentrality is not available on CRAN, we precompute reference
  # values:

  a2_1 <- 3.718206  # treeCentrality::computeEigenvector(focal_tree
  #                                      weight = TRUE))
  a2_2 <- 2.631605  # treeCentrality::computeEigenvector(focal_tree
  #                                      weight = FALSE))


  testthat::expect_equal(a1_1$eigenvalue, a2_1, tolerance = 1e-4)
  testthat::expect_equal(a1_2$eigenvalue, a2_2, tolerance = 1e-4)

  ltab <- treestats::phylo_to_l(focal_tree)
  testthat::expect_equal(
     treestats::eigen_centrality(focal_tree, weight = TRUE)$eigenvalue,
     treestats::eigen_centrality(ltab,       weight = TRUE)$eigenvalue)

  if (requireNamespace("igraph")) {
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
    df <- as.data.frame(cbind(focal_tree$edge,
                             weight = focal_tree$edge.length))
    g <- igraph::graph_from_data_frame(df, directed = FALSE)

    ref <- igraph::eigen_centrality(g)

    a1 <- treestats::eigen_centrality(focal_tree)

    testthat::expect_equal(a1$eigenvalue, ref$value)

    ltab <- treestats::phylo_to_l(focal_tree)
    a1_2 <- treestats::eigen_centrality(ltab)
    testthat::expect_equal(a1_2$eigenvalue, ref$value, tolerance = 0.01)
  }

  # compare namespaces
  if (requireNamespace("Matrix")) {
    focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
    a1_1 <- treestats::eigen_centrality(focal_tree,
                                        weight = TRUE,
                                        scale = FALSE)
    a2_1 <- treestats::eigen_centrality(focal_tree,
                                        weight = FALSE,
                                        scale = FALSE)

    testthat::with_mocked_bindings(
      {
        # Now `myfun()` should behave as if `data.tree` is not installed
        a1_2 <- treestats::eigen_centrality(focal_tree,
                                            weight = TRUE,
                                            scale = FALSE)
        testthat::expect_equal(a1_1, a1_2)

        a1_3 <- treestats::eigen_centrality(focal_tree,
                                            weight = FALSE,
                                            scale = FALSE)
        testthat::expect_equal(a2_1, a1_3)
      },
      requireNamespace = function(pkg, quietly = TRUE) {
        if (pkg == "Matrix") {
          return(FALSE)
        }
        if (pkg == "RSpectra") {
          return(FALSE)
        }
        # Call the real `requireNamespace` for other packages
        base::requireNamespace(pkg, quietly = TRUE)
      },
      .package = "base"
    )
  }
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::eigen_centrality(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::eigen_centrality(list()),
    "input object has to be phylo or ltable"
  )
})
