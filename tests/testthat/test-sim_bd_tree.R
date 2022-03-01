context("sim_bd_tree")

test_that("usage num lin", {
  set.seed(42)

  trees <- list()
  cnt <- 1
  methods <-  c("TreeSim", "ape", "TESS", "geiger", "phytools")
  testthat::expect_warning(
  for (m in methods) {
    trees1 <- treestats::sim_bd_tree(birth = 1, death = 0.0, num_trees = 10,
                                   max_lin = 100, max_t = NULL,
                                   method = m)
    trees[[cnt]] <- trees1
    testthat::expect_equal(10, length(trees1))
    for (i in 1:length(trees1)) {
      testthat::expect_equal(100, treestats::number_of_lineages(trees1[[i]]))
    }
    cnt <- cnt + 1
  }
  )

  results <- c()
  for (r in 1:length(trees)) {
    focal_trees <- trees[[r]]
    g_vals <- as.vector(unlist(lapply(focal_trees, treestats::gamma_statistic)))
    nltt_vals <- as.vector(unlist(lapply(focal_trees, treestats::nLTT_base)))
    to_add <- cbind(g_vals, nltt_vals, methods[r])
    results <- rbind(results, to_add)

  }

  ax <- colnames(results)
  for (i in 1:2) {

    to_plot <- as.numeric(results[, i])
    to_plot <- tibble::as_tibble(to_plot)
    to_plot$method <- results[, 3]

    A <- aov(to_plot$value~to_plot$method)
    p <- summary(A)[[1]][["Pr(>F)"]][1]
    testthat::expect_gt(p, 0.05)
  }
})

test_that("usage max t", {
  set.seed(42)

  trees <- list()
  cnt <- 1
  methods <-  c("TreeSim", "DDD", "TESS", "geiger", "phytools")
  for (m in methods) {
    trees1 <- treestats::sim_bd_tree(birth = 1, death = 0.0, num_trees = 10,
                                     max_t = 1,
                                     method = m)
    trees[[cnt]] <- trees1
    testthat::expect_equal(10, length(trees1))
    for (i in 1:length(trees1)) {
      testthat::expect_equal(1, treestats::crown_age(trees1[[i]]))
    }
    cnt <- cnt + 1
  }

  results <- c()
  for (r in 1:length(trees)) {
    focal_trees <- trees[[r]]
    g_vals <- as.vector(unlist(lapply(focal_trees, treestats::gamma_statistic)))
    nltt_vals <- as.vector(unlist(lapply(focal_trees, treestats::nLTT_base)))
    to_add <- cbind(g_vals, nltt_vals, methods[r])
    results <- rbind(results, to_add)

  }

  ax <- colnames(results)
  for (i in 1:2) {

    to_plot <- as.numeric(results[, i])
    to_plot <- tibble::as_tibble(to_plot)
    to_plot$method <- results[, 3]

    A <- aov(to_plot$value~to_plot$method)
    p <- summary(A)[[1]][["Pr(>F)"]][1]
    testthat::expect_gt(p, 0.05)
  }
})
