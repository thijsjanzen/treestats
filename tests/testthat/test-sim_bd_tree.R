context("sim_bd_tree")

test_that("usage", {
  set.seed(42)

  trees <- list()
  cnt <- 1
  methods <-  c("TreeSim", "ape", "TESS", "geiger", "phytools")
  for (m in methods) {
    trees1 <- treestats::sim_bd_tree(birth = 1, death = 0.0, num_trees = 10,
                                   max_lin = 100, max_t = NULL,
                                   method = m)
    trees[[cnt]] <- trees1
    testthat::expect_equal(10, length(trees1))
    cnt <- cnt + 1
  }

  results <- c()
  for (r in 1:length(trees)) {
    focal_trees <- trees[[r]]
    all_stats <- lapply(focal_trees, treestats::calc_all_stats)
    for (s in 1:length(all_stats)) {
      focal_stats <- all_stats[[s]]
      results <- rbind(results, c(unlist(focal_stats), methods[r]))
    }
  }

  ax <- colnames(results)
  for (i in 1:19) {

    to_plot <- as.numeric(results[, i])
    to_plot <- tibble::as_tibble(to_plot)
    to_plot$method <- results[, 20]

    A <- aov(to_plot$value~to_plot$method)
    p <- summary(A)[[1]][["Pr(>F)"]][1]
    testthat::expect_gt(p, 0.05)
  }
})
