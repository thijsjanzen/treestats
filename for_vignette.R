require(tidyverse)
require(TreeSim)

found <- c()
for (tree_size in ceiling(10^seq(1, 3, length.out = 10))) {
  cat(tree_size, "\n")
  stats <- c()
  num_repl <- min(100, ceiling(1 + 1000 / tree_size))
  speciation <- 0.5
  max_t <- (1 / speciation) * log(tree_size / 2)

  for (r in 1:num_repl) {
    tree <- TreeSim::sim.bd.taxa.age(n = tree_size,
                                     numbsim = 1,
                                     lambda = speciation,
                                     mu = 0,
                                     age = max_t,
                                     mrca = TRUE)[[1]]
    stats <- rbind(stats, treestats::calc_all_stats(tree, normalize = TRUE))
  }
  stats2 <- apply(stats, 2, median)
  found <- rbind(found, stats2)
}
found <- as_tibble(found)


norm_info <- read.table("https://raw.githubusercontent.com/thijsjanzen/treestats-scripts/main/datasets/normalize.txt", header = TRUE)

found2 <- found |>
  gather(key = "statistic", value = "value", -number_of_lineages)

found3 <- left_join(found2, norm_info)

found3 |>
  filter(normalize == TRUE & type == "Yule") %>%
  ggplot(aes(x = number_of_lineages, y = value)) +
  geom_point(col = "blue", alpha = 0.3) +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 3) +
  theme_minimal() +
  xlab("Extant tree size")
