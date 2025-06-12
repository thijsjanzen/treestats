tree_size <- 4

focal_tree1 <- ape::rphylo(n = tree_size, birth = 1, death = 0)
focal_tree2 <- ape::unroot(focal_tree1)

cbind(focal_tree1$edge, focal_tree1$edge.length)
cbind(focal_tree2$edge, focal_tree2$edge.length)

treestats::diameter(focal_tree1, weight = FALSE)
treestats::diameter(focal_tree2, weight = FALSE)


max(ape::cophenetic.phylo(focal_tree1))
max(ape::cophenetic.phylo(focal_tree2))




tree_size <- 100
set.seed(1)
focal_tree <- ape::rphylo(n = tree_size, birth = 1, death = 0)
num_nodes <- focal_tree$Nnode
found <- c()
for (r in 1:tree_size) {
  #cat(r, "\n")
  focal_tree2 <- ape::root.phylo(focal_tree, outgroup = r, resolve.root = TRUE)
  # focal_tree2 <- stats::reorder(focal_tree2)
  #treestats::pigot_rho(focal_tree2)

  #res <- treestats::avg_vert_depth(focal_tree2)

  res <- treestats::calc_all_stats(focal_tree2)
  found <- rbind(found, res)
}

for (r in 2:num_nodes) {
  #cat(r, "\n")
  focal_tree2 <- ape::root.phylo(focal_tree, node = tree_size + r, resolve.root = TRUE)
  # focal_tree2 <- stats::reorder(focal_tree2)
  #treestats::pigot_rho(focal_tree2)

  #res <- treestats::avg_vert_depth(focal_tree2)

  res <- treestats::calc_all_stats(focal_tree2)
  found <- rbind(found, res)
}




#hist(found)
require(tidyverse)
found <- as_tibble(found)
found %>%
  gather(key = "statistic", value = "value") %>%
  ggplot(aes(x = value)) +
    geom_histogram() +
    facet_wrap(~statistic, scales = "free")

