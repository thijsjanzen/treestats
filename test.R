set.seed(1)

a1 <- treestats::eigen_centrality(phy, use_rspectra = FALSE)
a2 <- treestats::eigen_centrality(phy, use_rspectra = TRUE)
all.equal(a1, a2)


library(tidyverse)
res <- microbenchmark::microbenchmark(treestats::eigen_centrality(phy, use_rspectra = FALSE),
                                      treestats::eigen_centrality(phy, use_rspectra = TRUE),
                                      times = 10L)

phy <- ape::rphylo(500, 1, 0)
profvis::profvis(
 for (r in 1:500) {
   treestats::eigen_centrality(phy, weight = TRUE, use_rspectra = TRUE)
 }

)


func1 <- function() {
  edge_for_mat <- rbind(phy$edge, cbind(phy$edge[, 2], phy$edge[, 1]))
  adj_matrix <- Matrix::sparseMatrix(i = edge_for_mat[, 1],
                                     j = edge_for_mat[, 2],
                                     x = c(phy$edge.length,
                                           phy$edge.length))
}

func2 <- function() {
  adj_matrix <- treestats::prep_adj_mat(as.vector(t(phy$edge)),
                             as.vector(phy$edge.length),
                             TRUE)
}

func3 <- function() {
  n_entries <- 2 * phy$Nnode + 1
  mat <- matrix(0, n_entries, n_entries)
  for (i in seq_along(nrow(phy$edge))) {
    x <- phy$edge[i, 1]
    y <- phy$edge[i, 2]
    z <- phy$edge.length[i]
    mat[x, y] <- mat[y, x] <- z
  }
}


phy <- ape::rphylo(500, 1, 0)

res <- microbenchmark::microbenchmark(func1(),
                                      func2(),
                                      func3())




res
require(tidyverse)
autoplot(res)

found <- c()
for (ntips in 10^(seq(1, 2.5, length.out = 10))) {
  cat(ntips, "\n")
  phy <- ape::rphylo(ntips, 1, 0)
  res <- microbenchmark::microbenchmark(treestats::laplacian_spectrum(phy, 0),
                                        treestats::laplacian_spectrum(phy, 1),
                                        treestats::laplacian_spectrum(phy, 2),
                                        times = 10)
  xvals <- as.numeric(res$expr)
  yvals <- as.numeric(res$time)
  to_add <- cbind(ntips, xvals, yvals)
  found <- rbind(found, to_add)
}

colnames(found) <- c("tree_size", "model", "time")
found <- as_tibble(found)
ggplot(found, aes(x = as.factor(tree_size), y = time, fill = as.factor(model))) +
  geom_boxplot()+
  scale_y_log10()

found %>%
  group_by(tree_size, model) %>%
  summarise("mean_time" = median(time)) %>%
  ggplot(aes(x = tree_size, y = mean_time, col = as.factor(model))) +
    geom_line() +
    scale_y_log10()




if (1 == 2) {

a2 <- treestats::laplacian_spectrum(phy, dens_n = 2000)

profvis::profvis(
for (r in 1:10) {
  treestats::laplacian_spectrum(phy, use_cpp = TRUE)
}
)

}
