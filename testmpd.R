vx <-ape::read.tree("/Users/thijsjanzen/Downloads/root_branch.txt")



treestats::mean_pair_dist(vx[[1]])
treestats::mean_pair_dist(vx[[2]])
ape::branching.times(vx[[1]])
vx[[1]]$edge.length
vx[[1]]$edge
cbind(vx[[1]]$edge, vx[[1]]$edge.length)
cbind(vx[[2]]$edge, vx[[2]]$edge.length)
a1 <- treestats::calc_all_stats(vx[[1]])
a2 <- treestats::calc_all_stats(vx[[2]])
rbind(a1, a2)
ape::is.binary(vx[[1]])
ape::is.binary(vx[[2]])
tabulate(vx[[1]]$edge)

phy <- vx[[1]]
ntips <- length(phy$tip.label)
if (!is.null(phy$root.edge)) return(TRUE)
tabulate(phy$edge[, 1])[ntips + 1] <= 2

ape::is.rooted(vx[[1]])
ape::is.rooted(vx[[2]])

n <- length(phy$tip.label)
m <- phy$Nnode
dgr <- tabulate(phy$edge, n + m)
ref <- c(rep.int(1L, n), rep.int(3L, m))
## can use identical() as long as tabulate() returns integers
if (ape::is.rooted(phy)) ref[n + 1L] <- 2L
identical(dgr, ref)
