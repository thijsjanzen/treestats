set.seed(5)
focal_tree <- ape::rphylo(n = 300, 1, 0)
plot(focal_tree)
ed <- focal_tree$edge
el <- focal_tree$edge.length

dat <- cbind(ed, el)

dat2 <- dat[order(dat[, 1]), ]

stat <- rep(0, max(dat2[, 1]))
total_diff <- rep(0, max(dat2[, 1]))
for (x in rev(unique(dat2[, 1]))) {
  a <- subset(dat2, dat2[, 1] == x)

  L <- a[1, 3] + stat[a[1, 2]]
  R <- a[2, 3] + stat[a[2, 2]]
  total_diff[x] <- abs(L - R)

  stat[x] <- stat[x] + L + R
}
sum(total_diff)
treestats::colless_branch(focal_tree)

# root
# L = 1.078
# R = 1.078 - 0.66 + 2 * 0.66 = 1.738
# diff = 0.66

# internal node
# L = 0.66
# R = 0.66
# diff = 0.0
# total stat = 0.66

