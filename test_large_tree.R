av <- ape::read.tree("/Users/thijsjanzen/Downloads/2trees/treeB4.txt")
av2 <- ape::root(av, outgroup = "otu17164", r = TRUE)
treestats::pigot_rho(av2)


av <- ape::read.tree("/Users/thijsjanzen/Downloads/2trees/treeAfter.txt")
av2 <- ape::root(av, outgroup = "otu17164", r = TRUE)
ape::is.rooted(av2)


treestats::mean_pair_dist(av2)
#av3 <- treestats::phylo_to_l(av2)
#treestats::mean_pair_dist(av3)
#treestats::psv(av2)
