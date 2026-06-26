require(ape)
require(treestats)

poly_tree <- ape::read.tree(text = "((t6:1,t2:1)1:0.5,(((t15:1,t8:1,t11:1)1:1,t12:1,t13:1,t10:1)1:1,((t9:1,(t4:1,t14:1)1:1)1:1,t1:1)1:1,(t3:1,t7:1,t5:1)1:1)1:0.5);")
set.seed(42)
#poly_tree <- ape::rphylo(n = 5, 1, 0)

#poly_tree <- drop.tip(poly_tree, c("t14", "t4", "t9", "t10", "t13", "t12", "t11", "t8", "t15"))

plot(poly_tree)

# a1 <- treestats::b1(poly_tree)
a2 <- treestats::calc_b1_cpp2(as.vector(t(poly_tree$edge)))
a3 <- treebalance::B1I(poly_tree)

cat( a2, a3, "\n")







