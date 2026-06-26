poly_tree <- ape::read.tree(text = "((t6:1,t2:1)1:0.5,(((t15:1,t8:1,t11:1)1:1,t12:1,t13:1,t10:1)1:1,((t9:1,(t4:1,t14:1)1:1)1:1,t1:1)1:1,(t3:1,t7:1,t5:1)1:1)1:0.5);") #nolint
set.seed(5)
#poly_tree <- ape::rphylo(n = 30, 1, 0)

a1 <- treestats::cherries(poly_tree)
a2 <- treebalance::cherryI(poly_tree)
cat("cherries: ", a1, a2, "\n")

a1 <- treestats::pitchforks(poly_tree)
a2 <- phyloTop::pitchforks(poly_tree)
cat("pitchforks: ", a1, a2, "\n")

a1 <- treestats::blum(poly_tree)
a2 <- treebalance::sShapeI(poly_tree)
cat("blum: ", a1, a2, "\n")

