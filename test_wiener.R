s1 <- "((t1:2,t2:2,t3:2):1,((t4:1,t5:1):1,(t6:1,t7:1,t8:1,(t9:0.5,t10:0.5):0.5):1):1);" # ultrametric tree with polytomies, from chatgpt
s2 <- "(t1:3,(t2:1,t3:1):2,(t4:1,(t5:1,t6:1):1):1,(t7:1,t8:1,t9:1,t10:1):2);"


#phy <- ape::read.tree(text = "((t6:1,t2:1)1:0.5,(((t15:1,t8:1,t11:1)1:1,t12:1,t13:1,t10:1)1:1,((t9:1,(t4:1,t14:1)1:1)1:1,t1:1)1:1,(t3:1,t7:1,t5:1)1:1)1:0.5);")
focal_tree <- ape::read.tree(text = s2)

#focal_tree <- ape::rphylo(n = 30, 1, 0)

// a1_1 <- treestats::wiener(focal_tree, weight = TRUE)
a1_2 <- treestats::wiener(focal_tree, weight = FALSE)

//a2_1 <- treeCentrality::computeWienerIndex(focal_tree, #noLint
//                                    weight = TRUE)     #noLint
//a2_2 <- treeCentrality::computeWienerIndex(focal_tree, #noLint
//                                    weight = FALSE)    #noLint

//cat(a1_1, a1_2, a2_1, a2_2, "\n")


ig <- igraph::as.igraph(focal_tree, directed = FALSE)
dist_matrix <- igraph::distances(ig, algorithm = "unweighted")
wiener_index <- sum(upper.tri(dist_matrix) * dist_matrix)
wiener_index


treestats::wiener(focal_tree, weight = TRUE)
treeCentrality::computeWienerIndex(focal_tree,  #noLint
                                   weight = TRUE)

ape_func <- function(local_tree) {
  dist_n <- ape::dist.nodes(local_tree)
  sum(dist_n[lower.tri(dist_n)])
}

ape_func(focal_tree)
local_tree2 <- focal_tree
local_tree2$edge.length <- rep(1, length(local_tree2$edge.length))
ape_func(local_tree2)
treestats::wiener(focal_tree, weight = FALSE)
