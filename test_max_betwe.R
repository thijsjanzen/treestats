poly_text <- c("((t6:1,t2:1)1:0.5,(((t15:1,t8:1,t11:1)1:1,t12:1,t13:1,t10:1)1:1,((t9:1,(t4:1,t14:1)1:1)1:1,t1:1)1:1,(t3:1,t7:1,t5:1)1:1)1:0.5);", #nolint
               "(((t9:1,t14:1)1:1,t4:1,t12:1)1:0.5,(((t6:1,t19:1,t10:1,t2:1)1:1,t8:1,(t5:1,t16:1)1:1)1:1,(((((t11:1,t17:1)1:1,(t7:1,t3:1)1:1)1:1,(t1:1,t18:1)1:1)1:1,t15:1)1:1,t13:1)1:1)1:0.5);", #nolint
               "((t1:2,t2:2,t3:2):1,((t4:1,t5:1):1,(t6:1,t7:1,t8:1,(t9:0.5,t10:0.5):0.5):1):1);") #nolint  ultrametric tree with polytomies

poly_trees <- list()
for (i in seq_along(poly_text)) {
  poly_trees[[i]] <- ape::read.tree(text = poly_text[i])
}

test_func <- function(phy) {
  df <- as.data.frame(cbind(phy$edge,
                            weight = phy$edge.length))
  g <- igraph::graph_from_data_frame(df, directed = FALSE)
  ref_betweenness <- igraph::betweenness(g)
  return(max(ref_betweenness))
}


for (focal_tree in poly_trees) {
  a1 <- treestats::max_betweenness(focal_tree)
  a2 <- test_func(focal_tree)
  #testthat::expect_equal(a1, a2)
  cat(a1, a2, a1 - a2, "\n")
}


test_func <- function(phy) {
  df <- as.data.frame(cbind(phy$edge,
                            weight = phy$edge.length))
  g <- igraph::graph_from_data_frame(df, directed = FALSE)
  ref_betweenness <- igraph::closeness(g)
  return(max(ref_betweenness))
}


for (focal_tree in poly_trees) {
  a1 <- treestats::max_closeness(focal_tree)
  a2 <- test_func(focal_tree)
  #testthat::expect_equal(a1, a2)
  cat(a1, a2, a1 - a2, "\n")
}

