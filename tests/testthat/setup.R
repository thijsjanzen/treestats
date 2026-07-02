poly_text <- c("((t6:1,t2:1)1:0.5,(((t15:1,t8:1,t11:1)1:1,t12:1,t13:1,t10:1)1:1,((t9:1,(t4:1,t14:1)1:1)1:1,t1:1)1:1,(t3:1,t7:1,t5:1)1:1)1:0.5);", #nolint
                "(((t9:1,t14:1)1:1,t4:1,t12:1)1:0.5,(((t6:1,t19:1,t10:1,t2:1)1:1,t8:1,(t5:1,t16:1)1:1)1:1,(((((t11:1,t17:1)1:1,(t7:1,t3:1)1:1)1:1,(t1:1,t18:1)1:1)1:1,t15:1)1:1,t13:1)1:1)1:0.5);", #nolint
                "((t1:2,t2:2,t3:2):1,((t4:1,t5:1):1,(t6:1,t7:1,t8:1,(t9:0.5,t10:0.5):0.5):1):1);") #nolint  ultrametric tree with polytomies

poly_trees <- list()
for (i in seq_along(poly_text)) {
  poly_trees[[i]] <- ape::read.tree(text = poly_text[i])
}

set.seed(42)
extant_tree   <- ape::rphylo(n = 100, birth = 1, death = 0)
complete_tree <- ape::rphylo(n = 100, birth = 1, death = 0.2, fossils = TRUE)

test_polytomies <- function(treestats_func, ref_func, tol = NA) {
  for (focal_tree in poly_trees) {
    a1 <- try_stat(focal_tree, treestats_func)
    a2 <- try_stat(focal_tree, ref_func)
    if (is.na(tol)) {
      testthat::expect_equal(a1, a2)
    } else {
      testthat::expect_equal(a1, a2, tolerance = tol)
    }
  }
}

standard_test <- function(treestats_func,
                          ref_func,
                          drop_complete = FALSE) {
  a1 <- treestats_func(extant_tree)
  a2 <- ref_func(extant_tree)
  testthat::expect_equal(a1, a2)

  ltab <- treestats::phylo_to_l(extant_tree)
  testthat::expect_equal(treestats_func(extant_tree),
                         treestats_func(ltab))

  if (drop_complete == FALSE) {
    # with extinct species:
    a1 <- treestats_func(complete_tree)
    a2 <- ref_func(complete_tree)
    testthat::expect_equal(a1, a2)

    ltab <- treestats::phylo_to_l(complete_tree)
    testthat::expect_equal(treestats_func(complete_tree),
                           treestats_func(ltab))
  }
}
