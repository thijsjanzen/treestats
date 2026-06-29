context("polytomies")

test_that("usage", {

  # Thanks to Fien Strijthaegen for providing a testing tree

  poly_trees <- c("(151:34,(152:2,((((153:10,154:16):1,155:24):2,156:2):6,(157:5,((158:0,(((159:0,160:5):14,161:1):0,((((162:4,163:0):6,164:0):0,(165:22,166:4):3):1,167:13):4):9):3,(168:25,169:3):29):2):0):2):0,170:11,171:4):0;",
                  "((t6:1,t2:1)1:0.5,(((t15:1,t8:1,t11:1)1:1,t12:1,t13:1,t10:1)1:1,((t9:1,(t4:1,t14:1)1:1)1:1,t1:1)1:1,(t3:1,t7:1,t5:1)1:1)1:0.5);",
                  "(((t9:1,t14:1)1:1,t4:1,t12:1)1:0.5,(((t6:1,t19:1,t10:1,t2:1)1:1,t8:1,(t5:1,t16:1)1:1)1:1,(((((t11:1,t17:1)1:1,(t7:1,t3:1)1:1)1:1,(t1:1,t18:1)1:1)1:1,t15:1)1:1,t13:1)1:1)1:0.5);") # nolint
  for (poly_text in poly_trees) {
    poly_tree <- ape::read.tree(text = poly_text)
    all_stats <- treestats::calc_all_stats(poly_tree)

    na_stats <- names(all_stats)[which(is.na(all_stats))]
    non_na_stats <- names(all_stats)[which(!is.na(all_stats))]

    names_na_stats <- c("avg_ladder", "beta",
                        "colless", "colless_corr", "colless_quad",
                        "diameter", "double_cherries",
                        "eigen_centrality", "eigen_centralityW",
                        "ew_colless", "four_prong", "gamma",
                        "i_stat", "imbalance_steps", "j_one",
                        "max_adj", "max_betweenness",
                        "max_ladder", "max_laplace", "min_adj",
                        "min_laplace", "rogers",
                        "stairs", "stairs2", "symmetry_nodes",
                        "wiener",
                        "root_imbalance",
                        "mean_inv_branch_dist", # TODO: make this doable
                        "nltt_base") # TODO: nLTT is NA because of ultrametric!

    a <- na_stats %in% names_na_stats
    b <- names_na_stats %in% na_stats
 # disable tests for now (29-06-2026)
    #   testthat::expect_all_true(a)
 #    testthat::expect_all_true(b)

    # test phylo_to_l and rebase_ltable on polytomies:
    testthat::expect_silent(
      ltab  <- treestats::phylo_to_l(poly_tree))
    testthat::expect_silent(
      ltab2 <- treestats::rebase_ltable(ltab)
    )
  }
})

test_that("manual trees", {
  # these trees were provided by Luise Haeuser.

  poly_tree <- ape::read.tree(text = "(((t9:1,t14:1)1:1,t4:1,t12:1)1:0.5,(((t6:1,t19:1,t10:1,t2:1)1:1,t8:1,(t5:1,t16:1)1:1)1:1,(((((t11:1,t17:1)1:1,(t7:1,t3:1)1:1)1:1,(t1:1,t18:1)1:1)1:1,t15:1)1:1,t13:1)1:1)1:0.5);") #nolint
  il_val <- treestats::ILnumber(poly_tree)
  testthat::expect_equal(il_val, 3) # reference value through solving by eye

  b1_val <- treestats::b1(poly_tree)
  b1_ref <- treebalance::B1I(poly_tree)
  testthat::expect_equal(b1_val, b1_ref)

  poly_tree <- ape::read.tree(text = "((t6:1,t2:1)1:0.5,(((t15:1,t8:1,t11:1)1:1,t12:1,t13:1,t10:1)1:1,((t9:1,(t4:1,t14:1)1:1)1:1,t1:1)1:1,(t3:1,t7:1,t5:1)1:1)1:0.5);") #nolint
  tip1 <- treestats::tot_internal_path(poly_tree)
  tip2 <- treebalance::totIntPathLen(poly_tree)
  testthat::expect_equal(tip1, tip2)

  b1_val <- treestats::b1(poly_tree)
  b1_ref <- treebalance::B1I(poly_tree)
  testthat::expect_equal(b1_val, b1_ref)
})
