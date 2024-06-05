context("hard_tree")

test_that("usage", {

  # Thanks to Fien Strijthaegen for providing a testing tree

  focal_tree <- ape::read.tree(text = "(151:34,(152:2,((((153:10,154:16):1,155:24):2,156:2):6,(157:5,((158:0,(((159:0,160:5):14,161:1):0,((((162:4,163:0):6,164:0):0,(165:22,166:4):3):1,167:13):4):9):3,(168:25,169:3):29):2):0):2):0,170:11,171:4):0;") #nolint

  all_stats <- treestats::calc_all_stats(focal_tree)
  testthat::expect_equal(sum(is.na(all_stats)), 27)

  brts_stats <- treestats::calc_brts_stats(focal_tree)
  testthat::expect_true(is.na(brts_stats$gamma))

  bal_stats <- treestats::calc_balance_stats(focal_tree)
  testthat::expect_equal(sum(is.na(bal_stats)), 22)
})
