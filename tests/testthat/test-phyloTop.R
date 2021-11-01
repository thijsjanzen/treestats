context("phyloTop")

test_that("usage", {
  set.seed(42)
  focal_tree <- TreeSim::sim.bd.taxa(n = 100,
                                     numbsim = 1,
                                     lambda = 1, mu = 0)[[1]]

  avglad <- treestats::avgLadder(focal_tree)
  avglad_check <- phyloTop::avgLadder(focal_tree)
  testthat::expect_equal(avglad, avglad_check)

  cherr <- treestats::cherries(focal_tree)
  cherr_check <- phyloTop::cherries(focal_tree)
  testthat::expect_equal(cherr, cherr_check)

  ILnum <- treestats::ILnumber(focal_tree) # nolint
  ILnum_check <- phyloTop::ILnumber(focal_tree) # nolint
  testthat::expect_equal(ILnum, ILnum_check)

  pitch <- treestats::pitchforks(focal_tree)
  pitch_check <- phyloTop::pitchforks(focal_tree)
  testthat::expect_equal(pitch, pitch_check)

  strs <- treestats::stairs(focal_tree)
  strs_check <- phyloTop::stairs(focal_tree)
  testthat::expect_equal(strs, strs_check)
})
