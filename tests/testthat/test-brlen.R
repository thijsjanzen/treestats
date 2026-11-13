context("branch_length_stats")

test_that("usage", {
  focal_tree <- ape::read.tree(text = "((t2:1,t1:1):2,(t4:2,t3:2):1);")

  m_br <- treestats::mean_branch_length(focal_tree)
  testthat::expect_equal(m_br, 1.5)

  v_br <- treestats::var_branch_length(focal_tree)
  testthat::expect_equal(v_br, 0.3)

  m_br_ext <- treestats::mean_branch_length_ext(focal_tree)
  testthat::expect_equal(m_br_ext, 1.5)

  m_br_int <- treestats::mean_branch_length_int(focal_tree)
  testthat::expect_equal(m_br_ext, 1.5)

  v_br_ext <-   treestats::var_branch_length_ext(focal_tree)
  testthat::expect_equal(v_br_ext, 1 / 3)

  v_br_int <-   treestats::var_branch_length_int(focal_tree)
  testthat::expect_equal(v_br_int, 1 / 2)


  ltab <- treestats::phylo_to_l(focal_tree)

  m_br <- treestats::mean_branch_length(ltab)
  testthat::expect_equal(m_br, 1.5)

  v_br <- treestats::var_branch_length(ltab)
  testthat::expect_equal(v_br, 0.3)

  m_br_ext <- treestats::mean_branch_length_ext(ltab)
  testthat::expect_equal(m_br_ext, 1.5)

  m_br_int <- treestats::mean_branch_length_int(ltab)
  testthat::expect_equal(m_br_ext, 1.5)

  v_br_ext <-   treestats::var_branch_length_ext(ltab)
  testthat::expect_equal(v_br_ext, 1 / 3)

  v_br_int <-   treestats::var_branch_length_int(ltab)
  testthat::expect_equal(v_br_int, 1 / 2)
})

test_that("wrong_object", {
  test_abuse(treestats::mean_branch_length)
  test_abuse(treestats::mean_branch_length_ext)
  test_abuse(treestats::mean_branch_length_int)
  test_abuse(treestats::var_branch_length)
  test_abuse(treestats::var_branch_length_ext)
  test_abuse(treestats::var_branch_length_int)
})
