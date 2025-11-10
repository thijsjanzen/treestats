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

test_that("unrooted", {
  set.seed(42)
  focal_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
  a1 <- treestats::mean_branch_length_ext(focal_tree)
  a2 <- treestats::mean_branch_length_ext(ape::unroot(focal_tree))
  testthat::expect_equal(a1, a2)
})

test_that("wrong_object", {
  testthat::expect_error(
    treestats::var_branch_length(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::var_branch_length(list()),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::var_branch_length_int(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::var_branch_length_int(list()),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::var_branch_length_ext(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::var_branch_length_ext(list()),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mean_branch_length_ext(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mean_branch_length_ext(list()),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mean_branch_length_int(10),
    "input object has to be phylo or ltable"
  )

  testthat::expect_error(
    treestats::mean_branch_length_int(list()),
    "input object has to be phylo or ltable"
  )
})
