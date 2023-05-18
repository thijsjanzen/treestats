focal_tree <- ape::rphylo(n = 1000, birth = 1, death = 0.0)

c1 <- treestats::colless(focal_tree)
c2 <- treestats::colless_test(focal_tree)
testthat::expect_equal(c1, c2)

ewc1 <- treestats::ew_colless(focal_tree)
ewc2 <- treestats::ew_colless_test(focal_tree)
testthat::expect_equal(ewc1, ewc2)

s1 <- treestats::stairs(focal_tree)
s2 <- treestats::stairs_test(focal_tree)
testthat::expect_equal(s1, s2)

measure_time <- FALSE
if (measure_time) {

    microbenchmark::microbenchmark(treestats::colless(focal_tree),
                                treestats::colless_test(focal_tree)
    )

    microbenchmark::microbenchmark(treestats::ew_colless(focal_tree),
                                treestats::ew_colless_test(focal_tree)
    )

    microbenchmark::microbenchmark(treestats::stairs(focal_tree),
                                treestats::stairs_test(focal_tree)
    )

}
