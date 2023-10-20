test_that("Package style", {
  testthat::skip_if_not_installed("lintr")
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  lintr::expect_lint_free(
    linters = lintr::linters_with_defaults(object_name_linter = NULL,
                                           indentation_linter = NULL))
})
