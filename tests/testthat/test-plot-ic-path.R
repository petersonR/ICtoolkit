# Tests for plot_ic_path()
#
# Plotting calls are wrapped in a null graphics device so the test run does
# not create stray Rplots.pdf files.

make_glmnet_fit <- function(n = 100, p = 20, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(NULL, paste0("x", seq_len(p))))
  y <- X[, 1] * 2 - X[, 3] + rnorm(n)
  glmnet::glmnet(X, y)
}

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

test_that("plot_ic_path returns ic, opt_lambda, and lambda invisibly", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()

  res <- with_null_device(plot_ic_path(fit, criteria = c("AIC", "BIC")))

  expect_type(res, "list")
  expect_named(res, c("ic", "opt_lambda", "lambda"))
  expect_named(res$ic, c("AIC", "BIC"))
  expect_length(res$ic$AIC, length(fit$lambda))
  expect_named(res$opt_lambda, c("AIC", "BIC"))
  expect_length(res$opt_lambda, 2)
  expect_equal(res$lambda, fit$lambda)
})

test_that("plot_ic_path accepts user-supplied ylim (regression: arg collision)", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()

  # Previously errored: "formal argument 'ylim' is matched by multiple
  # actual arguments" because ylim was both computed internally and
  # forwarded through `...`.
  expect_no_error(
    with_null_device(plot_ic_path(fit, criteria = c("AIC", "BIC"),
                                  ylim = c(0, 600)))
  )
})

test_that("plot_ic_path accepts user-supplied xlim", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()

  expect_no_error(
    with_null_device(plot_ic_path(fit, criteria = c("AIC", "BIC"),
                                  xlim = c(0, 4)))
  )
})

test_that("plot_ic_path accepts legend_pos as a keyword", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()

  expect_no_error(
    with_null_device(plot_ic_path(fit, criteria = c("AIC", "BIC"),
                                  legend_pos = "bottomleft"))
  )
})

test_that("plot_ic_path accepts legend_pos as numeric coordinates", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()

  expect_no_error(
    with_null_device(plot_ic_path(fit, criteria = c("AIC", "BIC"),
                                  legend_pos = c(1, 450)))
  )
})

test_that("plot_ic_path accepts ylim, xlim, and legend_pos together", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()

  res <- with_null_device(
    plot_ic_path(fit, criteria = c("AIC", "BIC", "EBIC"), P = 20,
                 ylim = c(0, 700), xlim = c(0, 5),
                 legend_pos = "topleft")
  )
  expect_named(res$ic, c("AIC", "BIC", "EBIC"))
})

test_that("plot_ic_path works with a single criterion and legend = FALSE", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()

  expect_no_error(
    with_null_device(plot_ic_path(fit, criteria = "BIC",
                                  legend = FALSE, mark_min = FALSE))
  )
})

test_that("plot_ic_path supports ncvreg fits", {
  skip_if_not_installed("ncvreg")
  set.seed(2)
  X <- matrix(rnorm(100 * 8), 100, 8,
              dimnames = list(NULL, paste0("x", 1:8)))
  y <- X[, 1] - X[, 2] + rnorm(100)
  fit <- ncvreg::ncvreg(X, y)

  res <- with_null_device(
    plot_ic_path(fit, criteria = c("AIC", "BIC"),
                 ylim = c(0, 600), legend_pos = "bottomright")
  )
  expect_named(res$ic, c("AIC", "BIC"))
})

test_that("plot_ic_path rejects non-glmnet/ncvreg fits", {
  fit <- lm(mpg ~ wt, data = mtcars)
  expect_error(plot_ic_path(fit), "glmnet, ncvreg, or ncvsurv")
})

test_that("plot_ic_path rejects unknown criteria", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()
  expect_error(
    with_null_device(plot_ic_path(fit, criteria = c("AIC", "NOPE"))),
    "Unknown criteria"
  )
})

test_that("plot_ic_path requires P_index for RBIC", {
  skip_if_not_installed("glmnet")
  fit <- make_glmnet_fit()
  expect_error(
    with_null_device(plot_ic_path(fit, criteria = c("BIC", "RBIC"))),
    "'P_index'"
  )
})
