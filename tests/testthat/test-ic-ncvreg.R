test_that("compute_aic.ncvreg returns vector of correct length", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(42)
  X <- matrix(rnorm(100 * 10), 100, 10)
  y <- X[, 1] * 2 + rnorm(100)
  fit <- ncvreg(X, y, family = "gaussian")

  result <- compute_aic(fit)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "fit_class"), "ncvreg")
  expect_equal(attr(result, "criterion"), "AIC")
})

test_that("compute_bic.ncvreg returns vector of correct length", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(1)
  X <- matrix(rnorm(80 * 8), 80, 8)
  y <- X[, 2] - X[, 5] + rnorm(80)
  fit <- ncvreg(X, y)

  result <- compute_bic(fit)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "criterion"), "BIC")
})

test_that("compute_aicc.ncvreg returns vector of correct length", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(2)
  X <- matrix(rnorm(60 * 5), 60, 5)
  y <- X[, 1] + rnorm(60)
  fit <- ncvreg(X, y)

  result <- compute_aicc(fit)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "criterion"), "AICc")
})

test_that("compute_ebic.ncvreg returns vector >= compute_bic.ncvreg", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(5)
  X <- matrix(rnorm(100 * 10), 100, 10)
  y <- X[, 1] + rnorm(100)
  fit <- ncvreg(X, y)

  ebic <- as.numeric(compute_ebic(fit, P = 10))
  bic  <- as.numeric(compute_bic(fit))
  expect_true(all(ebic >= bic - 1e-10))
  expect_equal(attr(compute_ebic(fit, P = 10), "criterion"), "EBIC")
})

test_that("compute_rbic.ncvreg validates P_index names", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  X   <- matrix(rnorm(60 * 4), 60, 4,
                dimnames = list(NULL, paste0("x", 1:4)))
  y   <- rnorm(60)
  fit <- ncvreg(X, y)
  bad <- list(g1 = c("x1", "x_missing"))
  expect_error(compute_rbic(fit, P_index = bad), "not found in rownames")
})

test_that("compute_rbic.ncvreg returns correct length and attributes", {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(9)
  X   <- matrix(rnorm(100 * 6), 100, 6,
                dimnames = list(NULL, paste0("x", 1:6)))
  y   <- X[, 1] - X[, 4] + rnorm(100)
  fit <- ncvreg(X, y)
  P_index <- list(g1 = paste0("x", 1:3), g2 = paste0("x", 4:6))

  result <- compute_rbic(fit, P_index = P_index)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "fit_class"), "ncvreg")
  expect_equal(attr(result, "criterion"), "RBIC")
})
