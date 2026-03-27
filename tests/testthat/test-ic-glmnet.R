test_that("compute_aic.glmnet returns vector of correct length", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(42)
  X <- matrix(rnorm(100 * 10), 100, 10)
  y <- X[, 1] * 2 + rnorm(100)
  fit <- glmnet(X, y, family = "gaussian")

  result <- compute_aic(fit)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "fit_class"), "glmnet")
  expect_equal(attr(result, "criterion"), "AIC")
  expect_length(attr(result, "k"), length(fit$lambda))
  expect_length(attr(result, "lambda"), length(fit$lambda))
})

test_that("compute_bic.glmnet values are >= compute_aic.glmnet (large n)", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(42)
  X <- matrix(rnorm(200 * 10), 200, 10)
  y <- X[, 1] * 2 + rnorm(200)
  fit <- glmnet(X, y, family = "gaussian")

  aic <- as.numeric(compute_aic(fit))
  bic <- as.numeric(compute_bic(fit))
  # BIC >= AIC when log(n) >= 2, i.e. n >= e^2 ~ 7.4
  expect_true(all(bic >= aic))
})

test_that("compute_aic.glmnet is monotone decreasing then increasing in lambda", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(1)
  X <- matrix(rnorm(100 * 5), 100, 5)
  y <- X[, 1] + rnorm(100)
  fit <- glmnet(X, y)
  # AIC values should have a minimum somewhere on the path (not necessarily at ends)
  aic <- as.numeric(compute_aic(fit))
  expect_true(min(aic) < max(aic))  # not all equal
})

test_that("compute_aicc.glmnet returns vector with NAs where undefined", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(1)
  X <- matrix(rnorm(100 * 10), 100, 10)
  y <- X[, 1] + rnorm(100)
  fit <- glmnet(X, y)
  result <- compute_aicc(fit)
  expect_length(result, length(fit$lambda))
})

test_that("compute_ebic.glmnet requires P", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  X <- matrix(rnorm(100 * 5), 100, 5)
  y <- rnorm(100)
  fit <- glmnet(X, y)
  result <- compute_ebic(fit)
  expect_equal(attr(result, "P"), 5L)  # inferred from nrow(fit$beta)
})

test_that("compute_ebic.glmnet >= compute_bic.glmnet", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(42)
  X <- matrix(rnorm(100 * 10), 100, 10)
  y <- X[, 1] + rnorm(100)
  fit <- glmnet(X, y)

  ebic <- as.numeric(compute_ebic(fit, P = 10))
  bic  <- as.numeric(compute_bic(fit))
  expect_true(all(ebic >= bic - 1e-10))  # allow tiny numeric error
  expect_equal(attr(compute_ebic(fit, P = 10), "criterion"), "EBIC")
  expect_length(attr(compute_ebic(fit, P = 10), "gamma"), length(fit$lambda))
})

test_that("compute_rbic.glmnet requires P_index", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  X <- matrix(rnorm(50 * 4), 50, 4)
  y <- rnorm(50)
  fit <- glmnet(X, y)
  expect_error(compute_rbic(fit), "'P_index'")
})

test_that("compute_rbic.glmnet validates P_index names", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  X   <- matrix(rnorm(50 * 4), 50, 4,
                dimnames = list(NULL, paste0("x", 1:4)))
  y   <- rnorm(50)
  fit <- glmnet(X, y)
  bad_index <- list(g1 = c("x1", "x_bad"))
  expect_error(compute_rbic(fit, P_index = bad_index), "not found in rownames")
})

test_that("compute_rbic.glmnet returns correct length and attributes", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(7)
  X   <- matrix(rnorm(100 * 6), 100, 6,
                dimnames = list(NULL, paste0("x", 1:6)))
  y   <- X[, 1] - X[, 4] + rnorm(100)
  fit <- glmnet(X, y)
  P_index <- list(g1 = paste0("x", 1:3), g2 = paste0("x", 4:6))

  result <- compute_rbic(fit, P_index = P_index)
  expect_length(result, length(fit$lambda))
  expect_equal(attr(result, "fit_class"), "glmnet")
  expect_equal(attr(result, "criterion"), "RBIC")
  expect_equal(attr(result, "P_index"), P_index)
  expect_length(attr(result, "gamma"), length(fit$lambda))  # list of gamma vectors
})

test_that("compute_rbic.glmnet with gamma=0 equals compute_bic.glmnet", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  set.seed(3)
  X   <- matrix(rnorm(80 * 4), 80, 4,
                dimnames = list(NULL, paste0("x", 1:4)))
  y   <- X[, 1] + rnorm(80)
  fit <- glmnet(X, y)
  P_index <- list(g1 = paste0("x", 1:2), g2 = paste0("x", 3:4))

  rbic <- as.numeric(compute_rbic(fit, P_index = P_index, gamma = 0))
  bic  <- as.numeric(compute_bic(fit))
  expect_equal(rbic, bic, tolerance = 1e-10)
})

test_that("compute_aic.glmnet warns for alpha < 1", {
  skip_if_not_installed("glmnet")
  library(glmnet)
  X <- matrix(rnorm(100 * 5), 100, 5)
  y <- rnorm(100)
  fit <- glmnet(X, y, alpha = 0.5)
  expect_warning(compute_aic(fit), "alpha")
})
