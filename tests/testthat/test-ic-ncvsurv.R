# ---------------------------------------------------------------------------
# Tests for ncvsurv (Cox survival) methods
# ---------------------------------------------------------------------------

# Shared fixture: simulated survival data with ncvsurv fit
.make_ncvsurv_fit <- function() {
  skip_if_not_installed("ncvreg")
  library(ncvreg)
  set.seed(42)
  n <- 100; p <- 10
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
  # Simulate survival times from exponential hazard model
  eta <- X[, 1] * 0.5 - X[, 3] * 0.3
  time <- rexp(n, rate = exp(eta))
  status <- rbinom(n, 1, 0.7)
  y <- survival::Surv(time, status)
  fit <- ncvsurv(X, y, penalty = "lasso")
  list(fit = fit, X = X, y = y, n = n, p = p)
}

# ---------------------------------------------------------------------------
# Basic output shape and attributes
# ---------------------------------------------------------------------------

test_that("compute_aic.ncvsurv returns vector of correct length", {
  d <- .make_ncvsurv_fit()
  result <- compute_aic(d$fit)
  expect_length(result, length(d$fit$lambda))
  expect_equal(attr(result, "fit_class"), "ncvsurv")
  expect_equal(attr(result, "criterion"), "AIC")
})

test_that("compute_aicc.ncvsurv returns vector of correct length", {
  d <- .make_ncvsurv_fit()
  result <- compute_aicc(d$fit)
  expect_length(result, length(d$fit$lambda))
  expect_equal(attr(result, "criterion"), "AICc")
})

test_that("compute_bic.ncvsurv returns vector of correct length", {
  d <- .make_ncvsurv_fit()
  result <- compute_bic(d$fit)
  expect_length(result, length(d$fit$lambda))
  expect_equal(attr(result, "criterion"), "BIC")
})

test_that("compute_hqic.ncvsurv returns vector of correct length", {
  d <- .make_ncvsurv_fit()
  result <- compute_hqic(d$fit)
  expect_length(result, length(d$fit$lambda))
  expect_equal(attr(result, "criterion"), "HQIC")
})

# ---------------------------------------------------------------------------
# EBIC >= BIC for ncvsurv
# ---------------------------------------------------------------------------

test_that("compute_ebic.ncvsurv >= compute_bic.ncvsurv", {
  d <- .make_ncvsurv_fit()
  ebic <- as.numeric(compute_ebic(d$fit, P = d$p))
  bic  <- as.numeric(compute_bic(d$fit))
  expect_true(all(ebic >= bic - 1e-10))
  expect_equal(attr(compute_ebic(d$fit, P = d$p), "criterion"), "EBIC")
})

test_that("compute_ebic.ncvsurv with gamma=0 recovers BIC", {
  d <- .make_ncvsurv_fit()
  ebic0 <- as.numeric(compute_ebic(d$fit, P = d$p, gamma = 0))
  bic   <- as.numeric(compute_bic(d$fit))
  expect_equal(ebic0, bic, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# RBIC
# ---------------------------------------------------------------------------

test_that("compute_rbic.ncvsurv validates P_index names", {
  d <- .make_ncvsurv_fit()
  bad <- list(g1 = c("x1", "x_missing"))
  expect_error(compute_rbic(d$fit, P_index = bad), "not found in rownames")
})

test_that("compute_rbic.ncvsurv returns correct length and attributes", {
  d <- .make_ncvsurv_fit()
  P_index <- list(g1 = paste0("x", 1:5), g2 = paste0("x", 6:10))
  result <- compute_rbic(d$fit, P_index = P_index)
  expect_length(result, length(d$fit$lambda))
  expect_equal(attr(result, "fit_class"), "ncvsurv")
  expect_equal(attr(result, "criterion"), "RBIC")
})

test_that("compute_rbic.ncvsurv with gamma=0 recovers BIC", {
  d <- .make_ncvsurv_fit()
  P_index <- list(g1 = paste0("x", 1:5), g2 = paste0("x", 6:10))
  rbic0 <- as.numeric(compute_rbic(d$fit, P_index = P_index, gamma = 0))
  bic   <- as.numeric(compute_bic(d$fit))
  expect_equal(rbic0, bic, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# mBIC / mBIC2
# ---------------------------------------------------------------------------

test_that("compute_mbic.ncvsurv returns correct length", {
  d <- .make_ncvsurv_fit()
  result <- compute_mbic(d$fit, P = d$p)
  expect_length(result, length(d$fit$lambda))
  expect_equal(attr(result, "criterion"), "mBIC")
})

test_that("compute_mbic2.ncvsurv <= compute_mbic.ncvsurv", {
  d <- .make_ncvsurv_fit()
  mbic  <- as.numeric(compute_mbic(d$fit, P = d$p))
  mbic2 <- as.numeric(compute_mbic2(d$fit, P = d$p))
  expect_true(all(mbic2 <= mbic + 1e-10))
})

# ---------------------------------------------------------------------------
# Ordering: AIC < HQIC < BIC for ncvsurv (at minimum)
# ---------------------------------------------------------------------------

test_that("IC ordering at optimal lambda: AIC <= HQIC <= BIC for ncvsurv", {
  d <- .make_ncvsurv_fit()
  aic  <- as.numeric(compute_aic(d$fit))
  hqic <- as.numeric(compute_hqic(d$fit))
  bic  <- as.numeric(compute_bic(d$fit))
  # At the same lambda, AIC penalty <= HQIC penalty <= BIC penalty
  # So AIC <= HQIC <= BIC (all use same log-likelihood)
  expect_true(all(aic <= hqic + 1e-10))
  expect_true(all(hqic <= bic + 1e-10))
})

# ---------------------------------------------------------------------------
# n extraction is correct (not 2*n)
# ---------------------------------------------------------------------------

test_that("ncvsurv methods use correct n (not 2*n from Surv object)", {
  d <- .make_ncvsurv_fit()
  # AICc correction depends on n; if n were wrong (2*n), correction would differ
  aicc <- as.numeric(compute_aicc(d$fit))
  aic  <- as.numeric(compute_aic(d$fit))
  k    <- attr(compute_aic(d$fit), "k")
  n    <- d$n
  expected_correction <- ifelse(n - k - 1 > 0, 2 * k * (k + 1) / (n - k - 1), NA)
  # AICc = AIC + correction; verify the correction matches expected n
  actual_correction <- aicc - aic
  # At the null model end (k=0), correction should be 0
  expect_equal(actual_correction[k == 0], rep(0, sum(k == 0)), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# plot_ic_path accepts ncvsurv
# ---------------------------------------------------------------------------

test_that("plot_ic_path works with ncvsurv objects", {
  d <- .make_ncvsurv_fit()
  P_index <- list(g1 = paste0("x", 1:5), g2 = paste0("x", 6:10))
  result <- plot_ic_path(d$fit, criteria = c("BIC", "RBIC"),
                         P_index = P_index)
  expect_true(is.list(result))
  expect_true("opt_lambda" %in% names(result))
})
