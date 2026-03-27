test_that("compute_aic.lm matches stats::AIC", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  result <- compute_aic(fit)
  expect_equal(as.numeric(result), stats::AIC(fit))
  expect_equal(attr(result, "fit_class"), "lm")
  expect_equal(attr(result, "k"), length(coef(fit)))
  expect_equal(attr(result, "criterion"), "AIC")
})

test_that("compute_bic.lm matches stats::BIC", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  result <- compute_bic(fit)
  expect_equal(as.numeric(result), stats::BIC(fit))
  expect_equal(attr(result, "fit_class"), "lm")
})

test_that("compute_aicc.lm is AIC + correction term", {
  fit  <- lm(mpg ~ wt + hp, data = mtcars)
  n    <- nobs(fit)
  k    <- length(coef(fit))
  corr <- 2 * k * (k + 1) / (n - k - 1)
  expect_equal(as.numeric(compute_aicc(fit)), stats::AIC(fit) + corr)
  expect_equal(attr(compute_aicc(fit), "criterion"), "AICc")
})

test_that("compute_aicc > compute_aic for small n", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_gt(as.numeric(compute_aicc(fit)), as.numeric(compute_aic(fit)))
})

test_that("compute_ebic.lm infers P when not supplied", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  result <- compute_ebic(fit)
  # P inferred as length(coef(fit)) - 1 = 2
  expect_equal(attr(result, "P"), 2L)
})

test_that("compute_ebic.lm >= compute_bic.lm", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_gte(as.numeric(compute_ebic(fit, P = 10)), as.numeric(compute_bic(fit)))
  expect_equal(attr(compute_ebic(fit, P = 10), "criterion"), "EBIC")
  expect_true(!is.null(attr(compute_ebic(fit, P = 10), "gamma")))
  expect_equal(attr(compute_ebic(fit, P = 10), "P"), 10)
})

test_that("compute_ebic.lm with gamma = 0 equals BIC", {
  fit    <- lm(mpg ~ wt + hp, data = mtcars)
  result <- compute_ebic(fit, P = 10, gamma = 0)
  expect_equal(as.numeric(result), as.numeric(compute_bic(fit)))
})

test_that("compute_ebic.lm with custom gamma function", {
  fit    <- lm(mpg ~ wt + hp, data = mtcars)
  result <- compute_ebic(fit, P = 10, gamma = function(P, k, n) 0)
  expect_equal(as.numeric(result), as.numeric(compute_bic(fit)))
})

test_that("compute_rbic.lm requires P_index", {
  fit <- lm(mpg ~ wt + hp, data = mtcars)
  expect_error(compute_rbic(fit), "'P_index'")
})

test_that("compute_rbic.lm uses P_index correctly", {
  fit     <- lm(mpg ~ wt + hp + cyl, data = mtcars)
  P_index <- list(group1 = c("wt", "hp"), group2 = "cyl")
  result  <- compute_rbic(fit, P_index = P_index)

  expect_equal(attr(result, "fit_class"), "lm")
  expect_equal(attr(result, "criterion"), "RBIC")
  expect_equal(attr(result, "P_index"), P_index)
  # Default gammafn returns a scalar (single gamma for all groups based on total P)
  expect_true(is.numeric(attr(result, "gamma")))
})

test_that("compute_rbic.lm with gamma=0 equals compute_bic.lm", {
  fit     <- lm(mpg ~ wt + hp + cyl, data = mtcars)
  P_index <- list(group1 = c("wt", "hp"), group2 = "cyl")
  result  <- compute_rbic(fit, P_index = P_index, gamma = 0)
  expect_equal(as.numeric(result), as.numeric(compute_bic(fit)))
})

test_that("compute_rbic.lm per_group gamma differs from ebic gamma", {
  fit     <- lm(mpg ~ wt + hp + cyl + disp + drat + qsec, data = mtcars)
  P_index <- list(small = c("wt", "hp"), large = c("cyl", "disp", "drat", "qsec"))
  ebic_res   <- compute_rbic(fit, P_index = P_index, gamma = "ebic")
  pg_res     <- compute_rbic(fit, P_index = P_index, gamma = "per_group")
  # These should generally differ because group sizes are unequal
  expect_false(identical(attr(ebic_res, "gamma"), attr(pg_res, "gamma")))
})

test_that("compute_aic.glm matches stats::AIC", {
  fit <- glm(am ~ wt + hp, data = mtcars, family = binomial)
  expect_equal(as.numeric(compute_aic(fit)), stats::AIC(fit))
  expect_equal(attr(compute_aic(fit), "fit_class"), "glm")
})

test_that("compute_bic.glm matches stats::BIC", {
  fit <- glm(am ~ wt + hp, data = mtcars, family = binomial)
  expect_equal(as.numeric(compute_bic(fit)), stats::BIC(fit))
})

test_that("intercept-only model: RBIC equals BIC when all k_g = 0", {
  fit     <- lm(mpg ~ 1, data = mtcars)
  P_index <- list(group1 = c("wt", "hp"), group2 = "cyl")
  result  <- compute_rbic(fit, P_index = P_index, gamma = 1)
  # All k_g = 0 so log(C(P_g, 0)) = 0 => penalty = 0 => RBIC = BIC
  expect_equal(as.numeric(result), as.numeric(compute_bic(fit)))
})
