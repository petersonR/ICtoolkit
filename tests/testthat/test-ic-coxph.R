# ---------------------------------------------------------------------------
# Tests for coxph methods (scalar IC) and ic_step with coxph
# ---------------------------------------------------------------------------

# Shared fixture: lung cancer data from survival package
.make_coxph_data <- function() {
  skip_if_not_installed("survival")
  library(survival)
  d <- lung
  d <- d[complete.cases(d[, c("time", "status", "age", "sex", "ph.ecog",
                                "ph.karno", "wt.loss")]), ]
  d$sex <- factor(d$sex)
  d
}

# ---------------------------------------------------------------------------
# Basic scalar IC for coxph
# ---------------------------------------------------------------------------

test_that("compute_aic.coxph returns scalar with correct attributes", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  result <- compute_aic(fit)
  expect_length(result, 1)
  expect_equal(attr(result, "fit_class"), "coxph")
  expect_equal(attr(result, "criterion"), "AIC")
  expect_equal(attr(result, "k"), 2L)
  # Should match stats::AIC
  expect_equal(as.numeric(result), stats::AIC(fit), tolerance = 1e-10)
})

test_that("compute_aicc.coxph returns scalar", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  result <- compute_aicc(fit)
  expect_length(result, 1)
  expect_equal(attr(result, "criterion"), "AICc")
  # AICc >= AIC
  expect_true(as.numeric(result) >= as.numeric(compute_aic(fit)) - 1e-10)
})

test_that("compute_bic.coxph returns scalar matching stats::BIC", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  result <- compute_bic(fit)
  expect_length(result, 1)
  expect_equal(attr(result, "criterion"), "BIC")
  expect_equal(as.numeric(result), stats::BIC(fit), tolerance = 1e-10)
})

test_that("compute_hqic.coxph returns scalar", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  result <- compute_hqic(fit)
  expect_length(result, 1)
  expect_equal(attr(result, "criterion"), "HQIC")
})

test_that("IC ordering for coxph: AIC <= HQIC <= BIC", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                          data = d)
  aic  <- as.numeric(compute_aic(fit))
  hqic <- as.numeric(compute_hqic(fit))
  bic  <- as.numeric(compute_bic(fit))
  expect_true(aic <= hqic + 1e-10)
  expect_true(hqic <= bic + 1e-10)
})

# ---------------------------------------------------------------------------
# EBIC / RBIC for coxph
# ---------------------------------------------------------------------------

test_that("compute_ebic.coxph >= compute_bic.coxph", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                          data = d)
  ebic <- as.numeric(compute_ebic(fit, P = 5))
  bic  <- as.numeric(compute_bic(fit))
  expect_true(ebic >= bic - 1e-10)
})

test_that("compute_ebic.coxph with gamma=0 recovers BIC", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  ebic0 <- as.numeric(compute_ebic(fit, P = 5, gamma = 0))
  bic   <- as.numeric(compute_bic(fit))
  expect_equal(ebic0, bic, tolerance = 1e-10)
})

test_that("compute_rbic.coxph returns scalar with correct attributes", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                          data = d)
  P_index <- list(demographic = c("age", "sex"), clinical = "ph.ecog")
  result <- compute_rbic(fit, P_index = P_index)
  expect_length(result, 1)
  expect_equal(attr(result, "criterion"), "RBIC")
})

test_that("compute_rbic.coxph with gamma=0 recovers BIC", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  P_index <- list(g1 = "age", g2 = "sex")
  rbic0 <- as.numeric(compute_rbic(fit, P_index = P_index, gamma = 0))
  bic   <- as.numeric(compute_bic(fit))
  expect_equal(rbic0, bic, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# mBIC / mBIC2 for coxph
# ---------------------------------------------------------------------------

test_that("compute_mbic.coxph returns scalar", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                          data = d)
  result <- compute_mbic(fit, P = 5)
  expect_length(result, 1)
  expect_equal(attr(result, "criterion"), "mBIC")
})

test_that("compute_mbic2.coxph <= compute_mbic.coxph", {
  d <- .make_coxph_data()
  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex + ph.ecog,
                          data = d)
  mbic  <- as.numeric(compute_mbic(fit, P = 5))
  mbic2 <- as.numeric(compute_mbic2(fit, P = 5))
  expect_true(mbic2 <= mbic + 1e-10)
})

# ---------------------------------------------------------------------------
# ic_step with coxph
# ---------------------------------------------------------------------------

test_that("ic_step forward BIC works with coxph", {
  d <- .make_coxph_data()
  fit0 <- survival::coxph(survival::Surv(time, status) ~ 1, data = d)
  scope_upper <- survival::Surv(time, status) ~ age + sex + ph.ecog +
    ph.karno + wt.loss
  result <- ic_step(fit0,
                    scope = list(lower = ~ 1, upper = scope_upper),
                    direction = "forward",
                    criterion = "BIC",
                    trace = 0)
  expect_s3_class(result, "coxph")
  path <- attr(result, "step_path")
  expect_s3_class(path, "data.frame")
  expect_true("BIC" %in% names(path))
})

test_that("ic_step backward AIC works with coxph", {
  d <- .make_coxph_data()
  fit_full <- survival::coxph(survival::Surv(time, status) ~ age + sex +
    ph.ecog + ph.karno + wt.loss, data = d)
  result <- ic_step(fit_full, direction = "backward", criterion = "AIC",
                    trace = 0)
  expect_s3_class(result, "coxph")
})

test_that("ic_step forward RBIC works with coxph", {
  d <- .make_coxph_data()
  fit0 <- survival::coxph(survival::Surv(time, status) ~ 1, data = d)
  scope_upper <- survival::Surv(time, status) ~ age + sex + ph.ecog +
    ph.karno + wt.loss
  P_index <- list(
    demographic = c("age", "sex"),
    clinical    = c("ph.ecog", "ph.karno", "wt.loss")
  )
  result <- ic_step(fit0,
                    scope = list(lower = ~ 1, upper = scope_upper),
                    direction = "forward",
                    criterion = "RBIC",
                    P_index = P_index,
                    trace = 0)
  expect_s3_class(result, "coxph")
  path <- attr(result, "step_path")
  expect_true("RBIC" %in% names(path))
})

test_that("ic_step forward EBIC works with coxph", {
  d <- .make_coxph_data()
  fit0 <- survival::coxph(survival::Surv(time, status) ~ 1, data = d)
  scope_upper <- survival::Surv(time, status) ~ age + sex + ph.ecog +
    ph.karno + wt.loss
  result <- ic_step(fit0,
                    scope = list(lower = ~ 1, upper = scope_upper),
                    direction = "forward",
                    criterion = "EBIC",
                    P = 5,
                    trace = 0)
  expect_s3_class(result, "coxph")
})

test_that("ic_step backward from null coxph returns immediately", {
  d <- .make_coxph_data()
  fit0 <- survival::coxph(survival::Surv(time, status) ~ 1, data = d)
  result <- ic_step(fit0, direction = "backward", criterion = "AIC", trace = 0)
  expect_equal(length(attr(terms(result), "term.labels")), 0L)
})

test_that("ic_step BIC selects sparser model than AIC for coxph", {
  d <- .make_coxph_data()
  fit0 <- survival::coxph(survival::Surv(time, status) ~ 1, data = d)
  scope_upper <- survival::Surv(time, status) ~ age + sex + ph.ecog +
    ph.karno + wt.loss
  scope <- list(lower = ~ 1, upper = scope_upper)
  res_aic <- ic_step(fit0, scope = scope, direction = "forward",
                     criterion = "AIC", trace = 0)
  res_bic <- ic_step(fit0, scope = scope, direction = "forward",
                     criterion = "BIC", trace = 0)
  # BIC should select at most as many terms as AIC
  expect_true(length(attr(terms(res_bic), "term.labels")) <=
              length(attr(terms(res_aic), "term.labels")))
})
