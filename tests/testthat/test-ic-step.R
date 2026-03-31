# Shared helper: covid_eol augmented with noise predictors and I(age^2).
# The quadratic age term creates path-dependent selection (backward retains
# age + I(age^2) while forward does not), and noise terms cause AIC/BIC
# divergence (AIC retains some noise; BIC does not).
.make_step_data <- function() {
  set.seed(42)
  d <- covid_eol
  for (i in 1:8) d[[paste0("noise", i)]] <- rnorm(nrow(d))
  d
}

# Full-model formulas are defined *inside* each test so that MASS::stepAIC
# (which evaluates the formula in its environment) can locate the data.

# ---------------------------------------------------------------------------
# Basic interface tests
# ---------------------------------------------------------------------------

test_that("ic_step returns the same class as the input", {
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  fit_full <- lm(ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8, data = d)
  result <- ic_step(fit_null,
                    scope     = list(lower = fit_null, upper = fit_full),
                    direction = "forward",
                    criterion = "AIC",
                    trace     = 0)
  expect_s3_class(result, "lm")
})

test_that("ic_step attaches step_path attribute", {
  d        <- .make_step_data()
  fit_full <- lm(ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8, data = d)
  result <- ic_step(fit_full, direction = "backward", criterion = "BIC", trace = 0)
  path   <- attr(result, "step_path")
  expect_s3_class(path, "data.frame")
  expect_true("step"   %in% names(path))
  expect_true("action" %in% names(path))
  expect_true("BIC"    %in% names(path))
})

test_that("ic_step errors without P for EBIC", {
  fit_null <- lm(ies_total ~ 1, data = covid_eol)
  fit_full <- lm(ies_total ~ age + sex, data = covid_eol)
  expect_error(
    ic_step(fit_null, scope = list(lower = fit_null, upper = fit_full),
            criterion = "EBIC", trace = 0),
    "'P'"
  )
})

test_that("ic_step errors without P_index for RBIC", {
  fit_null <- lm(ies_total ~ 1, data = covid_eol)
  fit_full <- lm(ies_total ~ age + sex, data = covid_eol)
  expect_error(
    ic_step(fit_null, scope = list(lower = fit_null, upper = fit_full),
            criterion = "RBIC", trace = 0),
    "'P_index'"
  )
})

test_that("ic_step errors without P for mBIC", {
  fit_null <- lm(ies_total ~ 1, data = covid_eol)
  fit_full <- lm(ies_total ~ age + sex, data = covid_eol)
  expect_error(
    ic_step(fit_null, scope = list(lower = fit_null, upper = fit_full),
            criterion = "mBIC", trace = 0),
    "'P'"
  )
})

test_that("ic_step errors without P for mBIC2", {
  fit_null <- lm(ies_total ~ 1, data = covid_eol)
  fit_full <- lm(ies_total ~ age + sex, data = covid_eol)
  expect_error(
    ic_step(fit_null, scope = list(lower = fit_null, upper = fit_full),
            criterion = "mBIC2", trace = 0),
    "'P'"
  )
})

test_that("ic_step forward HQIC runs and returns lm", {
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  result <- ic_step(fit_null, scope = scope, direction = "forward",
                    criterion = "HQIC", trace = 0)
  expect_s3_class(result, "lm")
  path <- attr(result, "step_path")
  expect_true("HQIC" %in% names(path))
})

test_that("ic_step forward mBIC runs and returns lm", {
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  result <- ic_step(fit_null, scope = scope, direction = "forward",
                    criterion = "mBIC", P = 13, kappa = 4, trace = 0)
  expect_s3_class(result, "lm")
  path <- attr(result, "step_path")
  expect_true("mBIC" %in% names(path))
})

test_that("ic_step forward mBIC2 runs and returns lm", {
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  result <- ic_step(fit_null, scope = scope, direction = "forward",
                    criterion = "mBIC2", P = 13, kappa = 4, trace = 0)
  expect_s3_class(result, "lm")
  path <- attr(result, "step_path")
  expect_true("mBIC2" %in% names(path))
})

test_that("ic_step backward from null model returns null model immediately", {
  fit_null <- lm(ies_total ~ 1, data = covid_eol)
  result <- ic_step(fit_null, direction = "backward", criterion = "AIC", trace = 0)
  expect_equal(length(attr(terms(result), "term.labels")), 0L)
})

test_that("ic_step respects lower bound", {
  d         <- .make_step_data()
  fit_start <- lm(ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8, data = d)
  fit_lower <- lm(ies_total ~ age, data = d)
  result <- ic_step(fit_start,
                    scope     = list(lower = fit_lower, upper = fit_start),
                    direction = "backward",
                    criterion = "BIC",
                    trace     = 0)
  expect_true("age" %in% attr(terms(result), "term.labels"))
})

# ---------------------------------------------------------------------------
# Agreement with MASS::stepAIC  -- lm (continuous IES score)
# ---------------------------------------------------------------------------

test_that("ic_step backward AIC matches MASS::stepAIC backward on lm", {
  d        <- .make_step_data()
  fit_full <- lm(ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8, data = d)
  res_ours <- ic_step(fit_full, direction = "backward", criterion = "AIC", trace = 0)
  res_mass <- MASS::stepAIC(fit_full, direction = "backward", trace = 0)
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

test_that("ic_step forward AIC matches MASS::stepAIC forward on lm", {
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  res_ours <- ic_step(fit_null, scope = scope, direction = "forward",
                      criterion = "AIC", trace = 0)
  res_mass <- MASS::stepAIC(fit_null, scope = scope, direction = "forward",
                            trace = 0)
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

test_that("ic_step both AIC matches MASS::stepAIC both on lm", {
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  res_ours <- ic_step(fit_null, scope = scope, direction = "both",
                      criterion = "AIC", trace = 0)
  res_mass <- MASS::stepAIC(fit_null, scope = scope, direction = "both",
                            trace = 0)
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

# ---------------------------------------------------------------------------
# Agreement with MASS::stepAIC  -- glm (binary PTSD outcome)
# ---------------------------------------------------------------------------

test_that("ic_step backward AIC matches MASS::stepAIC backward on glm", {
  d        <- .make_step_data()
  fit_full <- glm(death ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8,
    data = d, family = binomial)
  suppressWarnings({
    res_ours <- ic_step(fit_full, direction = "backward", criterion = "AIC", trace = 0)
    res_mass <- MASS::stepAIC(fit_full, direction = "backward", trace = 0)
  })
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

test_that("ic_step forward AIC matches MASS::stepAIC forward on glm", {
  d        <- .make_step_data()
  fit_null <- glm(death ~ 1, data = d, family = binomial)
  scope    <- list(lower = ~ 1, upper = death ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  suppressWarnings({
    res_ours <- ic_step(fit_null, scope = scope, direction = "forward",
                        criterion = "AIC", trace = 0)
    res_mass <- MASS::stepAIC(fit_null, scope = scope, direction = "forward",
                              trace = 0)
  })
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

test_that("ic_step both AIC matches MASS::stepAIC both on glm", {
  d        <- .make_step_data()
  fit_null <- glm(death ~ 1, data = d, family = binomial)
  scope    <- list(lower = ~ 1, upper = death ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  suppressWarnings({
    res_ours <- ic_step(fit_null, scope = scope, direction = "both",
                        criterion = "AIC", trace = 0)
    res_mass <- MASS::stepAIC(fit_null, scope = scope, direction = "both",
                              trace = 0)
  })
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

# ---------------------------------------------------------------------------
# Agreement with MASS::stepAIC using k=log(n) for BIC  -- lm
# ---------------------------------------------------------------------------

test_that("ic_step backward BIC matches MASS::stepAIC with k=log(n) on lm", {
  d        <- .make_step_data()
  fit_full <- lm(ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8, data = d)
  n        <- nobs(fit_full)
  res_ours <- ic_step(fit_full, direction = "backward", criterion = "BIC", trace = 0)
  res_mass <- MASS::stepAIC(fit_full, direction = "backward", k = log(n), trace = 0)
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

test_that("ic_step forward BIC matches MASS::stepAIC with k=log(n) on lm", {
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  n        <- nrow(d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  res_ours <- ic_step(fit_null, scope = scope, direction = "forward",
                      criterion = "BIC", trace = 0)
  res_mass <- MASS::stepAIC(fit_null, scope = scope, direction = "forward",
                            k = log(n), trace = 0)
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

test_that("ic_step both BIC matches MASS::stepAIC with k=log(n) on lm", {
  d        <- .make_step_data()
  fit_full <- lm(ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8, data = d)
  n        <- nobs(fit_full)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8)
  res_ours <- ic_step(fit_full, scope = scope, direction = "both",
                      criterion = "BIC", trace = 0)
  res_mass <- MASS::stepAIC(fit_full, scope = scope, direction = "both",
                            k = log(n), trace = 0)
  expect_equal(sort(attr(terms(res_ours), "term.labels")),
               sort(attr(terms(res_mass), "term.labels")))
})

# ---------------------------------------------------------------------------
# Parallel evaluation via cl argument
# ---------------------------------------------------------------------------

test_that("ic_step errors on invalid cl argument", {
  fit_null <- lm(ies_total ~ 1, data = covid_eol)
  fit_full <- lm(ies_total ~ age + sex, data = covid_eol)
  expect_error(
    ic_step(fit_null, scope = list(lower = fit_null, upper = fit_full),
            criterion = "AIC", cl = "bad", trace = 0),
    "'cl' must be an integer"
  )
})

test_that("ic_step with cl=integer gives same result as sequential (forward BIC)", {
  skip_on_cran()
  skip_if_not_installed("pbapply")
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race +
    ethnicity + sex + noise1 + noise2 + noise3 + noise4 + noise5 + noise6 +
    noise7 + noise8)
  res_seq <- ic_step(fit_null, scope = scope, direction = "forward",
                     criterion = "BIC", trace = 0)
  res_par <- ic_step(fit_null, scope = scope, direction = "forward",
                     criterion = "BIC", cl = 2L, trace = 0)
  expect_equal(sort(attr(terms(res_par), "term.labels")),
               sort(attr(terms(res_seq), "term.labels")))
})

test_that("ic_step with cl=integer gives same result as sequential (backward AIC)", {
  skip_on_cran()
  skip_if_not_installed("pbapply")
  d        <- .make_step_data()
  fit_full <- lm(ies_total ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8,
    data = d)
  res_seq <- ic_step(fit_full, direction = "backward", criterion = "AIC",
                     trace = 0)
  res_par <- ic_step(fit_full, direction = "backward", criterion = "AIC",
                     cl = 2L, trace = 0)
  expect_equal(sort(attr(terms(res_par), "term.labels")),
               sort(attr(terms(res_seq), "term.labels")))
})

test_that("ic_step with pre-made cluster object works", {
  skip_on_cran()
  skip_if_not_installed("pbapply")
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race +
    ethnicity + sex + noise1 + noise2 + noise3 + noise4 + noise5 + noise6 +
    noise7 + noise8)
  cl <- parallel::makeCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  res_seq <- ic_step(fit_null, scope = scope, direction = "forward",
                     criterion = "BIC", trace = 0)
  res_par <- ic_step(fit_null, scope = scope, direction = "forward",
                     criterion = "BIC", cl = cl, trace = 0)
  expect_equal(sort(attr(terms(res_par), "term.labels")),
               sort(attr(terms(res_seq), "term.labels")))
})

test_that("ic_step with cl for glm gives same result as sequential", {
  skip_on_cran()
  skip_if_not_installed("pbapply")
  d        <- .make_step_data()
  fit_full <- glm(death ~ age + I(age^2) + race + ethnicity + sex +
    noise1 + noise2 + noise3 + noise4 + noise5 + noise6 + noise7 + noise8,
    data = d, family = binomial)
  suppressWarnings({
    res_seq <- ic_step(fit_full, direction = "backward", criterion = "AIC",
                       trace = 0)
    res_par <- ic_step(fit_full, direction = "backward", criterion = "AIC",
                       cl = 2L, trace = 0)
  })
  expect_equal(sort(attr(terms(res_par), "term.labels")),
               sort(attr(terms(res_seq), "term.labels")))
})

test_that("ic_step parallel step_path matches sequential step_path", {
  skip_on_cran()
  skip_if_not_installed("pbapply")
  d        <- .make_step_data()
  fit_null <- lm(ies_total ~ 1, data = d)
  scope    <- list(lower = ~ 1, upper = ies_total ~ age + I(age^2) + race +
    ethnicity + sex + noise1 + noise2 + noise3 + noise4 + noise5 + noise6 +
    noise7 + noise8)
  res_seq <- ic_step(fit_null, scope = scope, direction = "forward",
                     criterion = "BIC", trace = 0)
  res_par <- ic_step(fit_null, scope = scope, direction = "forward",
                     criterion = "BIC", cl = 2L, trace = 0)
  expect_equal(attr(res_par, "step_path"), attr(res_seq, "step_path"))
})

test_that("ic_step parallel is faster than sequential on a wide dataset", {
  skip_on_cran()
  skip_on_os("windows")
  set.seed(1)
  n <- 5000
  p <- 300
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  d <- as.data.frame(X)
  d$y <- rbinom(n, 1, 0.5)

  upper_terms <- paste0("x", seq_len(p))
  fit_null <- glm(y ~ 1, data = d, family = binomial)
  scope <- list(lower = ~ 1, upper = upper_terms)

  # R CMD check enforces a 2-core limit via _R_CHECK_LIMIT_CORES_
  n_cores <- 2L

  t_seq <- system.time(
    suppressWarnings(ic_step(fit_null, scope = scope, direction = "forward",
            criterion = "BIC", trace = 0, steps = 3L))
  )[["elapsed"]]

  t_par <- system.time(
    suppressWarnings(ic_step(fit_null, scope = scope, direction = "forward",
            criterion = "BIC", cl = n_cores, trace = 0, steps = 3L))
  )[["elapsed"]]

  expect_lt(t_par, t_seq)
})

# ---------------------------------------------------------------------------
# Cox PH stepwise + warm starts
# ---------------------------------------------------------------------------

test_that("ic_step forward BIC works for coxph", {
  skip_if_not_installed("survival")
  lung <- survival::lung
  lung <- lung[complete.cases(lung[, c("time", "status", "age", "sex",
                                        "ph.ecog", "wt.loss")]), ]
  fit0 <- survival::coxph(survival::Surv(time, status) ~ 1, data = lung)
  scope <- list(lower = ~ 1,
                upper = survival::Surv(time, status) ~ age + sex + ph.ecog + wt.loss)
  result <- ic_step(fit0, scope = scope, direction = "forward",
                    criterion = "BIC", trace = 0)
  expect_s3_class(result, "coxph")
  path <- attr(result, "step_path")
  expect_s3_class(path, "data.frame")
  expect_true(nrow(path) > 1L)
})

