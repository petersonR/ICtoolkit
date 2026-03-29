# Internal utility functions for ICtoolkit.
# None of these are exported.

# ---------------------------------------------------------------------------
# glmnet helpers
# ---------------------------------------------------------------------------

# Number of nonzero coefficients (+1 for intercept) at each lambda.
# Returns integer vector of length nlambda.
.extract_k_glmnet <- function(fit) {
  # fit$df counts nonzero penalized coefficients only (excludes intercept).
  intercept <- !isFALSE(fit$call$intercept)
  fit$df + as.integer(intercept)
}

# Deviance at each lambda: nulldev * (1 - dev.ratio).
# For Gaussian this equals RSS; for other families it equals -2*logLik
# relative to the saturated model.
.extract_deviance_glmnet <- function(fit) {
  fit$nulldev * (1 - fit$dev.ratio)
}

# Proper log-likelihood for glmnet at each lambda.
# Gaussian: exact (matches logLik.lm).
# Other families: -deviance/2 (excludes the saturated-model constant;
#   correct for binomial since sat. logLik = 0, approximate for Poisson).
.glmnet_loglik <- function(fit) {
  n   <- fit$nobs
  dev <- .extract_deviance_glmnet(fit)

  if (.glmnet_family(fit) == "gaussian") {
    -n / 2 * (log(2 * pi) + 1 + log(dev / n))
  } else {
    -dev / 2
  }
}

# Detect the glmnet family as a character string.
.glmnet_family <- function(fit) {
  fam <- fit$call$family
  if (is.null(fam)) return("gaussian")
  fam_str <- as.character(fam)
  fam_str[1]   # handles calls like gaussian(), binomial()
}

# Warn when alpha < 1 (L2 shrinkage present).
.check_glmnet_alpha <- function(fit) {
  alpha <- fit$call$alpha
  if (is.null(alpha)) alpha <- 1   # glmnet default is lasso
  alpha <- suppressWarnings(as.numeric(alpha))
  if (!is.na(alpha) && alpha < 1) {
    warning(
      "glmnet model has alpha = ", alpha, " (< 1), introducing L2 shrinkage.\n",
      "Using the count of nonzero coefficients as k may be inappropriate ",
      "for information criteria under elastic-net or ridge penalties.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# ncvreg helpers
# ---------------------------------------------------------------------------

# Number of nonzero coefficients (including intercept) at each lambda.
# ncvreg$beta is (p+1) x nlambda with the intercept in the first row.
.extract_k_ncvreg <- function(fit) {
  beta <- fit$beta
  # Row 1 is the intercept; count nonzeros in rows 2:p+1, then add 1 for intercept.
  colSums(beta[-1, , drop = FALSE] != 0) + 1L
}

# ---------------------------------------------------------------------------
# ncvsurv helpers
# ---------------------------------------------------------------------------

# Sample size for ncvsurv objects.
# fit$y is a Surv object (matrix-like), so length(fit$y) gives 2*n.
.extract_n_ncvsurv <- function(fit) {
  as.integer(fit$n)
}

# Number of nonzero coefficients at each lambda.
# ncvsurv$beta is p x nlambda with NO intercept row (Cox models have no
# intercept).
.extract_k_ncvsurv <- function(fit) {
  colSums(fit$beta != 0)
}

# ---------------------------------------------------------------------------
# coxph helpers
# ---------------------------------------------------------------------------

# Number of estimated parameters for a coxph object (no intercept).
.extract_k_coxph <- function(fit) {
  length(stats::coef(fit))
}

# ---------------------------------------------------------------------------
# Shared penalty helpers
# ---------------------------------------------------------------------------

# Default gamma function for EBIC and RBIC (Chen & Chen 2008).
#
# For EBIC: P is a scalar (total candidate predictors).
# For RBIC: P is a VECTOR of group sizes [p_1, p_2, ..., p_K].
#   The default uses sum(P) -- the total number of candidates across all groups --
#   to compute a SINGLE gamma that is then applied uniformly across groups.
#   This matches the paper's default formula:
#     RBIC = BIC + 2 * gamma_EBIC * sum_k{ log C(p_k, m_k) }
#   where gamma_EBIC is based on P_total = sum(p_k).
#
# A user-supplied gammafn may return either a scalar (applied uniformly to all
# groups) or a vector of length K (group-specific gammas).
#
# k: vector of selected counts (unused in the default; available for custom rules).
# n: sample size.
.default_gammafn <- function(P, k, n) {
  P_total <- sum(P)
  if (P_total <= 1 || P_total < sqrt(n)) return(0)
  log(P_total / sqrt(n)) / log(P_total)
}

# log(C(n, k)) with graceful handling of boundary cases.
.log_choose <- function(n, k) {
  ifelse(k <= 0 | k >= n, 0, lchoose(n, k))
}

# ---------------------------------------------------------------------------
# mBIC / mBIC2 helpers
# ---------------------------------------------------------------------------

# Validate that P > kappa (required for log(P/kappa - 1) to be defined).
.check_kappa <- function(P, kappa) {
  if (!is.numeric(kappa) || length(kappa) != 1 || kappa <= 0)
    stop("'kappa' must be a positive numeric scalar.")
  if (P <= kappa)
    stop("'P' must be greater than 'kappa' (got P = ", P, ", kappa = ", kappa, ").")
}

# mBIC extra penalty beyond BIC: 2 * k_pred * log(P/kappa - 1).
# Vectorised over k_pred.
.mbic_penalty <- function(k_pred, P, kappa) {
  2 * k_pred * log(P / kappa - 1)
}

# mBIC2 extra penalty beyond BIC: mBIC penalty minus 2*log(C(P, k_pred)).
# Vectorised over k_pred.
.mbic2_penalty <- function(k_pred, P, kappa) {
  .mbic_penalty(k_pred, P, kappa) - 2 * .log_choose(P, k_pred)
}

# ---------------------------------------------------------------------------
# AICc helper
# ---------------------------------------------------------------------------

# AICc small-sample correction term: 2k(k+1)/(n-k-1).
# Vectorised over k. Returns NA where n <= k+1.
.aicc_correction <- function(k, n) {
  denom <- n - k - 1
  denom[denom <= 0] <- NA_real_
  2 * k * (k + 1) / denom
}

# ---------------------------------------------------------------------------
# Gamma resolution
# ---------------------------------------------------------------------------

# Resolve the user-facing `gamma` argument into an internal gammafn.
# `gamma` may be:
#   "ebic"       -> default Chen & Chen rule based on sum(P)
#   "per_group"  -> Chen & Chen rule applied to each group's P_g separately
#   numeric scalar/vector -> fixed gamma value(s), returned as-is
#   function     -> passed through unchanged (advanced use)
#
# The returned function always has signature f(P, k, n) where P and k may be
# scalar (EBIC) or vectors (RBIC), and returns a scalar or vector of gammas.
.resolve_gammafn <- function(gamma) {
  if (is.function(gamma)) return(gamma)

  if (is.character(gamma)) {
    gamma <- match.arg(gamma, c("ebic", "per_group"))
    if (gamma == "ebic") return(.default_gammafn)
    # per_group: apply the Chen & Chen formula to each group's P_g individually
    return(function(P, k, n) {
      sapply(P, function(p) {
        if (p <= 1 || p < sqrt(n)) return(0)
        log(p / sqrt(n)) / log(p)
      })
    })
  }

  if (is.numeric(gamma)) {
    gval <- gamma   # capture fixed value(s)
    return(function(P, k, n) gval)
  }

  stop("'gamma' must be \"ebic\", \"per_group\", a numeric value, or a function.")
}

# ---------------------------------------------------------------------------
# Attribute attachment
# ---------------------------------------------------------------------------

# Attach the standard ICtoolkit attributes to a numeric result.
# `extras` is a named list of additional attributes (e.g. gamma, P_index).
.ic_structure <- function(val, fit_class, k, criterion, extras = list()) {
  attr(val, "fit_class") <- fit_class
  attr(val, "k")         <- k
  attr(val, "criterion") <- criterion
  for (nm in names(extras)) attr(val, nm) <- extras[[nm]]
  val
}
