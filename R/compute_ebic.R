#' Compute EBIC for model objects
#'
#' Computes the Extended Bayesian Information Criterion (EBIC) of Chen & Chen
#' (2008).  EBIC adds a sparsity-encouraging penalty to BIC:
#'
#' \deqn{\mathrm{EBIC} = \mathrm{BIC} + 2\,\gamma\,\log\binom{P}{k}}
#'
#' where `P` is the total number of candidate predictors, `k` is the number
#' of selected predictors (nonzero coefficients, excluding the intercept),
#' and `gamma` is a tuning parameter.
#'
#' For `glmnet` and `ncvreg`, EBIC is computed at each point on the
#' regularization path (one value per lambda).
#'
#' @section Inferring P:
#' If `P` is not supplied, it is inferred from the model object:
#' \describe{
#'   \item{`glmnet`}{`nrow(fit$beta)` — the number of columns in the
#'     original design matrix `X`.}
#'   \item{`ncvreg`}{`nrow(fit$beta) - 1` — excluding the intercept row.}
#'   \item{`lm` / `glm`}{`length(coef(fit)) - 1` — the number of terms in
#'     the **fitted** model.  This is correct only when the fitted model
#'     contains all candidate predictors (e.g. the full model); if the model
#'     is a subset, you should supply `P` explicitly.}
#' }
#'
#' @param fit A fitted model object of class `lm`, `glm`, `glmnet`, or
#'   `ncvreg`.
#' @param P Integer. Total number of *candidate* predictors (excluding the
#'   intercept) considered during model selection.  If `NULL` (the default),
#'   inferred from `fit` (see Inferring P section).
#' @param gamma Controls the sparsity tuning parameter.  One of:
#'   \describe{
#'     \item{`"ebic"` (default)}{The data-adaptive rule of Chen & Chen (2008):
#'       \eqn{\gamma = \log(P/\sqrt{n}) / \log P \cdot \mathbf{1}(P \ge \sqrt{n})}.}
#'     \item{A numeric scalar}{Fixed gamma value, e.g. `gamma = 0.5`.}
#'     \item{A function `f(P, k, n)`}{Custom rule; advanced use only.}
#'   }
#' @param ... Currently unused; reserved for future use.
#'
#' @return A named numeric scalar (`lm`/`glm`) or vector (`glmnet`/`ncvreg`)
#'   with the following attributes:
#'   \describe{
#'     \item{`fit_class`}{Character. The primary class of `fit`.}
#'     \item{`k`}{Integer scalar or vector. Number of estimated parameters
#'       (including intercept).}
#'     \item{`gamma`}{Numeric scalar or vector. The computed gamma value(s).}
#'     \item{`P`}{Integer. The value of `P` used.}
#'     \item{`criterion`}{`"EBIC"`.}
#'   }
#'
#' @seealso [compute_bic()], [compute_rbic()], [ic_step()]
#'
#' @references Chen, J., & Chen, Z. (2008). Extended Bayesian information
#'   criteria for model selection with large model spaces. *Biometrika*,
#'   **95**(3), 759–771.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- X[, 1:3] %*% c(2, -1, 0.5) + rnorm(100)
#' fit <- lm(y ~ X)
#'
#' compute_ebic(fit, P = 20)               # explicit P
#' compute_ebic(fit)                        # infers P = 20 from model
#' compute_ebic(fit, P = 20, gamma = 0.5)  # fixed gamma
#' compute_ebic(fit, P = 20, gamma = 0)    # gamma = 0 recovers BIC
#'
#' @export
compute_ebic <- function(fit, P = NULL, gamma = "ebic", ...) UseMethod("compute_ebic")

#' @export
compute_ebic.lm <- function(fit, P = NULL, gamma = "ebic", ...) {
  k_total <- length(stats::coef(fit))
  if (is.null(P)) P <- k_total - 1L
  gammafn <- .resolve_gammafn(gamma)

  n      <- stats::nobs(fit)
  bic    <- stats::BIC(fit)
  k_pred <- k_total - 1L

  gval    <- gammafn(P, k_pred, n)
  penalty <- 2 * gval * .log_choose(P, k_pred)
  val     <- bic + penalty

  .ic_structure(val,
    fit_class = class(fit)[1], k = k_total, criterion = "EBIC",
    extras = list(gamma = gval, P = P)
  )
}

#' @export
compute_ebic.glmnet <- function(fit, P = NULL, gamma = "ebic", ...) {
  if (is.null(P)) P <- nrow(fit$beta)
  .check_glmnet_alpha(fit)
  gammafn <- .resolve_gammafn(gamma)

  n       <- fit$nobs
  bic     <- as.numeric(compute_bic(fit))
  k_total <- .extract_k_glmnet(fit)
  k_pred  <- k_total - 1L

  gvals   <- vapply(k_pred, function(ki) gammafn(P, ki, n), numeric(1))
  penalty <- 2 * gvals * .log_choose(P, k_pred)
  val     <- bic + penalty
  names(val) <- paste0("s", seq_along(fit$lambda) - 1L)

  .ic_structure(val,
    fit_class = "glmnet", k = k_total, criterion = "EBIC",
    extras = list(gamma = gvals, P = P, lambda = fit$lambda)
  )
}

#' @export
compute_ebic.ncvreg <- function(fit, P = NULL, gamma = "ebic", ...) {
  if (is.null(P)) P <- nrow(fit$beta) - 1L
  gammafn <- .resolve_gammafn(gamma)

  n       <- length(fit$y)
  bic     <- as.numeric(stats::BIC(fit))
  k_total <- .extract_k_ncvreg(fit)
  k_pred  <- k_total - 1L

  gvals   <- vapply(k_pred, function(ki) gammafn(P, ki, n), numeric(1))
  penalty <- 2 * gvals * .log_choose(P, k_pred)
  val     <- bic + penalty

  .ic_structure(val,
    fit_class = "ncvreg", k = k_total, criterion = "EBIC",
    extras = list(gamma = gvals, P = P, lambda = fit$lambda)
  )
}
