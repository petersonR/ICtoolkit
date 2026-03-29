#' Compute mBIC for model objects
#'
#' Computes the modified Bayesian Information Criterion (mBIC) of Bogdan,
#' Ghosh, & Doerge (2004).  mBIC adds a complexity penalty beyond BIC that
#' accounts for the total number of candidate predictors `P` and a prior
#' expected model size `kappa`:
#'
#' \deqn{\mathrm{mBIC} = \mathrm{BIC} + 2\,k_{\text{pred}}\,
#'   \log\!\bigl(P/\kappa - 1\bigr)}
#'
#' where `k_pred` is the number of selected predictors (nonzero coefficients,
#' excluding the intercept).  When `kappa` is small relative to `P`, the
#' penalty encourages sparser models than BIC.
#'
#' For `glmnet` and `ncvreg`, mBIC is computed at each point on the
#' regularization path (one value per lambda).
#'
#' @section Inferring P:
#' If `P` is not supplied, it is inferred identically to [compute_ebic()]:
#' \describe{
#'   \item{`glmnet`}{`nrow(fit$beta)`.}
#'   \item{`ncvreg`}{`nrow(fit$beta) - 1`.}
#'   \item{`lm` / `glm`}{`length(coef(fit)) - 1`.}
#' }
#'
#' @param fit A fitted model object of class `lm`, `glm`, `glmnet`, or
#'   `ncvreg`.
#' @param P Integer. Total number of *candidate* predictors (excluding the
#'   intercept).  If `NULL` (the default), inferred from `fit`.
#' @param kappa Positive numeric scalar controlling the prior expected number
#'   of true signals.  Must satisfy `P > kappa`.  Default is `4`.
#' @param ... Currently unused; reserved for future use.
#'
#' @return A named numeric scalar (`lm`/`glm`) or vector (`glmnet`/`ncvreg`)
#'   with the following attributes:
#'   \describe{
#'     \item{`fit_class`}{Character. The primary class of `fit`.}
#'     \item{`k`}{Integer scalar or vector. Number of estimated parameters
#'       (including intercept).}
#'     \item{`P`}{Integer. The value of `P` used.}
#'     \item{`kappa`}{Numeric. The value of `kappa` used.}
#'     \item{`criterion`}{`"mBIC"`.}
#'   }
#'
#' @seealso [compute_mbic2()], [compute_ebic()], [compute_bic()], [ic_step()]
#'
#' @references Bogdan, M., Ghosh, J. K., & Doerge, R. W. (2004). Modifying
#'   the Schwarz Bayesian information criterion to locate multiple interacting
#'   quantitative trait loci. *Genetics*, **167**(2), 989--999.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- X[, 1:3] %*% c(2, -1, 0.5) + rnorm(100)
#' fit <- lm(y ~ X)
#'
#' compute_mbic(fit, P = 20)                  # default kappa = 4
#' compute_mbic(fit, P = 20, kappa = 2)       # expect fewer signals
#' compute_mbic(fit, P = 20, kappa = 10)      # expect more signals
#'
#' @export
compute_mbic <- function(fit, P = NULL, kappa = 4, ...) UseMethod("compute_mbic")

#' @export
compute_mbic.lm <- function(fit, P = NULL, kappa = 4, ...) {
  k_total <- length(stats::coef(fit))
  if (is.null(P)) P <- k_total - 1L
  .check_kappa(P, kappa)

  bic    <- stats::BIC(fit)
  k_pred <- k_total - 1L
  val    <- bic + .mbic_penalty(k_pred, P, kappa)

  .ic_structure(val,
    fit_class = class(fit)[1], k = k_total, criterion = "mBIC",
    extras = list(P = P, kappa = kappa)
  )
}

#' @export
compute_mbic.glmnet <- function(fit, P = NULL, kappa = 4, ...) {
  if (is.null(P)) P <- nrow(fit$beta)
  .check_kappa(P, kappa)
  .check_glmnet_alpha(fit)

  bic     <- as.numeric(compute_bic(fit))
  k_total <- .extract_k_glmnet(fit)
  k_pred  <- k_total - 1L
  penalty <- .mbic_penalty(k_pred, P, kappa)
  val     <- bic + penalty
  names(val) <- paste0("s", seq_along(fit$lambda) - 1L)

  .ic_structure(val,
    fit_class = "glmnet", k = k_total, criterion = "mBIC",
    extras = list(P = P, kappa = kappa, lambda = fit$lambda)
  )
}

#' @export
compute_mbic.ncvreg <- function(fit, P = NULL, kappa = 4, ...) {
  if (is.null(P)) P <- nrow(fit$beta) - 1L
  .check_kappa(P, kappa)

  bic     <- as.numeric(stats::BIC(fit))
  k_total <- .extract_k_ncvreg(fit)
  k_pred  <- k_total - 1L
  penalty <- .mbic_penalty(k_pred, P, kappa)
  val     <- bic + penalty

  .ic_structure(val,
    fit_class = "ncvreg", k = k_total, criterion = "mBIC",
    extras = list(P = P, kappa = kappa, lambda = fit$lambda)
  )
}
