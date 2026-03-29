#' Compute BIC for model objects
#'
#' Computes the Bayesian Information Criterion (BIC), also known as the
#' Schwarz criterion.  For `lm` and `glm` models the result is a scalar;
#' for `glmnet` and `ncvreg` objects the result is a numeric vector of
#' length `length(lambda)`, one value per point on the regularization path.
#'
#' For `glmnet` with `family = "gaussian"`, the log-likelihood is computed
#' exactly (matching `BIC(lm(...))` on the same data).  For other families
#' the deviance-based approximation is used (see [compute_aic()] for details).
#'
#' @param fit A fitted model object of class `lm`, `glm`, `glmnet`, or
#'   `ncvreg`.
#' @param ... Currently unused; reserved for future use.
#'
#' @return A named numeric scalar (`lm`/`glm`) or vector (`glmnet`/`ncvreg`)
#'   with the following attributes:
#'   \describe{
#'     \item{`fit_class`}{Character. The primary class of `fit`.}
#'     \item{`k`}{Integer scalar or vector. Number of estimated parameters.}
#'     \item{`criterion`}{`"BIC"`.}
#'   }
#'   For `glmnet` and `ncvreg` results, a `lambda` attribute holds the
#'   corresponding lambda values.
#'
#' @seealso [compute_aic()], [compute_aicc()], [compute_ebic()],
#'   [compute_rbic()], [ic_step()]
#'
#' @references Schwarz, G. (1978). Estimating the dimension of a model.
#'   *The Annals of Statistics*, **6**(2), 461–464.
#'
#' @examples
#' fit <- lm(mpg ~ wt + hp, data = mtcars)
#' compute_bic(fit)
#'
#' fit_glm <- glm(am ~ wt + hp, data = mtcars, family = binomial)
#' compute_bic(fit_glm)
#'
#' @export
compute_bic <- function(fit, ...) UseMethod("compute_bic")

#' @export
compute_bic.lm <- function(fit, ...) {
  val <- stats::BIC(fit)
  k   <- length(stats::coef(fit))
  .ic_structure(val, fit_class = class(fit)[1], k = k, criterion = "BIC")
}

#' @export
compute_bic.glmnet <- function(fit, ...) {
  .check_glmnet_alpha(fit)
  n      <- fit$nobs
  k      <- .extract_k_glmnet(fit)
  loglik <- .glmnet_loglik(fit)
  val    <- -2 * loglik + log(n) * k
  names(val) <- paste0("s", seq_along(fit$lambda) - 1L)
  .ic_structure(val,
    fit_class = "glmnet", k = k, criterion = "BIC",
    extras = list(lambda = fit$lambda)
  )
}

#' @export
compute_bic.ncvreg <- function(fit, ...) {
  val <- stats::BIC(fit)   # relies on logLik.ncvreg
  k   <- .extract_k_ncvreg(fit)
  .ic_structure(val,
    fit_class = "ncvreg", k = k, criterion = "BIC",
    extras = list(lambda = fit$lambda)
  )
}

#' @export
compute_bic.coxph <- function(fit, ...) {
  val <- stats::BIC(fit)
  k   <- .extract_k_coxph(fit)
  .ic_structure(val, fit_class = "coxph", k = k, criterion = "BIC")
}

#' @export
compute_bic.ncvsurv <- function(fit, ...) {
  val <- stats::BIC(fit)   # relies on logLik.ncvsurv
  k   <- .extract_k_ncvsurv(fit)
  .ic_structure(val,
    fit_class = "ncvsurv", k = k, criterion = "BIC",
    extras = list(lambda = fit$lambda)
  )
}
