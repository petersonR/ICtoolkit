#' Compute HQIC for model objects
#'
#' Computes the Hannan-Quinn Information Criterion (HQIC).  HQIC uses a penalty
#' term that lies between AIC and BIC:
#'
#' \deqn{\mathrm{HQIC} = -2\,\log L + 2\,k\,\log(\log(n))}
#'
#' where `k` is the number of estimated parameters and `n` is the sample size.
#' HQIC is strongly consistent (like BIC) yet less conservative than BIC for
#' moderate sample sizes.
#'
#' For `lm` and `glm` models the result is a scalar; for `glmnet` and `ncvreg`
#' objects the result is a numeric vector of length `length(lambda)`, one value
#' per point on the regularization path.
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
#'     \item{`criterion`}{`"HQIC"`.}
#'   }
#'   For `glmnet` and `ncvreg` results, a `lambda` attribute holds the
#'   corresponding lambda values.
#'
#' @seealso [compute_aic()], [compute_bic()], [compute_ebic()], [ic_step()]
#'
#' @references Hannan, E. J., & Quinn, B. G. (1979). The determination of the
#'   order of an autoregression. *Journal of the Royal Statistical Society:
#'   Series B*, **41**(2), 190--195.
#'
#' @examples
#' fit <- lm(mpg ~ wt + hp, data = mtcars)
#' compute_hqic(fit)
#'
#' fit_glm <- glm(am ~ wt + hp, data = mtcars, family = binomial)
#' compute_hqic(fit_glm)
#'
#' @export
compute_hqic <- function(fit, ...) UseMethod("compute_hqic")

#' @export
compute_hqic.lm <- function(fit, ...) {
  n      <- stats::nobs(fit)
  k      <- length(stats::coef(fit))
  ll     <- stats::logLik(fit)
  df     <- attr(ll, "df")   # includes sigma for lm; matches stats::AIC/BIC convention
  val    <- -2 * as.numeric(ll) + 2 * df * log(log(n))
  .ic_structure(val, fit_class = class(fit)[1], k = k, criterion = "HQIC")
}

#' @export
compute_hqic.glmnet <- function(fit, ...) {
  .check_glmnet_alpha(fit)
  n      <- fit$nobs
  k      <- .extract_k_glmnet(fit)
  loglik <- .glmnet_loglik(fit)
  val    <- -2 * loglik + 2 * k * log(log(n))
  names(val) <- paste0("s", seq_along(fit$lambda) - 1L)
  .ic_structure(val,
    fit_class = "glmnet", k = k, criterion = "HQIC",
    extras = list(lambda = fit$lambda)
  )
}

#' @export
compute_hqic.ncvreg <- function(fit, ...) {
  n      <- length(fit$y)
  k      <- .extract_k_ncvreg(fit)
  loglik <- as.numeric(stats::logLik(fit))
  val    <- -2 * loglik + 2 * k * log(log(n))
  .ic_structure(val,
    fit_class = "ncvreg", k = k, criterion = "HQIC",
    extras = list(lambda = fit$lambda)
  )
}

#' @export
compute_hqic.coxph <- function(fit, ...) {
  n   <- fit$n
  k   <- .extract_k_coxph(fit)
  ll  <- stats::logLik(fit)
  df  <- attr(ll, "df")
  val <- -2 * as.numeric(ll) + 2 * df * log(log(n))
  .ic_structure(val, fit_class = "coxph", k = k, criterion = "HQIC")
}

#' @export
compute_hqic.ncvsurv <- function(fit, ...) {
  .check_ncvreg_ns()
  n      <- .extract_n_ncvsurv(fit)
  k      <- .extract_k_ncvsurv(fit)
  loglik <- as.numeric(stats::logLik(fit))
  val    <- -2 * loglik + 2 * k * log(log(n))
  .ic_structure(val,
    fit_class = "ncvsurv", k = k, criterion = "HQIC",
    extras = list(lambda = fit$lambda)
  )
}
