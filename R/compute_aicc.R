#' Compute AICc for model objects
#'
#' Computes the corrected Akaike Information Criterion (AICc), which adjusts
#' AIC for small sample sizes.  The correction term is `2k(k+1)/(n-k-1)`,
#' where `k` is the number of estimated parameters and `n` is the sample size.
#' AICc converges to AIC as `n -> infinity`.
#'
#' For `glmnet` with `family = "gaussian"`, the log-likelihood is computed
#' exactly (matching `AIC(lm(...))` on the same data).  For other families
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
#'     \item{`criterion`}{`"AICc"`.}
#'   }
#'   For `glmnet` and `ncvreg` results, a `lambda` attribute holds the
#'   corresponding lambda values.  Values for which `n <= k + 1` (undefined
#'   correction) are returned as `NA`.
#'
#' @seealso [compute_aic()], [compute_bic()], [compute_ebic()],
#'   [compute_rbic()], [ic_step()]
#'
#' @references Hurvich, C. M., & Tsai, C.-L. (1989). Regression and time
#'   series model selection in small samples. *Biometrika*, **76**(2), 297–307.
#'
#' @examples
#' fit <- lm(mpg ~ wt + hp, data = mtcars)
#' compute_aicc(fit)
#'
#' fit_glm <- glm(am ~ wt + hp, data = mtcars, family = binomial)
#' compute_aicc(fit_glm)
#'
#' @export
compute_aicc <- function(fit, ...) UseMethod("compute_aicc")

#' @export
compute_aicc.lm <- function(fit, ...) {
  n   <- stats::nobs(fit)
  k   <- length(stats::coef(fit))
  val <- stats::AIC(fit) + .aicc_correction(k, n)
  .ic_structure(val, fit_class = class(fit)[1], k = k, criterion = "AICc")
}

#' @export
compute_aicc.glmnet <- function(fit, ...) {
  .check_glmnet_alpha(fit)
  n      <- fit$nobs
  k      <- .extract_k_glmnet(fit)
  loglik <- .glmnet_loglik(fit)
  aic    <- -2 * loglik + 2 * k
  val    <- aic + .aicc_correction(k, n)
  names(val) <- paste0("s", seq_along(fit$lambda) - 1L)
  .ic_structure(val,
    fit_class = "glmnet", k = k, criterion = "AICc",
    extras = list(lambda = fit$lambda)
  )
}

#' @export
compute_aicc.ncvreg <- function(fit, ...) {
  n   <- length(fit$y)
  k   <- .extract_k_ncvreg(fit)
  val <- stats::AIC(fit) + .aicc_correction(k, n)
  .ic_structure(val,
    fit_class = "ncvreg", k = k, criterion = "AICc",
    extras = list(lambda = fit$lambda)
  )
}
