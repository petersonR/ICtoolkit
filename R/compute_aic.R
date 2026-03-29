#' Compute AIC for model objects
#'
#' Computes the Akaike Information Criterion (AIC) for fitted model objects.
#' For `lm` and `glm` models the result is a scalar; for `glmnet` and `ncvreg`
#' objects the result is a numeric vector of length `length(lambda)`,
#' one value per point on the regularization path.
#'
#' For `glmnet` with `family = "gaussian"`, the log-likelihood is computed
#' exactly (matching `AIC(lm(...))` on the same data).  For other families
#' the deviance-based approximation `-deviance/2` is used; this is exact for
#' `binomial` and approximate for `poisson` and `multinomial`.
#'
#' @param fit A fitted model object of class `lm`, `glm`, `glmnet`, or
#'   `ncvreg`.
#' @param ... Currently unused; reserved for future use.
#'
#' @return A named numeric scalar (`lm`/`glm`) or vector (`glmnet`/`ncvreg`)
#'   with the following attributes:
#'   \describe{
#'     \item{`fit_class`}{Character. The primary class of `fit`.}
#'     \item{`k`}{Integer scalar or vector. Number of estimated parameters
#'       (including the intercept).}
#'     \item{`criterion`}{`"AIC"`.}
#'   }
#'   For `glmnet` and `ncvreg` results, the vector is named by lambda index
#'   and an additional `lambda` attribute holds the corresponding lambda values.
#'
#' @seealso [compute_aicc()], [compute_bic()], [compute_ebic()],
#'   [compute_rbic()], [ic_step()]
#'
#' @references Akaike, H. (1974). A new look at the statistical model
#'   identification. *IEEE Transactions on Automatic Control*, **19**(6),
#'   716–723.
#'
#' @examples
#' fit <- lm(mpg ~ wt + hp, data = mtcars)
#' compute_aic(fit)
#'
#' fit_glm <- glm(am ~ wt + hp, data = mtcars, family = binomial)
#' compute_aic(fit_glm)
#'
#' @export
compute_aic <- function(fit, ...) UseMethod("compute_aic")

#' @export
compute_aic.lm <- function(fit, ...) {
  val <- stats::AIC(fit)
  k   <- length(stats::coef(fit))
  .ic_structure(val, fit_class = class(fit)[1], k = k, criterion = "AIC")
}

#' @export
compute_aic.glmnet <- function(fit, ...) {
  .check_glmnet_alpha(fit)
  k      <- .extract_k_glmnet(fit)
  loglik <- .glmnet_loglik(fit)
  val    <- -2 * loglik + 2 * k
  names(val) <- paste0("s", seq_along(fit$lambda) - 1L)
  .ic_structure(val,
    fit_class = "glmnet", k = k, criterion = "AIC",
    extras = list(lambda = fit$lambda)
  )
}

#' @export
compute_aic.ncvreg <- function(fit, ...) {
  val <- stats::AIC(fit)   # relies on logLik.ncvreg
  k   <- .extract_k_ncvreg(fit)
  .ic_structure(val,
    fit_class = "ncvreg", k = k, criterion = "AIC",
    extras = list(lambda = fit$lambda)
  )
}

#' @export
compute_aic.coxph <- function(fit, ...) {
  val <- stats::AIC(fit)
  k   <- .extract_k_coxph(fit)
  .ic_structure(val, fit_class = "coxph", k = k, criterion = "AIC")
}

#' @export
compute_aic.ncvsurv <- function(fit, ...) {
  .check_ncvreg_ns()
  val <- stats::AIC(fit)   # relies on logLik.ncvsurv
  k   <- .extract_k_ncvsurv(fit)
  .ic_structure(val,
    fit_class = "ncvsurv", k = k, criterion = "AIC",
    extras = list(lambda = fit$lambda)
  )
}
