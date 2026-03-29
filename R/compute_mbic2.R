#' Compute mBIC2 for model objects
#'
#' Computes the modified BIC2 criterion (mBIC2) of Frommlet & Bogdan (2013).
#' mBIC2 extends mBIC by subtracting a combinatorial correction that accounts
#' for the number of models of a given size:
#'
#' \deqn{\mathrm{mBIC2} = \mathrm{BIC} + 2\,k_{\text{pred}}\,
#'   \log\!\bigl(P/\kappa - 1\bigr) - 2\,\log\binom{P}{k_{\text{pred}}}}
#'
#' The additional `log C(P, k)` term partially offsets the mBIC penalty,
#' making mBIC2 less conservative than mBIC when multiple predictors are
#' selected.
#'
#' For `glmnet` and `ncvreg`, mBIC2 is computed at each point on the
#' regularization path (one value per lambda).
#'
#' @inheritParams compute_mbic
#'
#' @return A named numeric scalar (`lm`/`glm`) or vector (`glmnet`/`ncvreg`)
#'   with the following attributes:
#'   \describe{
#'     \item{`fit_class`}{Character. The primary class of `fit`.}
#'     \item{`k`}{Integer scalar or vector. Number of estimated parameters
#'       (including intercept).}
#'     \item{`P`}{Integer. The value of `P` used.}
#'     \item{`kappa`}{Numeric. The value of `kappa` used.}
#'     \item{`criterion`}{`"mBIC2"`.}
#'   }
#'
#' @seealso [compute_mbic()], [compute_ebic()], [compute_bic()], [ic_step()]
#'
#' @references Frommlet, F., & Bogdan, M. (2013). Some optimality properties
#'   of FDR controlling rules under sparsity. *Electronic Journal of
#'   Statistics*, **7**, 1328--1368.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- X[, 1:3] %*% c(2, -1, 0.5) + rnorm(100)
#' fit <- lm(y ~ X)
#'
#' compute_mbic2(fit, P = 20)                 # default kappa = 4
#' compute_mbic2(fit, P = 20, kappa = 2)      # expect fewer signals
#'
#' # Compare mBIC and mBIC2
#' compute_mbic(fit, P = 20)                  # more conservative
#' compute_mbic2(fit, P = 20)                 # less conservative
#'
#' @export
compute_mbic2 <- function(fit, P = NULL, kappa = 4, ...) UseMethod("compute_mbic2")

#' @export
compute_mbic2.lm <- function(fit, P = NULL, kappa = 4, ...) {
  k_total <- length(stats::coef(fit))
  if (is.null(P)) P <- k_total - 1L
  .check_kappa(P, kappa)

  bic    <- stats::BIC(fit)
  k_pred <- k_total - 1L
  val    <- bic + .mbic2_penalty(k_pred, P, kappa)

  .ic_structure(val,
    fit_class = class(fit)[1], k = k_total, criterion = "mBIC2",
    extras = list(P = P, kappa = kappa)
  )
}

#' @export
compute_mbic2.glmnet <- function(fit, P = NULL, kappa = 4, ...) {
  if (is.null(P)) P <- nrow(fit$beta)
  .check_kappa(P, kappa)
  .check_glmnet_alpha(fit)

  bic     <- as.numeric(compute_bic(fit))
  k_total <- .extract_k_glmnet(fit)
  k_pred  <- k_total - 1L
  penalty <- .mbic2_penalty(k_pred, P, kappa)
  val     <- bic + penalty
  names(val) <- paste0("s", seq_along(fit$lambda) - 1L)

  .ic_structure(val,
    fit_class = "glmnet", k = k_total, criterion = "mBIC2",
    extras = list(P = P, kappa = kappa, lambda = fit$lambda)
  )
}

#' @export
compute_mbic2.ncvreg <- function(fit, P = NULL, kappa = 4, ...) {
  if (is.null(P)) P <- nrow(fit$beta) - 1L
  .check_kappa(P, kappa)

  bic     <- as.numeric(stats::BIC(fit))
  k_total <- .extract_k_ncvreg(fit)
  k_pred  <- k_total - 1L
  penalty <- .mbic2_penalty(k_pred, P, kappa)
  val     <- bic + penalty

  .ic_structure(val,
    fit_class = "ncvreg", k = k_total, criterion = "mBIC2",
    extras = list(P = P, kappa = kappa, lambda = fit$lambda)
  )
}

#' @export
compute_mbic2.coxph <- function(fit, P = NULL, kappa = 4, ...) {
  k_total <- .extract_k_coxph(fit)
  if (is.null(P)) P <- k_total
  .check_kappa(P, kappa)

  bic    <- stats::BIC(fit)
  k_pred <- k_total
  val    <- bic + .mbic2_penalty(k_pred, P, kappa)

  .ic_structure(val,
    fit_class = "coxph", k = k_total, criterion = "mBIC2",
    extras = list(P = P, kappa = kappa)
  )
}

#' @export
compute_mbic2.ncvsurv <- function(fit, P = NULL, kappa = 4, ...) {
  if (is.null(P)) P <- nrow(fit$beta)
  .check_kappa(P, kappa)

  bic     <- as.numeric(compute_bic(fit))
  k_pred  <- .extract_k_ncvsurv(fit)
  penalty <- .mbic2_penalty(k_pred, P, kappa)
  val     <- bic + penalty

  .ic_structure(val,
    fit_class = "ncvsurv", k = k_pred, criterion = "mBIC2",
    extras = list(P = P, kappa = kappa, lambda = fit$lambda)
  )
}
