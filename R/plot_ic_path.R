#' Plot information criterion paths for glmnet or ncvreg objects
#'
#' Computes one or more information criteria across the full regularisation
#' path of a `glmnet` or `ncvreg` object and plots each as a line against
#' a transform of lambda.  A vertical dashed line marks the lambda that
#' minimises each criterion.
#'
#' @param fit A fitted `glmnet` or `ncvreg` object.
#' @param criteria Character vector of criteria to plot.  One or more of
#'   `"AIC"`, `"AICc"`, `"BIC"`, `"EBIC"`, `"RBIC"`.  Default `c("AIC", "BIC")`.
#' @param P Integer. Total candidate predictors (excluding intercept).
#'   Required when `"EBIC"` is in `criteria`.
#' @param P_index Named list of variable groups.  Required when `"RBIC"` is
#'   in `criteria`.  See [compute_rbic()] for format details.
#' @param gamma Gamma specification for EBIC / RBIC.  Default `"ebic"`.
#'   See [compute_ebic()] for accepted values.
#' @param xvar What to place on the x-axis.  One of `"-log(lambda)"` (default),
#'   `"log(lambda)"`, or `"lambda"`.
#' @param mark_min Logical.  Draw a vertical dashed line at the
#'   IC-minimising lambda for each criterion.  Default `TRUE`.
#' @param legend Logical.  Draw a legend when more than one criterion is
#'   plotted.  Default `TRUE`.
#' @param cols Named or unnamed character/integer vector of colours, one per
#'   criterion.  Defaults to `seq_along(criteria)` (R's standard palette).
#' @param lwd Line width.  Default `2`.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns a list with:
#'   \describe{
#'     \item{`ic`}{Named list of IC vectors, one element per criterion.}
#'     \item{`opt_lambda`}{Named numeric vector of IC-minimising lambda values.}
#'     \item{`lambda`}{The full lambda sequence from `fit`.}
#'   }
#'
#' @seealso [compute_aic()], [compute_bic()], [compute_ebic()],
#'   [compute_rbic()]
#'
#' @examples
#' \dontrun{
#' library(glmnet)
#' set.seed(1)
#' X <- matrix(rnorm(100 * 20), 100, 20,
#'             dimnames = list(NULL, paste0("x", 1:20)))
#' y <- X[, 1] * 2 - X[, 3] + rnorm(100)
#' fit <- glmnet(X, y)
#'
#' plot_ic_path(fit, criteria = c("AIC", "BIC", "EBIC"), P = 20)
#'
#' # RBIC with two predictor groups
#' P_index <- list(first10 = paste0("x", 1:10), second10 = paste0("x", 11:20))
#' plot_ic_path(fit, criteria = c("BIC", "RBIC"), P_index = P_index)
#' }
#'
#' @export
plot_ic_path <- function(fit,
                          criteria  = c("AIC", "BIC"),
                          P         = NULL,
                          P_index   = NULL,
                          gamma     = "ebic",
                          xvar      = c("-log(lambda)", "log(lambda)", "lambda"),
                          mark_min  = TRUE,
                          legend    = TRUE,
                          cols      = NULL,
                          lwd       = 2,
                          ...) {

  if (!inherits(fit, c("glmnet", "ncvreg")))
    stop("'fit' must be a glmnet or ncvreg object.")

  valid_criteria <- c("AIC", "AICc", "BIC", "EBIC", "RBIC")
  bad <- setdiff(criteria, valid_criteria)
  if (length(bad))
    stop("Unknown criteria: ", paste(bad, collapse = ", "),
         ". Must be one of: ", paste(valid_criteria, collapse = ", "))

  if ("RBIC" %in% criteria && is.null(P_index))
    stop("'P_index' (named list of variable groups) is required for RBIC.")

  xvar   <- match.arg(xvar)
  lambda <- fit$lambda

  # ---- compute each criterion -----------------------------------------------
  ic_vals <- lapply(criteria, function(crit) {
    suppressWarnings(as.numeric(switch(crit,
      AIC  = compute_aic(fit),
      AICc = compute_aicc(fit),
      BIC  = compute_bic(fit),
      EBIC = compute_ebic(fit, P = P, gamma = gamma),   # P=NULL is OK: inferred
      RBIC = compute_rbic(fit, P_index = P_index, gamma = gamma)
    )))
  })
  names(ic_vals) <- criteria

  # ---- x-axis values --------------------------------------------------------
  x <- switch(xvar,
    "-log(lambda)" = -log(lambda),
    "log(lambda)"  =  log(lambda),
    "lambda"       =  lambda
  )

  # ---- colours --------------------------------------------------------------
  if (is.null(cols)) {
    cols <- seq_along(criteria)
  } else if (length(cols) < length(criteria)) {
    cols <- rep_len(cols, length(criteria))
  }
  cols <- stats::setNames(cols, criteria)

  # ---- optimal lambda per criterion -----------------------------------------
  opt_idx    <- sapply(ic_vals, which.min)
  opt_x      <- x[opt_idx]
  opt_lambda <- lambda[opt_idx]
  names(opt_lambda) <- criteria

  # ---- plot -----------------------------------------------------------------
  ylim <- range(unlist(ic_vals), na.rm = TRUE)

  graphics::plot(x, ic_vals[[1]],
    type = "l", col = cols[1], lwd = lwd,
    xlab = xvar, ylab = "IC value", ylim = ylim, ...
  )

  if (length(criteria) > 1) {
    for (i in seq_along(criteria)[-1])
      graphics::lines(x, ic_vals[[i]], col = cols[i], lwd = lwd)
  }

  if (mark_min) {
    for (i in seq_along(criteria))
      graphics::abline(v = opt_x[i], col = cols[i], lty = 2, lwd = 1)
  }

  if (legend && length(criteria) > 1)
    graphics::legend("topright", legend = criteria,
                     col = cols, lty = 1, lwd = lwd, bty = "n")

  invisible(list(ic = ic_vals, opt_lambda = opt_lambda, lambda = lambda))
}
