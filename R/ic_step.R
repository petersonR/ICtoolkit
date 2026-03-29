#' Stepwise model selection using information criteria
#'
#' Performs forward, backward, or bidirectional stepwise model selection for
#' `lm` and `glm` objects using any of the information criteria provided
#' by ICtoolkit: AIC, AICc, BIC, HQIC, EBIC, RBIC, mBIC, or mBIC2.
#'
#' At each step, all single-term additions (forward) and/or deletions
#' (backward) are evaluated.  The move that gives the greatest IC reduction
#' is accepted.  Stepping halts when no candidate improves the current IC.
#'
#' @section Scope:
#' `scope` follows the same convention as [MASS::stepAIC()]:
#' * A single formula specifies the **upper** bound (largest allowable model).
#' * A list with elements `lower` and/or `upper` specifies both bounds.
#' The starting model is always `object`; terms in `lower` can never be
#' dropped.
#'
#' @section P_index format for RBIC:
#' `P_index` is a named list whose elements are character vectors of term
#' labels **exactly as they appear** in `attr(terms(fit), "term.labels")`.
#' For example, if the upper model contains `x1`, `I(x1^2)`, and `x1:x2`,
#' a valid `P_index` might be:
#' ```r
#' list(main = "x1", squared = "I(x1^2)", interactions = "x1:x2")
#' ```
#' All terms in the upper-model scope that are not covered by any group in
#' `P_index` receive no additional RBIC penalty (their BIC contribution still
#' applies).
#'
#' @param object A fitted `lm` or `glm` object (the starting model).
#' @param scope A formula or list with `lower` and/or `upper` elements
#'   specifying the model search space.  See [MASS::stepAIC()] for details.
#' @param direction One of `"both"` (default), `"backward"`, or `"forward"`.
#' @param criterion One of `"AIC"` (default), `"AICc"`, `"BIC"`, `"HQIC"`,
#'   `"EBIC"`, `"RBIC"`, `"mBIC"`, or `"mBIC2"`.
#' @param P Integer. Total number of candidate predictors.  Required only
#'   for `criterion = "EBIC"`.
#' @param P_index Named list of character vectors defining variable groups.
#'   Required only for `criterion = "RBIC"`.  See Details.
#' @param kappa Positive numeric scalar controlling the prior expected number
#'   of true signals.  Used only for `criterion = "mBIC"` or `"mBIC2"`.
#'   Default `4`.  See [compute_mbic()].
#' @param gamma Gamma specification passed to [compute_ebic()] or
#'   [compute_rbic()].  Ignored for AIC / AICc / BIC / HQIC / mBIC / mBIC2.
#'   See those functions for accepted values (`"ebic"`, `"per_group"`, numeric,
#'   or a function).
#' @param trace Integer controlling verbosity.  `0` = silent, `1` = one line
#'   per accepted step (default), `2` = all candidates at each step.
#' @param steps Maximum number of steps.  Default is `1000`.
#'
#' @return The final selected model (same class as `object`), with an
#'   additional attribute `"step_path"`: a data frame recording each accepted
#'   step with columns `step`, `action` (the term added/dropped), and the IC
#'   value after that step.
#'
#' @seealso [compute_aic()], [compute_aicc()], [compute_bic()],
#'   [compute_hqic()], [compute_ebic()], [compute_rbic()],
#'   [compute_mbic()], [compute_mbic2()], [MASS::stepAIC()]
#'
#' @examples
#' ## Backward BIC selection
#' fit_full <- lm(mpg ~ ., data = mtcars)
#' ic_step(fit_full, direction = "backward", criterion = "BIC")
#'
#' ## Forward EBIC selection
#' fit_null <- lm(mpg ~ 1, data = mtcars)
#' fit_full <- lm(mpg ~ ., data = mtcars)
#' ic_step(fit_null,
#'         scope     = list(lower = fit_null, upper = fit_full),
#'         direction = "forward",
#'         criterion = "EBIC",
#'         P         = ncol(mtcars) - 1L)
#'
#' @importFrom stats update
#' @export
ic_step <- function(object,
                    scope,
                    direction = c("both", "backward", "forward"),
                    criterion = c("AIC", "AICc", "BIC", "HQIC", "EBIC", "RBIC",
                                  "mBIC", "mBIC2"),
                    P         = NULL,
                    P_index   = NULL,
                    kappa     = 4,
                    gamma     = "ebic",
                    trace     = 1L,
                    steps     = 1000L) {

  direction  <- match.arg(direction)
  criterion  <- match.arg(criterion)
  caller_env <- parent.frame()

  # ---- validate criterion-specific arguments --------------------------------
  if (criterion == "EBIC" && is.null(P))
    stop("'P' (total candidate predictors) must be supplied for criterion = 'EBIC'.")
  if (criterion == "RBIC" && is.null(P_index))
    stop("'P_index' (named list of variable groups) must be supplied for criterion = 'RBIC'.")
  if (criterion %in% c("mBIC", "mBIC2") && is.null(P))
    stop("'P' (total candidate predictors) must be supplied for criterion = '", criterion, "'.")

  # ---- parse scope ----------------------------------------------------------
  lower_terms <- character(0)
  upper_terms <- attr(stats::terms(object), "term.labels")

  if (!missing(scope) && !is.null(scope)) {
    if (is.list(scope)) {
      if (!is.null(scope$lower))
        lower_terms <- attr(stats::terms(stats::formula(scope$lower)), "term.labels")
      if (!is.null(scope$upper))
        upper_terms <- attr(stats::terms(stats::formula(scope$upper)), "term.labels")
    } else {
      upper_terms <- attr(stats::terms(stats::formula(scope)), "term.labels")
    }
  }

  # ---- IC evaluation helper -------------------------------------------------
  .ic <- function(fit) {
    as.numeric(switch(criterion,
      AIC   = compute_aic(fit),
      AICc  = compute_aicc(fit),
      BIC   = compute_bic(fit),
      HQIC  = compute_hqic(fit),
      EBIC  = compute_ebic(fit, P = P, gamma = gamma),
      RBIC  = compute_rbic(fit, P_index = P_index, gamma = gamma),
      mBIC  = compute_mbic(fit, P = P, kappa = kappa),
      mBIC2 = compute_mbic2(fit, P = P, kappa = kappa)
    ))
  }

  # ---- helper: rebuild formula from term labels ----------------------------
  # Avoids paste-based manipulation which breaks terms like I(x^2).
  .make_formula <- function(terms_vec, base_fit) {
    resp      <- deparse(stats::formula(base_fit)[[2]])
    intercept <- as.logical(attr(stats::terms(base_fit), "intercept"))
    if (length(terms_vec) == 0L) {
      stats::formula(paste(resp, "~ 1"))
    } else {
      stats::reformulate(terms_vec, response = resp, intercept = intercept)
    }
  }

  # ---- initialise -----------------------------------------------------------
  current    <- object
  current_ic <- .ic(current)
  path_rows  <- list(list(step = 0L, action = "<start>", ic = current_ic))

  if (trace > 0L)
    message(sprintf("Start:  %s = %.4f", criterion, current_ic))

  # ---- main loop ------------------------------------------------------------
  for (step in seq_len(steps)) {

    current_terms <- attr(stats::terms(current), "term.labels")
    droppable     <- setdiff(current_terms, lower_terms)
    addable       <- setdiff(upper_terms,   current_terms)

    candidates <- list()

    if (direction %in% c("backward", "both")) {
      for (trm in droppable) {
        new_terms <- setdiff(current_terms, trm)
        candidates[[paste0("- ", trm)]] <- .make_formula(new_terms, current)
      }
    }
    if (direction %in% c("forward", "both")) {
      for (trm in addable) {
        new_terms <- c(current_terms, trm)
        candidates[[paste0("+ ", trm)]] <- .make_formula(new_terms, current)
      }
    }

    if (length(candidates) == 0L) {
      if (trace > 0L) message("No candidates available. Stopping.")
      break
    }

    # evaluate all candidates
    best_ic   <- current_ic
    best_name <- NULL
    best_fit  <- NULL

    for (cname in names(candidates)) {
      cfit <- eval(update(current, candidates[[cname]], evaluate = FALSE), caller_env)
      cic  <- .ic(cfit)

      if (trace > 1L)
        message(sprintf("  %-40s %s = %.4f", cname, criterion, cic))

      if (cic < best_ic) {
        best_ic   <- cic
        best_name <- cname
        best_fit  <- cfit
      }
    }

    if (is.null(best_name)) {
      if (trace > 0L)
        message(sprintf("No improvement at step %d. Stopping.", step))
      break
    }

    if (trace > 0L)
      message(sprintf("Step %2d: %-35s %s = %.4f",
                      step, best_name, criterion, best_ic))

    current    <- best_fit
    current_ic <- best_ic
    path_rows[[length(path_rows) + 1L]] <- list(step = step, action = best_name,
                                                 ic = best_ic)
  }

  path <- do.call(rbind.data.frame, path_rows)
  names(path) <- c("step", "action", criterion)
  attr(current, "step_path") <- path
  current
}
