#' Stepwise model selection using information criteria
#'
#' Performs forward, backward, or bidirectional stepwise model selection for
#' `lm`, `glm`, and `coxph` objects using any of the information criteria
#' provided by ICtoolkit: AIC, AICc, BIC, HQIC, EBIC, RBIC, mBIC, or mBIC2.
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
#' @param object A fitted `lm`, `glm`, or `coxph` object (the starting model).
#'   For `coxph`, the response should be a `Surv()` object; a null Cox model
#'   can be specified as `coxph(Surv(time, status) ~ 1, data = ...)`.
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
#' @param cl Cluster for parallel evaluation of candidates at each step.
#'   Either an integer (number of cores) or a `parallel::makeCluster()` object.
#'   When an integer is given on Unix/macOS, `parallel::mclapply()` is used
#'   (fork-based, zero-copy); on Windows a temporary PSOCK cluster is created.
#'   A user-supplied cluster object is used with `parallel::parLapply()` on
#'   all platforms.
#'   If the \pkg{pbapply} package is installed, a progress bar is shown
#'   when using a cluster object.
#'
#'   **Windows note:** PSOCK clusters must serialise data to each worker over
#'   sockets, so parallelisation overhead may outweigh the benefit for fast
#'   model fits (e.g., small `lm` problems). On Windows, parallel evaluation
#'   is most effective with expensive per-candidate fits such as large-\eqn{n}
#'   `glm` or `coxph` models, or when the number of candidates is very large.
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
#' ## Forward BIC selection for Cox models
#' \dontrun{
#' library(survival)
#' fit0 <- coxph(Surv(time, status) ~ 1, data = lung)
#' scope_upper <- Surv(time, status) ~ age + sex + ph.ecog + wt.loss
#' ic_step(fit0,
#'         scope     = list(lower = ~ 1, upper = scope_upper),
#'         direction = "forward",
#'         criterion = "BIC")
#' }
#'
#' @importFrom stats update getCall
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
                    cl        = NULL,
                    trace     = 1L,
                    steps     = 1000L) {

  direction  <- match.arg(direction)
  criterion  <- match.arg(criterion)
  caller_env <- parent.frame()

  # ---- set up parallel / progress-bar back-end ------------------------------
  mc_cores    <- NULL
  temp_cluster <- FALSE
  if (!is.null(cl)) {
    if (is.numeric(cl)) {
      if (.Platform$OS.type == "unix") {
        ## On Unix/macOS, mclapply forks and shares memory — no serialisation.
        mc_cores <- as.integer(cl)
        cl       <- NULL
      } else {
        cl <- parallel::makeCluster(as.integer(cl))
        temp_cluster <- TRUE
      }
    } else if (!inherits(cl, "cluster")) {
      stop("'cl' must be an integer (number of cores) or a 'cluster' object ",
           "from parallel::makeCluster().", call. = FALSE)
    }
  }

  if (temp_cluster) on.exit(parallel::stopCluster(cl), add = TRUE)

  has_pbapply <- requireNamespace("pbapply", quietly = TRUE)
  use_cluster <- !is.null(cl)

  # ---- validate criterion-specific arguments --------------------------------
  if (criterion == "EBIC" && is.null(P))
    stop("'P' (total candidate predictors) must be supplied for criterion = 'EBIC'.")
  if (criterion == "RBIC" && is.null(P_index))
    stop("'P_index' (named list of variable groups) must be supplied for criterion = 'RBIC'.")
  if (criterion %in% c("mBIC", "mBIC2") && is.null(P))
    stop("'P' (total candidate predictors) must be supplied for criterion = '", criterion, "'.")

  # ---- parse scope ----------------------------------------------------------
  # scope$lower / scope$upper (or a bare scope) may be a formula OR a
  # character vector of term labels.  The latter avoids R's formula parser,
  # which hits a stack overflow when the number of terms exceeds ~10 000.
  .terms_from_scope <- function(x) {
    if (is.character(x)) return(x)
    attr(stats::terms(stats::formula(x)), "term.labels")
  }

  lower_terms <- character(0)
  upper_terms <- attr(stats::terms(object), "term.labels")

  if (!missing(scope) && !is.null(scope)) {
    if (is.list(scope)) {
      if (!is.null(scope$lower))
        lower_terms <- .terms_from_scope(scope$lower)
      if (!is.null(scope$upper))
        upper_terms <- .terms_from_scope(scope$upper)
    } else {
      upper_terms <- .terms_from_scope(scope)
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

  ## Set up PSOCK cluster workers (only when a cluster object is used).
  ## On Unix with an integer cl, we use mclapply instead (no setup needed).
  if (use_cluster) {
    ## Build a minimal environment for workers containing only the data
    ## referenced by the model call (not the caller's full frame).
    model_call  <- getCall(object)
    .worker_env <- new.env(parent = globalenv())
    data_args   <- c("data", "weights", "subset", "offset")
    for (da in data_args) {
      sym <- model_call[[da]]
      if (!is.null(sym) && is.symbol(sym))
        .worker_env[[as.character(sym)]] <- eval(sym, caller_env)
    }

    parallel::clusterExport(cl,
      c(".worker_env", "criterion", "P", "P_index", "kappa", "gamma"),
      envir = environment())
    parallel::clusterCall(cl, function() {
      requireNamespace("ICtoolkit", quietly = TRUE)
      requireNamespace("survival", quietly = TRUE)
    })

    ## Lightweight worker function for PSOCK cluster evaluation.
    ## Created with globalenv() as its environment so that serialisation sends
    ## only the function body (a reference), not ic_step()'s entire frame.
    ## On workers, .worker_env / criterion / P / etc. are found in .GlobalEnv
    ## (placed there by clusterExport above), and ICtoolkit:: qualifies the
    ## compute functions so they resolve via the loaded namespace.
    .par_eval <- eval(quote(function(ucall) {
      cfit <- eval(ucall, .worker_env)
      as.numeric(switch(criterion,
        AIC   = ICtoolkit::compute_aic(cfit),
        AICc  = ICtoolkit::compute_aicc(cfit),
        BIC   = ICtoolkit::compute_bic(cfit),
        HQIC  = ICtoolkit::compute_hqic(cfit),
        EBIC  = ICtoolkit::compute_ebic(cfit, P = P, gamma = gamma),
        RBIC  = ICtoolkit::compute_rbic(cfit, P_index = P_index, gamma = gamma),
        mBIC  = ICtoolkit::compute_mbic(cfit, P = P, kappa = kappa),
        mBIC2 = ICtoolkit::compute_mbic2(cfit, P = P, kappa = kappa)
      ))
    }), envir = globalenv())
  }

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
    if (trace > 0L && length(candidates) > 100L)
      message(sprintf("  Evaluating %d candidates...", length(candidates)))

    ## Build the unevaluated update calls (lightweight, no data copying)
    update_calls <- lapply(candidates, function(fml) {
      update(current, fml, evaluate = FALSE)
    })

    ## Evaluate all candidates (parallel or sequential).
    ## When pbapply is available a progress bar with ETA is shown.
    .seq_eval <- function(ucall) .ic(eval(ucall, caller_env))

    if (!is.null(mc_cores)) {
      ## Fork-based (Unix): no serialisation, workers inherit parent memory.
      ## pbapply::pblapply with an integer cl uses mclapply internally on Unix.
      if (has_pbapply) {
        ic_values <- unlist(pbapply::pblapply(update_calls, .seq_eval,
                                              cl = mc_cores))
      } else {
        ic_values <- unlist(parallel::mclapply(update_calls, .seq_eval,
                                               mc.cores = mc_cores))
      }
    } else if (use_cluster) {
      ## PSOCK / user-supplied cluster
      if (has_pbapply) {
        ic_values <- unlist(pbapply::pblapply(update_calls, .par_eval, cl = cl))
      } else {
        ic_values <- unlist(parallel::parLapply(cl, update_calls, .par_eval))
      }
    } else {
      ## Sequential
      if (has_pbapply) {
        ic_values <- unlist(pbapply::pblapply(update_calls, .seq_eval))
      } else {
        ic_values <- unlist(lapply(update_calls, .seq_eval))
      }
    }

    if (trace > 1L) {
      for (j in seq_along(ic_values))
        message(sprintf("  %-40s %s = %.4f",
                        names(candidates)[j], criterion, ic_values[j]))
    }

    ## Find the best candidate
    best_idx <- which.min(ic_values)
    best_ic  <- ic_values[best_idx]

    if (best_ic >= current_ic) {
      best_name <- NULL
    } else {
      best_name <- names(candidates)[best_idx]
      ## Refit only the winner locally (workers may not return the fit object)
      best_fit <- eval(update_calls[[best_idx]], caller_env)
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
