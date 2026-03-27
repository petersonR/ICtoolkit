#' Compute RBIC for model objects
#'
#' Computes the Regularized Bayesian Information Criterion (RBIC), an
#' extension of EBIC that accounts for structured (grouped) predictors.
#' RBIC applies group-specific EBIC penalties:
#'
#' \deqn{\mathrm{RBIC} = \mathrm{BIC} + 2 \sum_{g=1}^{G}
#'   \gamma_g \log\binom{P_g}{k_g}}
#'
#' where `G` is the number of groups, `P_g` is the total number of candidate
#' predictors in group `g`, `k_g` is the number selected from group `g`, and
#' `gamma_g` is a group-level tuning parameter.
#'
#' @section P_index format:
#' `P_index` is a named list whose elements identify which predictors belong
#' to each group:
#'
#' * **`lm`/`glm`**: character vectors of term labels as returned by
#'   `attr(terms(fit), "term.labels")` (e.g. `"x1"`, `"I(x1^2)"`,
#'   `"x1:x2"`).
#' * **`glmnet`**: character vectors matching `rownames(fit$beta)`.
#' * **`ncvreg`**: character vectors matching `rownames(fit$beta[-1, ])`.
#'
#' @section The `gamma` argument:
#' `gamma` controls how the group-level tuning parameters are determined.
#' Accepted values:
#' \describe{
#'   \item{`"ebic"` (default)}{Single gamma based on the **total** number of
#'     candidates `P_total = sum(p_1, ..., p_K)`, applied uniformly across
#'     all groups.  Matches the paper's default formula:
#'     \deqn{\gamma = \frac{\log(P_{total}/\sqrt{n})}{\log P_{total}}
#'       \cdot \mathbf{1}(P_{total} \ge \sqrt{n})}}
#'   \item{`"per_group"`}{Separate gamma for each group, each computed from
#'     its own `p_k`.  Stronger penalisation of larger groups.}
#'   \item{A numeric scalar}{Fixed gamma applied uniformly to all groups,
#'     e.g. `gamma = 0.5`.}
#'   \item{A numeric vector of length `K`}{Fixed group-specific gammas.}
#'   \item{A function `f(P, k, n)`}{Advanced: called once with the full
#'     vectors of group sizes and selected counts; must return a scalar or
#'     length-`K` vector.}
#' }
#'
#' @param fit A fitted model object of class `lm`, `glm`, `glmnet`, or
#'   `ncvreg`.
#' @param P_index A named list of character vectors defining variable groups.
#'   See the P_index format section.
#' @param gamma Gamma specification.  See the `gamma` argument section.
#'   Default `"ebic"`.
#' @param ... Currently unused; reserved for future use.
#'
#' @return A named numeric scalar (`lm`/`glm`) or vector (`glmnet`/`ncvreg`)
#'   with the following attributes:
#'   \describe{
#'     \item{`fit_class`}{Character. The primary class of `fit`.}
#'     \item{`k`}{Integer scalar or vector. Total estimated parameters
#'       (including intercept).}
#'     \item{`gamma`}{The computed gamma value(s): a scalar or vector for
#'       `lm`/`glm`; a list (one entry per lambda) for `glmnet`/`ncvreg`.}
#'     \item{`P_index`}{The `P_index` supplied by the user.}
#'     \item{`criterion`}{`"RBIC"`.}
#'   }
#'
#' @seealso [compute_ebic()], [compute_bic()], [ic_step()]
#'
#' @references Chen, J., & Chen, Z. (2008). Extended Bayesian information
#'   criteria for model selection with large model spaces. *Biometrika*,
#'   **95**(3), 759–771.
#'
#' @examples
#' set.seed(1)
#' n  <- 100
#' x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
#' y  <- 2 * x1 - x2 + rnorm(n)
#' df <- data.frame(y, x1, x2, x3, x1sq = x1^2, x2sq = x2^2)
#' fit <- lm(y ~ x1 + x2 + x3 + x1sq + x2sq, data = df)
#'
#' P_index <- list(main = c("x1", "x2", "x3"), squared = c("x1sq", "x2sq"))
#'
#' compute_rbic(fit, P_index)                          # default (total-P gamma)
#' compute_rbic(fit, P_index, gamma = "per_group")     # per-group gammas
#' compute_rbic(fit, P_index, gamma = 0.5)             # fixed gamma
#' compute_rbic(fit, P_index, gamma = 0)               # gamma=0 recovers BIC
#'
#' @export
compute_rbic <- function(fit, P_index, gamma = "ebic", ...) UseMethod("compute_rbic")

# ---------------------------------------------------------------------------
# Internal helpers shared across methods
# ---------------------------------------------------------------------------

# Validate that P_index was supplied.
.check_P_index <- function(P_index) {
  if (missing(P_index) || is.null(P_index))
    stop("'P_index' (named list of variable groups) must be supplied.")
}

# Compute the RBIC penalty:
#   P_g     : integer vector of all group sizes  [p_1, ..., p_K]
#   k_g     : integer vector of selected counts  [m_1, ..., m_K]
#   n       : sample size
#   gammafn : resolved function; called once with (P_g, k_g, n)
# Returns list(gamma = computed value(s), penalty = scalar).
.rbic_penalty <- function(P_g, k_g, n, gammafn) {
  gval    <- gammafn(P_g, k_g, n)           # scalar or length-K vector
  lc_g    <- mapply(.log_choose, P_g, k_g)  # length-K vector
  penalty <- 2 * sum(gval * lc_g)
  list(gamma = gval, penalty = penalty)
}

# Count selected predictors per group from a character vector of active terms.
.count_selected <- function(P_index, active_terms) {
  vapply(P_index, function(grp) sum(grp %in% active_terms), integer(1))
}

# ---------------------------------------------------------------------------
# lm method
# ---------------------------------------------------------------------------

#' @export
compute_rbic.lm <- function(fit, P_index, gamma = "ebic", ...) {
  .check_P_index(P_index)
  gammafn <- .resolve_gammafn(gamma)

  n   <- stats::nobs(fit)
  bic <- stats::BIC(fit)
  k   <- length(stats::coef(fit))

  active_terms <- attr(stats::terms(fit), "term.labels")
  P_g <- vapply(P_index, length, integer(1))
  k_g <- .count_selected(P_index, active_terms)

  res <- .rbic_penalty(P_g, k_g, n, gammafn)
  val <- bic + res$penalty

  .ic_structure(val,
    fit_class = class(fit)[1], k = k, criterion = "RBIC",
    extras = list(gamma = res$gamma, P_index = P_index)
  )
}

# ---------------------------------------------------------------------------
# glmnet method
# ---------------------------------------------------------------------------

#' @export
compute_rbic.glmnet <- function(fit, P_index, gamma = "ebic", ...) {
  if (missing(P_index) || is.null(P_index))
    stop("'P_index' (named list of variable groups) must be supplied.")
  .check_glmnet_alpha(fit)
  gammafn <- .resolve_gammafn(gamma)

  n       <- fit$nobs
  bic     <- as.numeric(compute_bic(fit))
  k_total <- .extract_k_glmnet(fit)
  nlambda <- length(fit$lambda)

  beta      <- as.matrix(fit$beta)
  var_names <- rownames(beta)

  bad <- setdiff(unlist(P_index), var_names)
  if (length(bad) > 0)
    stop("P_index contains variables not found in rownames(fit$beta): ",
         paste(bad, collapse = ", "))

  P_g <- vapply(P_index, length, integer(1))

  val        <- numeric(nlambda)
  gamma_list <- vector("list", nlambda)

  for (l in seq_len(nlambda)) {
    active          <- var_names[beta[, l] != 0]
    k_g             <- .count_selected(P_index, active)
    res             <- .rbic_penalty(P_g, k_g, n, gammafn)
    val[l]          <- bic[l] + res$penalty
    gamma_list[[l]] <- res$gamma
  }

  names(val) <- paste0("s", seq_len(nlambda) - 1L)

  .ic_structure(val,
    fit_class = "glmnet", k = k_total, criterion = "RBIC",
    extras = list(gamma = gamma_list, P_index = P_index, lambda = fit$lambda)
  )
}

# ---------------------------------------------------------------------------
# ncvreg method
# ---------------------------------------------------------------------------

#' @export
compute_rbic.ncvreg <- function(fit, P_index, gamma = "ebic", ...) {
  .check_P_index(P_index)
  gammafn <- .resolve_gammafn(gamma)

  n       <- length(fit$y)
  bic     <- as.numeric(compute_bic(fit))
  k_total <- .extract_k_ncvreg(fit)
  nlambda <- length(fit$lambda)

  beta      <- fit$beta[-1, , drop = FALSE]
  var_names <- rownames(beta)

  bad <- setdiff(unlist(P_index), var_names)
  if (length(bad) > 0)
    stop("P_index contains variables not found in rownames(fit$beta[-1,]): ",
         paste(bad, collapse = ", "))

  P_g <- vapply(P_index, length, integer(1))

  val        <- numeric(nlambda)
  gamma_list <- vector("list", nlambda)

  for (l in seq_len(nlambda)) {
    active          <- var_names[beta[, l] != 0]
    k_g             <- .count_selected(P_index, active)
    res             <- .rbic_penalty(P_g, k_g, n, gammafn)
    val[l]          <- bic[l] + res$penalty
    gamma_list[[l]] <- res$gamma
  }

  .ic_structure(val,
    fit_class = "ncvreg", k = k_total, criterion = "RBIC",
    extras = list(gamma = gamma_list, P_index = P_index, lambda = fit$lambda)
  )
}
