#' Expand a base formula with interactions and polynomials
#'
#' Given a base formula and data, generates an expanded formula containing
#' interaction terms (up to order `k`) and polynomial terms (up to degree
#' `poly`) for numeric predictors. Returns the expanded formula along with a
#' `P_index` list suitable for use with [compute_rbic()] or [ic_step()].
#'
#' Polynomial terms are added as `I(x^2)`, `I(x^3)`, etc. for each numeric
#' predictor in the base formula. Factor predictors are included in
#' interactions but do not receive polynomial terms. Interactions are formed
#' among the **base** terms only (polynomials do not interact with each other
#' or with base terms).
#'
#' @param formula A formula specifying the response and base (main-effect)
#'   predictors, e.g. `y ~ x1 + x2 + g`.
#' @param data A data frame used to determine which predictors are numeric
#'   (eligible for polynomials) vs. factor/character.
#' @param k Integer. Maximum interaction order among base terms.
#'   `k = 1` means no interactions (default), `k = 2` means pairwise, etc.
#' @param poly Integer. Maximum polynomial degree for numeric predictors.
#'   `poly = 1` means no polynomials (default), `poly = 2` adds squared
#'   terms, `poly = 3` adds squared and cubic, etc.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`formula`}{The expanded formula.}
#'     \item{`P_index`}{A named list of character vectors grouping term labels
#'       by type: `"main"` for base terms, `"poly_2"`, `"poly_3"`, ... for
#'       polynomial terms, and `"interact_2"`, `"interact_3"`, ... for
#'       interaction terms.  Only non-empty groups are included.}
#'   }
#'
#' @seealso [ic_step()], [compute_rbic()]
#'
#' @examples
#' es <- expand_scope(ies_total ~ age + race + hispanic + female,
#'                    data = covid_eol, k = 2, poly = 2)
#' es$formula
#' es$P_index
#'
#' # Use directly with ic_step
#' fit_null <- lm(ies_total ~ 1, data = covid_eol)
#' ic_step(fit_null, scope = list(lower = ~ 1, upper = es$formula),
#'         direction = "forward", criterion = "RBIC",
#'         P_index = es$P_index, trace = 0)
#'
#' @export
expand_scope <- function(formula, data, k = 1L, poly = 1L) {
  k    <- as.integer(k)
  poly <- as.integer(poly)
  if (k < 1L)    stop("'k' must be >= 1.")
  if (poly < 1L) stop("'poly' must be >= 1.")

  tt         <- stats::terms(formula, data = data)
  resp       <- deparse(tt[[2L]])
  base_terms <- attr(tt, "term.labels")

  if (length(base_terms) == 0L)
    stop("Formula must contain at least one predictor.")

  # Identify numeric predictors (eligible for polynomials)
  numeric_terms <- base_terms[vapply(base_terms, function(t) {
    is.numeric(data[[t]])
  }, logical(1))]

  # --- Build polynomial terms ---
  poly_terms <- list()
  if (poly >= 2L && length(numeric_terms) > 0L) {
    for (d in 2:poly) {
      poly_terms[[paste0("poly_", d)]] <- paste0("I(", numeric_terms, "^", d, ")")
    }
  }

  # --- Build interaction formula ---
  interact_terms <- list()
  if (k >= 2L && length(base_terms) >= 2L) {
    ia_formula <- stats::as.formula(
      paste("~", paste0("(", paste(base_terms, collapse = " + "), ")^", k))
    )
    ia_labels <- attr(stats::terms(ia_formula, data = data), "term.labels")
    # Group by interaction order (number of colons + 1)
    for (lab in ia_labels) {
      order <- lengths(regmatches(lab, gregexpr(":", lab))) + 1L
      if (order >= 2L) {
        grp <- paste0("interact_", order)
        interact_terms[[grp]] <- c(interact_terms[[grp]], lab)
      }
    }
  }

  # --- Assemble expanded formula ---
  all_extra <- c(unlist(poly_terms, use.names = FALSE),
                 unlist(interact_terms, use.names = FALSE))
  if (length(all_extra) > 0L) {
    expanded_formula <- stats::as.formula(
      paste(resp, "~", paste(c(base_terms, all_extra), collapse = " + "))
    )
  } else {
    expanded_formula <- formula
  }

  # --- Build P_index ---
  P_index <- c(list(main = base_terms), poly_terms, interact_terms)

  list(formula = expanded_formula, P_index = P_index)
}
