#' ICtoolkit: Information Criteria for Penalized and Stepwise Model Selection
#'
#' @description
#' ICtoolkit provides a unified interface for computing five information
#' criteria across four model classes:
#'
#' | Criterion | Function          | Notes                                      |
#' |-----------|-------------------|--------------------------------------------|
#' | AIC       | [compute_aic()]   | Akaike (1974)                              |
#' | AICc      | [compute_aicc()]  | Small-sample correction; Hurvich & Tsai (1989) |
#' | BIC       | [compute_bic()]   | Schwarz (1978)                             |
#' | EBIC      | [compute_ebic()]  | Extended BIC; Chen & Chen (2008)           |
#' | RBIC      | [compute_rbic()]  | Grouped EBIC for structured predictors     |
#'
#' **Supported model classes:**
#' * `lm`, `glm` — scalar output.
#' * `glmnet`, `ncvreg` — vector output of length `length(lambda)`, enabling
#'   IC-based selection over the full regularization path without stepwise
#'   fitting.
#'
#' Stepwise model selection for `lm`/`glm` objects is available via
#' [ic_step()], which generalises [MASS::stepAIC()] to all five criteria.
#'
#' @section Getting started:
#' ```r
#' library(ICtoolkit)
#'
#' # IC for a linear model
#' fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)
#' compute_aic(fit)
#' compute_bic(fit)
#'
#' # Stepwise selection with RBIC
#' fit_null <- lm(mpg ~ 1, data = mtcars)
#' fit_full <- lm(mpg ~ ., data = mtcars)
#' P_index  <- list(continuous = c("wt", "hp", "disp"),
#'                  categorical = c("cyl", "gear", "carb", "am", "vs"))
#' ic_step(fit_null,
#'         scope     = list(lower = fit_null, upper = fit_full),
#'         direction = "forward",
#'         criterion = "RBIC",
#'         P_index   = P_index)
#' ```
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' *IEEE Transactions on Automatic Control*, **19**(6), 716–723.
#'
#' Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for
#' model selection with large model spaces. *Biometrika*, **95**(3), 759–771.
#'
#' Hurvich, C. M., & Tsai, C.-L. (1989). Regression and time series model
#' selection in small samples. *Biometrika*, **76**(2), 297–307.
#'
#' Schwarz, G. (1978). Estimating the dimension of a model.
#' *The Annals of Statistics*, **6**(2), 461–464.
#'
#' @keywords internal
"_PACKAGE"
