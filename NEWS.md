# ICtoolkit 0.1.0

## Initial release

* `compute_aic()`, `compute_aicc()`, `compute_bic()`: Standard information
  criteria with S3 dispatch for `lm`, `glm`, `glmnet`, and `ncvreg` objects.
  For `glmnet`/`ncvreg` the full regularization path is evaluated, returning
  a vector of length `length(lambda)`.

* `compute_ebic()`: Extended BIC (Chen & Chen 2008) with a data-adaptive
  default gamma and support for user-supplied `gammafn`.

* `compute_rbic()`: Grouped EBIC for structured predictors, applying
  group-specific penalties defined by a `P_index` list.  Designed for
  multimodal or hierarchically organized predictor sets.

* `ic_step()`: Stepwise model selection for `lm`/`glm` objects supporting
  all five criteria, forward/backward/both directions, and a `step_path`
  attribute recording the selection trajectory.

* All `compute_*` functions return numeric vectors/scalars with informative
  attributes: `fit_class`, `k`, `criterion`, and (where applicable) `gamma`,
  `P`, `P_index`, `lambda`.
