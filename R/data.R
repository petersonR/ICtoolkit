#' Synthetic ICU end-of-life PTSD survey data
#'
#' Synthetic data modeled after a 2021--2022 study of family members of
#' ICU patients who were isolated from loved ones during the COVID-19 pandemic.
#' Respondents were assessed for PTSD using the Impact of Events Scale-6
#' (IES-6); a score of 10 or higher indicates likely PTSD.
#'
#' @format A data frame with 279 rows and 6 variables:
#' \describe{
#'   \item{ies_total}{Integer. IES-6 total score (1--30).}
#'   \item{ptsd}{Integer. Binary indicator: 1 if `ies_total >= 10`, 0 otherwise.}
#'   \item{age}{Integer. Respondent age in years.}
#'   \item{race}{Factor with levels `White`, `Black`, `Asian`, `Other`.}
#'   \item{hispanic}{Factor with levels `No`, `Yes`.}
#'   \item{female}{Factor with levels `No`, `Yes`.}
#' }
#'
#' @source Synthetic data generated to resemble a real ICU family-member survey.
#'   Not from actual patient records.
"covid_eol"
