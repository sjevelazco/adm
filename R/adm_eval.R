#' Calculate different model performance metrics
#'
#' @description This function, when supplied with observed and predicted values, calculates
#'  accuracy, discrimination, and precision between the two and returns their values as a tibble
#' table. The accuracy is evaluated through mean absolute error. The discrimination is calculated
#' using Spearman correlation, Pearson correlation, intercept and slope of a linear regression
#' between observed and predicted values. The precision is obtained from the standard deviations
#' of predicted and observed values.
#'
#' @param obs numeric. Observed abundance
#' @param pred numeric. Predicted abundance
#'
#' @importFrom dplyr as_tibble
#' @importFrom stats cor lm sd
#'
#' @return a tibble with next columns: corr_spear, corr_pear, mae, inter, slope, and pdisp(see details)
#'
#' @export
#'
#' @details
#' This function calculate metric related to the accuracy, discrimination, and precision of a model:
#' \itemize{
#'   \item Accuracy: mean absolute error (mae)
#'   \item Discrimination: Spearman<U+2019>s rank correlation (corr_spear)
#'   \item Discrimination: Pearson<U+2019>s correlation (corr_pear)
#'   \item Discrimination: regression intercept between observed and predicted values (inter)
#'   \item Discrimination: regression slope between observed and predicted values (slope)
#'   \item Precision: ratio between predicted and observed standard deviation (pdisp)
#'   }
#'   Further details see Waldock et al. (2022)
#'
#' @references
#' \itemize{
#'   \item Waldock, C., Stuart-Smith, R.D., Albouy, C., Cheung, W.W.L., Edgar, G.J., Mouillot, D., Tjiputra, J., Pellissier, L., 2022. A quantitative review of abundance-based species distribution models. Ecography https://doi.org/10.1111/ecog.05694 }
#'
#' @examples
#' \dontrun{
#' pred_a <- c(
#'   3, 2, 0, 0, 2, 5, 1, 3, 1, 2, 1, 1, 2, 5, 4,
#'   1, 2, 5, 3, 3, 4, 3, 2, 0, 2, 1, 2, 2, 1, 4,
#'   4, 2, 2, 1, 6, 1, 1, 3, 5, 0, 1, 1, 0, 1, 2
#' )
#' obs_a <- c(
#'   3, 1, 1, 3, 2, 3, 0, 3, 5, 3, 4, 2, 0, 5, 2,
#'   1, 2, 2, 3, 6, 3, 2, 4, 2, 1, 2, 3, 5, 0, 3,
#'   3, 2, 1, 2, 3, 2, 2, 1, 2, 3, 3, 1, 2, 1, 4
#' )
#'
#' adm_eval(obs = obs_a, pred = pred_a)
#' }
adm_eval <- function(obs, pred) {
  # Discrimination
  corr_spear <- stats::cor(obs, pred, method = "spearman")
  corr_pear <- stats::cor(obs, pred, method = "pearson")

  # Accuracy
  mae <- mean(abs(obs - pred))
  lm_out <- stats::lm(obs ~ pred, data = data.frame(obs, pred))
  inter <- lm_out$coefficients[1]
  slope <- lm_out$coefficients[2]

  # Precision
  pdisp <- stats::sd(pred) / stats::sd(obs)

  result <- data.frame(
    mae = mae,
    corr_spear = corr_spear,
    corr_pear = corr_pear,
    inter = inter,
    slope = slope,
    pdisp = pdisp
  ) %>% dplyr::as_tibble()

  return(result)
}
