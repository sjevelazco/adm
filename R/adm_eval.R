#' Calculate different model performance metrics
#'
#' @description Calculate different model performance related to model accuracy, discrimination, and precision.
#'
#' @param obs numeric. Observed abundance
#' @param pred numeric. Predicted abundance
#'
#' @importFrom stats cor lm sd
#' @importFrom dplyr as_tibble
#'
#' @return a tibble with next columns: corr_spear, corr_pear, mae, inter, slope, and pdisp(see details)
#'
#' @export
#'
#' @details
#' This function calculate metric related to:
#' \itemize{
#'   \item Accuracy: mean absolute error (mae)
#'   TODO explain inter = inter slope = slope,
#'   \item Discrimination: Spearman’s rank correlation (corr_spear) and Pearson’s correlation (corr_pear)
#'   \item Precision:  Ratio between SD of predicted and observed abundance (pdisp),
#'   }
#'   Further details see Waldock et al. (2022)
#'
#' @references
#' \itemize{
#'   \item Waldock, C., Stuart-Smith, R.D., Albouy, C., Cheung, W.W.L., Edgar, G.J., Mouillot, D., Tjiputra, J., Pellissier, L., 2022. A quantitative review of abundance-based species distribution models. Ecography https://doi.org/10.1111/ecog.05694 }
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' TODO
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
    corr_spear = corr_spear,
    corr_pear = corr_pear,
    mae = mae,
    inter = inter,
    slope = slope,
    pdisp = pdispersion
  ) %>% dplyr::as_tibble()

  return(result)
}
