#' Calculate different model performance metrics
#'
#' @param obs numeric. Observed abundance
#' @param pred umeric. Predicted abundance
#'
#' @importFrom stats cor lm sd
#' @importFrom tibble as_tibble
#'
#' @return a tibble with next columns
#' @export
#'
#' @examples
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
  pdispersion <- stats::sd(pred) / stats::sd(obs)

  result <- data.frame(
    corr_spear = corr_spear,
    corr_pear = corr_pear,
    mae = mae,
    inter = inter,
    slope = slope,
    pdispersion = pdispersion
  ) %>% tibble::as_tibble()

  return(result)
}
