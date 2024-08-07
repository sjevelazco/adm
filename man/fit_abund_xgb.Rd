% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_abund_xgb.R
\name{fit_abund_xgb}
\alias{fit_abund_xgb}
\title{Fit and validate Extreme Gradient Boosting models}
\usage{
fit_abund_xgb(
  data,
  response,
  predictors,
  predictors_f,
  fit_formula = NULL,
  partition,
  predict_part = FALSE,
  nrounds = 100,
  max_depth = 5,
  eta = 0.1,
  gamma = 1,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 0.5,
  objective = "reg:squarederror",
  verbose = TRUE
)
}
\arguments{
\item{data}{tibble or data.frame. Database with response, predictors, and partition values}

\item{response}{character. Column name with species abundance.}

\item{predictors}{character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")}

\item{predictors_f}{character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")}

\item{fit_formula}{formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL}

\item{partition}{character. Column name with training and validation partition groups.}

\item{predict_part}{logical. Save predicted abundance for testing data. Default = FALSE.}

\item{nrounds}{integer. Max number of boosting iterations. Default is 100.}

\item{max_depth}{integer. The maximum depth of each tree. Default 5}

\item{eta}{numeric. The learning rate of the algorithm. Default 0.1}

\item{gamma}{numeric. Minimum loss reduction required to make a further partition on a leaf node of the tree. Default is 1.}

\item{colsample_bytree}{numeric. Subsample ratio of columns when constructing each tree. Default is 1.}

\item{min_child_weight}{numeric. Minimum sum of instance weight needed in a child. Default is 1.}

\item{subsample}{numeric. Subsample ratio of the training instance. Default is 0.5.}

\item{objective}{character. The learning task and the corresponding learning objective. Default is "reg:squarederror", regression with squared loss.}

\item{verbose}{logical. If FALSE, disables all console messages. Default TRUE.}
}
\value{
A list object with:
\itemize{
\item model: A "xgb.Booster" class object from xgboost package. This object can be used for predicting.
\item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
\item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
\item performance_part: Performance metrics for each partition.
\item predicted_part: Observed and predicted abundance for each test partition.
}
}
\description{
Fit and validate Extreme Gradient Boosting models
}
