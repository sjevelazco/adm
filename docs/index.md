# adm - Abundance-based species distribution models ![](reference/figures/adm_logo.png)

[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.md)
[![R-CMD-check](https://github.com/sjevelazco/adm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjevelazco/adm/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/sjevelazco/adm/graph/badge.svg?token=cKRmbNhn0A)](https://codecov.io/gh/sjevelazco/adm)
[![](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://img.shields.io/badge/DOI-10.1111%2F2041--210X.70074-blue)](https://doi.org/10.1111/2041-210X.70074)
[![Ask
DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/sjevelazco/adm)

------------------------------------------------------------------------

### Overview

This package aims to support the construction of Abundance-based species
distribution models, including data preparation, model fitting,
prediction, and model exploration. The package offers several modeling
approaches (i.e., algorithms) that users can fine-tune and customize.
Models can be predicted in geographic space and explored regarding
performance and response curves. Because modeling workflows in ***adm***
are constructed based on a combination of distinct functions and simple
outputs, ***adm*** can be easily integrated into other packages.

### Structure of adm

***adm*** functions are grouped in three categories: modeling,
post-modeling, and miscellaneous tools

[![](https://raw.githubusercontent.com/sjevelazco/adm/main/man/figures/adm.png)](https://sjevelazco.github.io/adm)

#### i) **modeling**

Functions to tune, fit, and validate models with nine different
algorithms, with a suite of possible model-specific hyperparameters

***Fit and validate models without hyperparameters tuning***

- [`fit_abund_cnn()`](https://sjevelazco.github.io/adm/reference/fit_abund_cnn.md)
  Fit and validate Convolutional Neural Network Model

- [`fit_abund_dnn()`](https://sjevelazco.github.io/adm/reference/fit_abund_dnn.md)
  Fit and validate Deep Neural Network model

- [`fit_abund_gam()`](https://sjevelazco.github.io/adm/reference/fit_abund_gam.md)
  Fit and validate Generalized Additive Models

- [`fit_abund_gbm()`](https://sjevelazco.github.io/adm/reference/fit_abund_gbm.md)
  Fit and validate Generalized Boosted Regression models

- [`fit_abund_glm()`](https://sjevelazco.github.io/adm/reference/fit_abund_glm.md)
  Fit and validate Generalized Linear Models

- [`fit_abund_net()`](https://sjevelazco.github.io/adm/reference/fit_abund_net.md)
  Fit and validate Artificial Neural Network models

- [`fit_abund_raf()`](https://sjevelazco.github.io/adm/reference/fit_abund_raf.md)
  Fit and validate Random Forests models

- [`fit_abund_svm()`](https://sjevelazco.github.io/adm/reference/fit_abund_svm.md)
  Fit and validate Support Vector Machine models

- [`fit_abund_xgb()`](https://sjevelazco.github.io/adm/reference/fit_abund_xgb.md)
  Fit and validate Extreme Gradient Boosting models

***Fit and validate models with hyperparameters tuning***

- [`tune_abund_cnn()`](https://sjevelazco.github.io/adm/reference/tune_abund_cnn.md)
  Fit and validate Convolutional Neural Network with exploration of
  hyper-parameters that optimize performance

- [`tune_abund_dnn()`](https://sjevelazco.github.io/adm/reference/tune_abund_dnn.md)
  Fit and validate Deep Neural Network model with exploration of
  hyper-parameters that optimize performance

- [`tune_abund_gam()`](https://sjevelazco.github.io/adm/reference/tune_abund_gam.md)
  Fit and validate Generalized Additive Models with exploration of
  hyper-parameters that optimize performance

- [`tune_abund_gbm()`](https://sjevelazco.github.io/adm/reference/tune_abund_gbm.md)
  Fit and validate Generalized Boosted Regression models with
  exploration of hyper-parameters that optimize performance

- [`tune_abund_glm()`](https://sjevelazco.github.io/adm/reference/tune_abund_glm.md)
  Fit and validate Generalized Linear Models with exploration of
  hyper-parameters that optimize performance

- [`tune_abund_net()`](https://sjevelazco.github.io/adm/reference/tune_abund_net.md)
  Fit and validate Shallow Neural Networks models with exploration of
  hyper-parameters that optimize performance

- [`tune_abund_raf()`](https://sjevelazco.github.io/adm/reference/tune_abund_raf.md)
  Fit and validate Random Forest models with exploration of
  hyperparameters that optimize performance

- [`tune_abund_svm()`](https://sjevelazco.github.io/adm/reference/tune_abund_svm.md)
  Fit and validate Support Vector Machine models with exploration of
  hyper-parameters that optimize performance

- [`tune_abund_xgb()`](https://sjevelazco.github.io/adm/reference/tune_abund_xgb.md)
  Fit and validate Extreme Gradient Boosting models with exploration of
  hyper-parameters that optimize performance

Modeling evaluation

- [`adm_eval()`](https://sjevelazco.github.io/adm/reference/adm_eval.md)
  Calculate different model performance metrics

#### ii) **post-modeling**

Functions to predict abundance across space and construct partial
dependence plots to explore the relationships between abundance and
environmental predictors

- [`adm_predict()`](https://sjevelazco.github.io/adm/reference/adm_predict.md)
  Spatial predictions from individual and ensemble models

- [`p_abund_bpdp()`](https://sjevelazco.github.io/adm/reference/p_abund_bpdp.md)
  Bivariate partial dependence plots for abundance-based distribution
  models

- [`p_abund_pdp()`](https://sjevelazco.github.io/adm/reference/p_abund_pdp.md)
  Partial dependent plots for abundance-based distribution models

- [`data_abund_bpdp()`](https://sjevelazco.github.io/adm/reference/data_abund_bpdp.md)
  Calculate data to construct bivariate partial dependence plots

- [`data_abund_pdp()`](https://sjevelazco.github.io/adm/reference/data_abund_pdp.md)
  Calculate data to construct partial dependence plots

#### iii) **miscellaneous tools**

Extra functions to support the modeling workflow, including data
handling, transformations, and hyperparameter selection.

- [`adm_extract()`](https://sjevelazco.github.io/adm/reference/adm_extract.md)
  Extract values from a spatial raster based on x and y coordinates

- [`adm_summarize()`](https://sjevelazco.github.io/adm/reference/adm_summarize.md)
  Merge model performance tables

- [`adm_transform()`](https://sjevelazco.github.io/adm/reference/adm_transform.md)
  Performs data transformation on a variable based on the specified
  method.

- [`balance_dataset()`](https://sjevelazco.github.io/adm/reference/balance_dataset.md)
  Balance database at a given absence-presence ratio

- [`cnn_make_samples()`](https://sjevelazco.github.io/adm/reference/cnn_make_samples.md)
  Creates sample data for Convolutional Neural Network

- [`croppin_hood()`](https://sjevelazco.github.io/adm/reference/croppin_hood.md)
  Crop rasters around a point (for Convolutional Neural Networks)

- [`family_selector()`](https://sjevelazco.github.io/adm/reference/family_selector.md)
  Select probability distributions for GAM and GLM

- [`generate_arch_list()`](https://sjevelazco.github.io/adm/reference/generate_arch_list.md)
  Generate architecture list for Deep Neural Network and Convolutional
  Neural Network

- [`generate_cnn_architecture()`](https://sjevelazco.github.io/adm/reference/generate_cnn_architecture.md)
  Generate architectures for Convolutional Neural Network

- [`generate_dnn_architecture()`](https://sjevelazco.github.io/adm/reference/generate_dnn_architecture.md)
  Generate architectures for Deep Neural Network

- [`model_selection()`](https://sjevelazco.github.io/adm/reference/model_selection.md)
  Best hyper-parameters selection

- [`res_calculate()`](https://sjevelazco.github.io/adm/reference/res_calculate.md)
  Calculate the output resolution of a layer

- [`select_arch_list()`](https://sjevelazco.github.io/adm/reference/select_arch_list.md)
  Select architectures for Convolutional Neural Network or Deep Neural
  Network

### Installation

You can install the development version of ***adm*** from
[github](https://github.com/sjevelazco/adm)

``` r

# For Windows and Mac OS operating systems
remotes::install_github("sjevelazco/adm")
```

### Package website

See the package website (<https://sjevelazco.github.io/adm/>) for
functions explanation and vignettes.

### Package citation

de Oliveira Junior, A.C., Velazco, S.J.E., 2025. adm: An R package for
constructing abundance-based species distribution models. *Methods in
Ecology and Evolution*. 16, 1404–1412.
<https://doi.org/10.1111/2041-210X.70074>
