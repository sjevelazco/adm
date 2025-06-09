# adm - Abundance-based species distribution models <img src="man/figures/adm_logo.svg" align="right" style="height:168px;"/>

[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html) [![R-CMD-check](https://github.com/sjevelazco/adm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sjevelazco/adm/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/gh/sjevelazco/adm/graph/badge.svg?token=cKRmbNhn0A)](https://codecov.io/gh/sjevelazco/adm) [![](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

------------------------------------------------------------------------

### Overview

This package aims to support the construction of Abundance-based species distribution models, including data preparation, model fitting, prediction, and model exploration. The package offers several modeling approaches (i.e., algorithms) that users can fine-tune and customize. Models can be predicted in geographic space and explored regarding performance and response curves. Because modeling workflows in ***adm*** are constructed based on a combination of distinct functions and simple outputs, ***adm*** can be easily integrated into other packages.

### Structure of adm

***adm*** functions are grouped in three categories: modeling, post-modeling, and miscellaneous tools

<a href='https://sjevelazco.github.io/adm'><img src="https://raw.githubusercontent.com/sjevelazco/adm/main/man/figures/adm.png" align="centre" height="450"/></a>

#### i) **modeling**

Functions to tune, fit, and validate models with nine different algorithms, with a suite of possible model-specific hyperparameters

***Fit and validate models without hyperparameters tuning***

-   `fit_abund_cnn()` Fit and validate Convolutional Neural Network Model

-   `fit_abund_dnn()` Fit and validate Deep Neural Network model

-   `fit_abund_gam()` Fit and validate Generalized Additive Models

-   `fit_abund_gbm()` Fit and validate Generalized Boosted Regression models

-   `fit_abund_glm()` Fit and validate Generalized Linear Models

-   `fit_abund_net()` Fit and validate Artificial Neural Network models

-   `fit_abund_raf()` Fit and validate Random Forests models

-   `fit_abund_svm()` Fit and validate Support Vector Machine models

-   `fit_abund_xgb()` Fit and validate Extreme Gradient Boosting models

***Fit and validate models with hyperparameters tuning***

-   `tune_abund_cnn()` Fit and validate Convolutional Neural Network with exploration of hyper-parameters that optimize performance

-   `tune_abund_dnn()` Fit and validate Deep Neural Network model with exploration of hyper-parameters that optimize performance

-   `tune_abund_gam()` Fit and validate Generalized Additive Models with exploration of hyper-parameters that optimize performance

-   `tune_abund_gbm()` Fit and validate Generalized Boosted Regression models with exploration of hyper-parameters that optimize performance

-   `tune_abund_glm()` Fit and validate Generalized Linear Models with exploration of hyper-parameters that optimize performance

-   `tune_abund_net()` Fit and validate Shallow Neural Networks models with exploration of hyper-parameters that optimize performance

-   `tune_abund_raf()` Fit and validate Random Forest models with exploration of hyperparameters that optimize performance

-   `tune_abund_svm()` Fit and validate Support Vector Machine models with exploration of hyper-parameters that optimize performance

-   `tune_abund_xgb()` Fit and validate Extreme Gradient Boosting models with exploration of hyper-parameters that optimize performance

Modeling evaluation

-   `adm_eval()` Calculate different model performance metrics

#### ii) **post-modeling**

Functions to predict abundance across space and construct partial dependence plots to explore the relationships between abundance and environmental predictors

-   `adm_predict()` Spatial predictions from individual and ensemble models

-   `p_abund_bpdp()` Bivariate partial dependence plots for abundance-based distribution models

-   `p_abund_pdp()` Partial dependent plots for abundance-based distribution models

-   `data_abund_bpdp()` Calculate data to construct bivariate partial dependence plots

-   `data_abund_pdp()` Calculate data to construct partial dependence plots

#### iii) **miscellaneous tools**

Extra functions to support the modeling workflow, including data handling, transformations, and hyperparameter selection.

-   `adm_extract()` Extract values from a spatial raster based on x and y coordinates

-   `adm_summarize()` Merge model performance tables

-   `adm_transform()` Performs data transformation on a variable based on the specified method.

-   `balance_dataset()` Balance database at a given absence-presence ratio

-   `cnn_make_samples()` Creates sample data for Convolutional Neural Network

-   `croppin_hood()` Crop rasters around a point (for Convolutional Neural Networks)

-   `family_selector()` Select probability distributions for GAM and GLM

-   `generate_arch_list()` Generate architecture list for Deep Neural Network and Convolutional Neural Network

-   `generate_cnn_architecture()` Generate architectures for Convolutional Neural Network

-   `generate_dnn_architecture()` Generate architectures for Deep Neural Network

-   `model_selection()` Best hyper-parameters selection

-   `res_calculate()` Calculate the output resolution of a layer

-   `select_arch_list()` Select architectures for Convolutional Neural Network or Deep Neural Network

### Installation

You can install the development version of ***adm*** from [github](https://github.com/sjevelazco/adm)

``` r
# For Windows and Mac OS operating systems
remotes::install_github("sjevelazco/adm")
```

### Package website

See the package website (<https://sjevelazco.github.io/adm/>) for functions explanation and vignettes.

### Package citation

Oliveira Junior A.C., Velazco S.J.E. (2025). adm: an R package for constructing abundance-based species distribution models. Methods in Ecology and Evolution, x(x) x–xx. <https://doi.org/10.1111/2041-210X70074>
