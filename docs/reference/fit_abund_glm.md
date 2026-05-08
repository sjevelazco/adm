# Fit and validate Generalized Linear Models

Fit and validate Generalized Linear Models

## Usage

``` r
fit_abund_glm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  sigma_formula = ~1,
  nu_formula = ~1,
  tau_formula = ~1,
  partition = NULL,
  hold_out_set = NULL,
  predict_part = FALSE,
  distribution = NULL,
  poly = 0,
  inter_order = 0,
  control_gamlss = gamlss::gamlss.control(trace = FALSE),
  verbose = TRUE
)
```

## Arguments

- data:

  tibble or data.frame. Database with response, predictors, and
  partition values

- response:

  character. Column name with species abundance.

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). Usage predictors = c("temp",
  "precipt", "sand")

- predictors_f:

  character. Vector with the column names of qualitative predictor
  variables (i.e. ordinal or nominal variables type). Usage predictors_f
  = c("landform")

- fit_formula:

  formula. A formula object with response and predictor variables (e.g.
  formula(abund ~ temp + precipt + sand + landform)). Note that the
  variables used here must be consistent with those used in response,
  predictors, and predictors_f arguments. Default NULL

- sigma_formula:

  formula. formula for fitting a model to the nu parameter. Usage
  sigma_formula = ~ precipt + temp

- nu_formula:

  formula. formula for fitting a model to the nu parameter. Usage
  nu_formula = ~ precipt + temp

- tau_formula:

  formula. formula for fitting a model to the tau parameter. Usage
  tau_formula = ~ precipt + temp

- partition:

  character. Column name with training and validation partition groups.

- predict_part:

  logical. Save predicted abundance for testing data. Default is FALSE.

- distribution:

  character. A string specifying the distribution to be used. See
  [gamlss.family](https://rdrr.io/pkg/gamlss.dist/man/gamlss.family.html)
  documentation for details. Use distribution = gamlss.dist::NO().
  Default NULL

- poly:

  integer \>= 2. If used with values \>= 2 model will use polynomials
  for those continuous variables (i.e. used in predictors argument).
  Default is 0.

- inter_order:

  integer \>= 0. The interaction order between explanatory variables.
  Default is 0.

- control_gamlss:

  function. control parameters of the outer iterations algorithm in
  gamlss See
  [gamlss.control](https://rdrr.io/pkg/gamlss/man/gamlss.control.html)
  documentation for details. Default gamlss.control()

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "gamlss" class object from gamlss package. This object can be
  used for predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Averaged performance metrics (see
  [`adm_eval`](https://sjevelazco.github.io/adm/reference/adm_eval.md)).

- performance_part: Performance metrics for each replica and partition.

- predicted_part: Observed and predicted abundance for each test
  partition.

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)
require(dplyr)
require(gamlss)

# Database with species abundance and x and y coordinates
data("sppabund")

# Extract data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

# Explore reponse variables
some_sp$ind_ha %>% range()
some_sp$ind_ha %>% hist()

# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# Explore different family distributions
family_selector(data = some_sp, response = "ind_ha") %>% tail()

# Fit a GLM model
glm_1 <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part",
  distribution = "ZAGA",
  poly = 0,
  inter_order = 0,
  predict_part = TRUE
)

glm_1

# Using second order polynomials and first order interaction terms
glm_2 <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part",
  distribution = "ZAGA",
  poly = 2,
  inter_order = 1,
  predict_part = TRUE
)

glm_2

# Using third order polynomials and second order interaction terms
glm_3 <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part",
  distribution = "ZAGA",
  poly = 3,
  inter_order = 2,
  predict_part = TRUE
)

glm_3

# Setting formulas for different distribution parameters
glm_4 <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part",
  distribution = "ZAGA",
  fit_formula = ind_ha ~ bio12 + elevation + sand + eco,
  sigma_formula = ind_ha ~ bio12 + elevation + sand,
  poly = 0,
  inter_order = 0,
  predict_part = TRUE
)

glm_4
} # }
```
