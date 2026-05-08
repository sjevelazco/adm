# Fit and validate Generalized Additive Models

Fit and validate Generalized Additive Models

## Usage

``` r
fit_abund_gam(
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
  inter = "automatic",
  verbose = TRUE,
  control_gamlss = gamlss::gamlss.control(trace = FALSE)
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

  logical. Save predicted abundance for testing data. Default = FALSE

- distribution:

  character. A string specifying the distribution to be used. See
  [gamlss.family](https://rdrr.io/pkg/gamlss.dist/man/gamlss.family.html)
  documentation for details. Use distribution = gamlss.dist::NO().
  Default NULL

- inter:

  integer. Number of knots in x-axis. Default "automatic"

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

- control_gamlss:

  function. control parameters of the outer iterations algorithm in
  gamlss See
  [gamlss.control](https://rdrr.io/pkg/gamlss/man/gamlss.control.html)
  documentation for details. Default gamlss.control()

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

# Fit a GAM model
mgam <- fit_abund_gam(
  data = some_sp,
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  sigma_formula = ~ elevation + bio3 + bio12,
  predictors_f = NULL,
  partition = ".part",
  distribution = gamlss.dist::ZAGA()
)

mgam
} # }
```
