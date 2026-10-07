# 1. Modeling workflow with adm

## Introduction

Abundance-based Distribution Models (ADM) are promising approach to
construct spatially explicit correlative models of species abundance.
The *adm* package allow users to easily construct and validate those
models, meeting researches specific needs. In this vignette, users will
learn how to support a modeling workflow with *adm*, from data
preparation to model fitting and prediction.

## Installation

\
[`require`](https://rdrr.io/r/base/library.html)`(`[`adm`](https://sjevelazco.github.io/adm/)`)`\
[`require`](https://rdrr.io/r/base/library.html)`(`[`terra`](https://rspatial.org/)`)`\
[`require`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`)`

## Preparaing data

In this tutorial, we will model abundace of *Cynophalla retusa*
(Griseb.) Cornejo & Iltis (Capparaceae), a dry biome shrub native to
northeastern Argentina, Paraguay, Bolivia, and central Brazil. As
predictor variables, we will use the first 7 PC (cumulative viarance \>
90%) of a PCA performed with 35 climatic and edaphic variables. We can
load all needed data with:

\
`# Load species abundance data`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"cretusa_data"``)`\
\
\
`# Load raster with environmental variables`\
`cretusa_predictors`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(`\
`  ``"external/cretusa_predictors.tif"``,`\
`  package ``=`` ``"adm"`\
`)`\
`cretusa_predictors`` ``<-`` ``terra``::`[`rast`](https://rspatial.github.io/terra/reference/rast.html)`(``cretusa_predictors``)`\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``cretusa_predictors``)`\
`#> [1] "PC1" "PC2" "PC3" "PC4" "PC5" "PC6" "PC7"`\
\
`# Species training area`\
`sp_train_a`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"external/cretusa_calib_area.gpkg"``, package ``=`` ``"adm"``)`\
`sp_train_a`` ``<-`` ``terra``::`[`vect`](https://rspatial.github.io/terra/reference/vect.html)`(``sp_train_a``)`

Let’s explore these data

\
`# Species data`\
`# ?cretusa_data`\
`cretusa_data`` ``# species dat`\
`#> ``# A tibble: 366 × 5`\
`#>    species           ind_ha     x     y .part`\
`#>    ``<chr>``              ``<int>`` ``<dbl>`` ``<dbl>`` ``<int>`\
`#> `` 1`` Cynophalla retusa     10 -``64.5`` -``22.7``     1`\
`#> `` 2`` Cynophalla retusa     10 -``64.1`` -``22.7``     1`\
`#> `` 3`` Cynophalla retusa     20 -``64.7`` -``23.1``     2`\
`#> `` 4`` Cynophalla retusa      0 -``62.6`` -``23.2``     3`\
`#> `` 5`` Cynophalla retusa      0 -``61.7`` -``24.5``     3`\
`#> `` 6`` Cynophalla retusa      0 -``61.6`` -``25.0``     3`\
`#> `` 7`` Cynophalla retusa      0 -``61.2`` -``24.5``     3`\
`#> `` 8`` Cynophalla retusa      0 -``64.9`` -``23.8``     2`\
`#> `` 9`` Cynophalla retusa      0 -``65.3`` -``24.4``     3`\
`#> ``10`` Cynophalla retusa      0 -``64.7`` -``24.8``     1`\
`#> ``# ℹ 356 more rows`

\
`# Environmental predictors`\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``cretusa_predictors``)`\
`#> [1] "PC1" "PC2" "PC3" "PC4" "PC5" "PC6" "PC7"`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``cretusa_predictors``)`

![](v01_modelling_workflow_files/figure-html/Environmental%20predictors-1.png)

\
`# Training area`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``sp_train_a``)`

![](v01_modelling_workflow_files/figure-html/Training%20area-1.png)

\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``cretusa_predictors``[[``1``]``]``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``sp_train_a``, add ``=`` ``TRUE``)`\
[`points`](https://rspatial.github.io/terra/reference/lines.html)`(``cretusa_data`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` ``dplyr``::`[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``x``, ``y``)``, col ``=`` ``"red"``, pch ``=`` ``20``)`

![](v01_modelling_workflow_files/figure-html/unnamed-chunk-2-1.png)

In `cretusa_data` we have 366 georeferenced points for *C. retusa*.
“ind_ha” column contains abundance measured in individuals per hectare.
“x” and “y” are decimals for longitude and latitude, respectively.
“.part” has folds of a spatial block partitioning used for
cross-validation. Now we need to extract environmental data from
predictors raster. For that, we will use *adm_extract*, and columns with
x and y coordinates will be important.

\
`species_data`` ``<-`` `[`adm_extract`](https://sjevelazco.github.io/adm/reference/adm_extract.md)`(`\
`  data ``=`` ``cretusa_data``, ``# georeferenced dataframe`\
`  x ``=`` ``"x"``, ``# spatial x coordinates`\
`  y ``=`` ``"y"``, ``# spatial y coordinates`\
`  env_layer ``=`` ``cretusa_predictors``, ``# raster with environmental variables`\
`  variables ``=`` ``NULL``, ``# return data for all layers`\
`  filter_na ``=`` ``TRUE`\
`)`\
\
`species_data`\
`#> ``# A tibble: 366 × 12`\
`#>    species    ind_ha     x     y .part    PC1    PC2    PC3   PC4    PC5     PC6`\
`#>    ``<chr>``       ``<int>`` ``<dbl>`` ``<dbl>`` ``<int>``  ``<dbl>``  ``<dbl>``  ``<dbl>`` ``<dbl>``  ``<dbl>``   ``<dbl>`\
`#> `` 1`` Cynophall…     10 -``64.5`` -``22.7``     1  1.50  -``1.28``   0.773 0.784 -``0.249`` -``0.559`` `\
`#> `` 2`` Cynophall…     10 -``64.1`` -``22.7``     1  0.750 -``1.22``   1.36  0.359 -``0.434`` -``0.742`` `\
`#> `` 3`` Cynophall…     20 -``64.7`` -``23.1``     2  1.21  -``1.82``   1.36  0.336 -``0.862``  0.087``8`\
`#> `` 4`` Cynophall…      0 -``62.6`` -``23.2``     3 -``1.71``  -``2.85``  -``1.34``  0.762 -``0.745`` -``1.02``  `\
`#> `` 5`` Cynophall…      0 -``61.7`` -``24.5``     3 -``1.65``  -``2.17``  -``1.68``  0.175 -``0.839`` -``0.637`` `\
`#> `` 6`` Cynophall…      0 -``61.6`` -``25.0``     3 -``1.29``  -``2.35``  -``1.94``  0.614 -``0.944`` -``0.805`` `\
`#> `` 7`` Cynophall…      0 -``61.2`` -``24.5``     3 -``0.904`` -``2.73``  -``2.30``  1.05  -``0.604`` -``0.830`` `\
`#> `` 8`` Cynophall…      0 -``64.9`` -``23.8``     2  0.676 -``1.57``   0.401 1.54  -``0.346`` -``0.417`` `\
`#> `` 9`` Cynophall…      0 -``65.3`` -``24.4``     3  0.460 -``0.795``  0.409 2.19  -``0.172`` -``0.249`` `\
`#> ``10`` Cynophall…      0 -``64.7`` -``24.8``     1  1.18   0.201  0.606 2.65  -``0.356`` -``0.775`` `\
`#> ``# ℹ 356 more rows`\
`#> ``# ℹ 1 more variable: PC7 <dbl>`

Notice that the new dataframe have one column for each environmental
variable (layers in raster). It is possible to extract data from
specific layers using “variables” argument in *adm_extract*. To
stabilize Deep Neural Networks training, we will transform response data
using “zscore” method.

\
`species_data`` ``<-`` `[`adm_transform`](https://sjevelazco.github.io/adm/reference/adm_transform.md)`(`\
`  data ``=`` ``species_data``,`\
`  variable ``=`` ``"ind_ha"``,`\
`  method ``=`` ``"zscore"`\
`)`\
\
`species_data`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` ``dplyr``::`[`select`](https://dplyr.tidyverse.org/reference/select.html)`(``ind_ha``, ``ind_ha_zscore``)`\
`#> ``# A tibble: 366 × 2`\
`#>    ind_ha ind_ha_zscore`\
`#>     ``<int>``         ``<dbl>`\
`#> `` 1``     10        0.075``4`\
`#> `` 2``     10        0.075``4`\
`#> `` 3``     20        0.802 `\
`#> `` 4``      0       -``0.651`` `\
`#> `` 5``      0       -``0.651`` `\
`#> `` 6``      0       -``0.651`` `\
`#> `` 7``      0       -``0.651`` `\
`#> `` 8``      0       -``0.651`` `\
`#> `` 9``      0       -``0.651`` `\
`#> ``10``      0       -``0.651`` `\
`#> ``# ℹ 356 more rows`

It creates a new column called “ind_ha_zscore”, which can be used as
response variable.

## Tuning models

With all set, we are good to proceed with ADM construction. In this
tutorial, we will fine-tune, fit and validate Deep Neural Network (DNN),
Generalized Linear Models (GLM) and Random Forest (RAF) models, using
*tune_abund\_* family functions.

### RAF

Starting with RAF, the first thing we need to do is to determine values
for hyperparameters to be tested. This could be done for all or just
part of the hyperparameters. In this example, we will set values for
all. This is done creating a grid which will guide the values
exploration, what can be easily constructed with *expand.grid* base
function:

\
`raf_grid`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  mtry ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``4``)``,`\
`  ntree ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``100``, ``300``)`\
`)`\
[`head`](https://rspatial.github.io/terra/reference/headtail.html)`(``raf_grid``)`\
`#>   mtry ntree`\
`#> 1    2   100`\
`#> 2    4   100`\
`#> 3    2   300`\
`#> 4    4   300`\
[`nrow`](https://rspatial.github.io/terra/reference/dimensions.html)`(``raf_grid``)`` ``# 4 combinations of these two hyper-parameters`\
`#> [1] 4`

For RAF, “mtry” determines the number of variables randomly sampled as
candidates at each split. We set its values to {2, 4}. “ntree” determine
the number of decision trees to grow. We set its values to {100, 300}.
The grid combines every possible pair of values, totalizing 4
combinations (number of rows in “raf_grid”). Now we can use this grid
with *tune_abund_raf* to tune and validate a RAF model:

\
`mraf`` ``<-`` `[`tune_abund_raf`](https://sjevelazco.github.io/adm/reference/tune_abund_raf.md)`(`\
`  data ``=`` ``species_data``,`\
`  response ``=`` ``"ind_ha"``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``, ``"PC2"``, ``"PC3"``, ``"PC4"``, ``"PC5"``)``,`\
`  partition ``=`` ``".part"``,`\
`  predict_part ``=`` ``TRUE``, ``# predictions for every partition will be returned`\
`  grid ``=`` ``raf_grid``,`\
`  metrics ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"corr_pear"``, ``"mae"``)``, ``# metrics to select the best model`\
`  n_cores ``=`` ``2``, ``# number of cores to be used in parallel processing`\
`  verbose ``=`` ``FALSE`\
`)`\
`#> Using provided grid.`\
`#> Searching for optimal hyperparameters...`\
`#> `\
`#> Fitting the best model...`\
`#> The best model was achieved with: `\
`#>  mtry = 2 and ntree = 100`

The function returns a list with the following elements:

\
[`names`](https://rspatial.github.io/terra/reference/names.html)`(``mraf``)`\
`#> [1] "model"               "predictors"          "performance"        `\
`#> [4] "performance_part"    "predicted_part"      "metadata"           `\
`#> [7] "optimal_combination" "all_combinations"`

“model”: a “randomForest” class object.

\
[`class`](https://rdrr.io/r/base/class.html)`(``mraf``$``model``)`\
`#> [1] "randomForest.formula" "randomForest"`\
`mraf``$``model`\
`#> `\
`#> Call:`\
`#>  randomForest(formula = formula1, data = data, mtry = mtry, ntree = ntree,      importance = FALSE) `\
`#>                Type of random forest: regression`\
`#>                      Number of trees: 100`\
`#> No. of variables tried at each split: 2`\
`#> `\
`#>           Mean of squared residuals: 150.5632`\
`#>                     % Var explained: 20.37`

“predictors”: a tibble containing relevant informations about the model
fitted.

\
`mraf``$``predictors`\
`#> ``# A tibble: 1 × 7`\
`#>   model response c1    c2    c3    c4    c5   `\
`#>   ``<chr>`` ``<chr>``    ``<chr>`` ``<chr>`` ``<chr>`` ``<chr>`` ``<chr>`\
`#> ``1`` raf   ind_ha   PC1   PC2   PC3   PC4   PC5`

“performance”: a tibble containing the best models’ performance.

\
`mraf``$``performance`\
`#> ``# A tibble: 1 × 13`\
`#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean`\
`#>   ``<chr>``    ``<dbl>``  ``<dbl>``           ``<dbl>``         ``<dbl>``          ``<dbl>`\
`#> ``1`` raf       7.84   1.83           0.384         0.200          0.316`\
`#> ``# ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,`\
`#> ``#   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>`

“performance_part”: a tibble with the performance of each partition.

\
`mraf``$``performance_part`\
`#> ``# A tibble: 3 × 9`\
`#>   replica partition model   mae corr_spear corr_pear    inter slope pdisp`\
`#>   ``<chr>``   ``<chr>``     ``<chr>`` ``<dbl>``      ``<dbl>``     ``<dbl>``    ``<dbl>`` ``<dbl>`` ``<dbl>`\
`#> ``1`` 1       1         raf    9.85      0.594     0.478 -``0.001``10`` 1.22  0.393`\
`#> ``2`` 1       2         raf    7.39      0.362     0.281  3.24    0.445 0.631`\
`#> ``3`` 1       3         raf    6.27      0.195     0.190  3.15    0.429 0.443`

“predicted_part”: predictions for each partition.

\
`mraf``$``predicted_part`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rspatial.github.io/terra/reference/headtail.html)`(``)`\
`#> ``# A tibble: 6 × 4`\
`#>   replica partition observed predicted`\
`#>   ``<chr>``   ``<chr>``        ``<int>``     ``<dbl>`\
`#> ``1`` 1       1               10     2.47 `\
`#> ``2`` 1       1               10     3.02 `\
`#> ``3`` 1       1                0     0.627`\
`#> ``4`` 1       1               10    10.4  `\
`#> ``5`` 1       1               10    10.8  `\
`#> ``6`` 1       1               10     9.13`

“optimal_combination”: the set of hyperparameters values considered the
best given the metrics and its performance.

\
`mraf``$``optimal_combination`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` ``dplyr``::`[`glimpse`](https://pillar.r-lib.org/reference/glimpse.html)`(``)`\
`#> Rows: 1`\
`#> Columns: 16`\
`#> $ comb_id         ``<chr>`` "comb_1"`\
`#> $ mtry            ``<dbl>`` 2`\
`#> $ ntree           ``<dbl>`` 100`\
`#> $ model           ``<chr>`` "raf"`\
`#> $ mae_mean        ``<dbl>`` 7.835086`\
`#> $ mae_sd          ``<dbl>`` 1.828213`\
`#> $ corr_spear_mean ``<dbl>`` 0.383678`\
`#> $ corr_spear_sd   ``<dbl>`` 0.2003848`\
`#> $ corr_pear_mean  ``<dbl>`` 0.3163061`\
`#> $ corr_pear_sd    ``<dbl>`` 0.1470123`\
`#> $ inter_mean      ``<dbl>`` 2.130081`\
`#> $ inter_sd        ``<dbl>`` 1.846114`\
`#> $ slope_mean      ``<dbl>`` 0.6965336`\
`#> $ slope_sd        ``<dbl>`` 0.4493598`\
`#> $ pdisp_mean      ``<dbl>`` 0.4891603`\
`#> $ pdisp_sd        ``<dbl>`` 0.1252488`

“all_combinations”: performance for every hyper-parameter combination.

\
`mraf``$``all_combinations`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`head`](https://rspatial.github.io/terra/reference/headtail.html)`(``)`\
`#>   comb_id mtry ntree model mae_mean   mae_sd corr_spear_mean corr_spear_sd`\
`#> 1  comb_1    2   100   raf 7.835086 1.828213       0.3836780     0.2003848`\
`#> 2  comb_2    4   100   raf 7.910588 1.988897       0.3592797     0.2033749`\
`#> 3  comb_3    2   300   raf 7.834701 1.875159       0.3757081     0.1983759`\
`#> 4  comb_4    4   300   raf 7.753951 1.980398       0.3789300     0.1930393`\
`#>   corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean  slope_sd`\
`#> 1      0.3163061    0.1470123   2.130081 1.846114  0.6965336 0.4493598`\
`#> 2      0.2707336    0.1664333   3.664373 1.007601  0.5104461 0.3514737`\
`#> 3      0.3135883    0.1440350   2.372848 1.609151  0.6766354 0.4353365`\
`#> 4      0.2970138    0.1489508   3.184253 0.830947  0.5748996 0.3176335`\
`#>   pdisp_mean  pdisp_sd`\
`#> 1  0.4891603 0.1252488`\
`#> 2  0.5496524 0.1252156`\
`#> 3  0.4988350 0.1264250`\
`#> 4  0.5307853 0.1113435`

### GLM

To tune a GLM we perform basically the same steps to tune a RAF, paying
attention to GLM singularities. First, we need to construct a grid, just
like before. However, GLM needs a “distribution” hyper-parameter that
specifies the probability distribution family to be used. Choosing a
distribution family needs attention and must be done wisely, but *adm*
provides help via *family_selector* function. This function compares the
response variable range to the *gamlss* compatible families:

\
`suitable_families`` ``<-`` `[`family_selector`](https://sjevelazco.github.io/adm/reference/family_selector.md)`(`\
`  data ``=`` ``species_data``,`\
`  response ``=`` ``"ind_ha"`\
`)`\
`#> Response variable is discrete. Both continuous and discrete families will be tested.`\
`#> Selected 61 suitable families for the data.`

The function returns a tibble with suitable families information. The
column “family_call” can be directly used in a grid.

If you are interested in exploring the attributes of the families for
GLM and GAM, you can use the *families_bank* database. For further
details about family distributions see
[`?gamlss.dist::gamlss.family`](https://rdrr.io/pkg/gamlss.dist/man/gamlss.family.html).

\
`fm`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"external/families_bank.txt"``, package ``=`` ``"adm"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``utils``::`[`read.delim`](https://rdrr.io/r/utils/read.table.html)`(``.``, header ``=`` ``TRUE``, quote ``=`` ``"\t"``)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  ``dplyr``::`[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)`(``)`\
`fm`\
`#> ``# A tibble: 87 × 9`\
`#>    family_name             family_call range no_parameters discrete accepts_zero`\
`#>    ``<chr>``                   ``<chr>``       ``<chr>``         ``<int>``    ``<int>``        ``<int>`\
`#> `` 1`` Beta                    BE          (0,1)             2        0            0`\
`#> `` 2`` Beta one inflated       BEOI        (0,1]             3        0            0`\
`#> `` 3`` Box-Cox Cole and Green  BCCG        (0, …             3        0            0`\
`#> `` 4`` Box-Cox Power Exponent… BCPE        (0, …             4        0            0`\
`#> `` 5`` Box-Cox-t               BCT         (0, …             4        0            0`\
`#> `` 6`` Exponential             EXP         (0, …             1        0            0`\
`#> `` 7`` Gamma                   GA          (0, …             2        0            0`\
`#> `` 8`` Generalized Beta type 1 GB1         (0,1)             4        0            0`\
`#> `` 9`` Generalized Beta type 2 GB2         (0, …             4        0            0`\
`#> ``10`` Generalized Gamma       GG          (0, …             3        0            0`\
`#> ``# ℹ 77 more rows`\
`#> ``# ℹ 3 more variables: one_restricted <int>, accepts_one <int>,`\
`#> ``#   accepts_negatives <int>`

In this example, we selected a few suitable distributions for use (to
keep this vignette fast; in a real analysis, explore more of them). Now,
we can construct the grid:

\
`glm_grid`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(`\
`  distribution ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"NO"``,`\
`    ``"TF"``,`\
`    ``"PIG"`\
`  ``)``,`\
`  poly ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``2``)``,`\
`  inter_order ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``1``)`\
`)`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)\
`  `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(``)`\
`` # Note that in `distribution` argument it is necessary use the acronyms of `family_call` column ``

For GLM, “poly” refers to the polynomials used in model formula, and
“inter_order” refers to the interaction order between the variables.
Tuning the GLM with the grid:

\
`mglm`` ``<-`` `[`tune_abund_glm`](https://sjevelazco.github.io/adm/reference/tune_abund_glm.md)`(`\
`  data ``=`` ``species_data``,`\
`  response ``=`` ``"ind_ha"``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``, ``"PC2"``, ``"PC3"``, ``"PC4"``, ``"PC5"``)``,`\
`  predictors_f ``=`` ``NULL``,`\
`  partition ``=`` ``".part"``,`\
`  predict_part ``=`` ``TRUE``,`\
`  grid ``=`` ``glm_grid``,`\
`  metrics ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"corr_pear"``, ``"mae"``)``,`\
`  n_cores ``=`` ``2``,`\
`  verbose ``=`` ``FALSE`\
`)`\
`#> Using provided grid.`\
`#> Searching for optimal hyperparameters...`\
`#> `\
`#> Fitting the best model...`\
`#> The best model was achieved with:`\
`#>  distribution = TF`\
`#>  poly = 2`\
`#>  inter_order = 0`

The output is a list with basically the same elements as
*tune_abund_raf*, so we don’t need to go over it again. The only
difference here is the “model”, which is now is a “gamlss” class object.

\
[`class`](https://rdrr.io/r/base/class.html)`(``mglm``$``model``)`\
`#> [1] "gamlss" "gam"    "glm"    "lm"`\
`mglm``$``model`\
`#> `\
`#> Family:  c("TF", "t Family") `\
`#> Fitting method: RS() `\
`#> `\
`#> Call:  gamlss::gamlss(formula = formula1, sigma.formula = sigma_formula,  `\
`#>     nu.formula = nu_formula, tau.formula = tau_formula,  `\
`#>     family = family, data = data, control = control_gamlss,      trace = FALSE) `\
`#> `\
`#> `\
`#> Mu Coefficients:`\
`#> (Intercept)          PC1          PC2          PC3          PC4          PC5  `\
`#>      1.1361       1.1419       3.0295       0.4442       0.2818      -2.2983  `\
`#>    I(PC1^2)     I(PC2^2)     I(PC3^2)     I(PC4^2)     I(PC5^2)  `\
`#>     -0.7504       2.2560      -0.4292      -0.5493       0.7783  `\
`#> Sigma Coefficients:`\
`#> (Intercept)  `\
`#>       1.596  `\
`#> Nu Coefficients:`\
`#> (Intercept)  `\
`#>      0.5979  `\
`#> `\
`#>  Degrees of Freedom for the fit: 13 Residual Deg. of Freedom   353 `\
`#> Global Deviance:     2644.64 `\
`#>             AIC:     2670.64 `\
`#>             SBC:     2721.37`

### DNN

For DNN, in addition to a grid, the function *tune_abund_dnn* also needs
a list of architectures to test, or a single one. We will use
*generate_arch_list* function for this purpose, . However, here things
could get complicated. As DNN architecture aspects such the number and
size of layers, batch normalization, and dropout are customizable in
*adm*, many architectures could be generated at once, with all possible
combinations between these parameters, resulting in very large lists. To
filter this list, users can use *select_arch_list* to reduce it,
sampling the list using the number of parameters as a net complexity
measurement. This is highly recommended. Let’s create and select some
architectures:

\
`archs`` ``<-`` ``adm``::`[`generate_arch_list`](https://sjevelazco.github.io/adm/reference/generate_arch_list.md)`(`\
`  type ``=`` ``"dnn"``,`\
`  number_of_features ``=`` ``5``, ``# input/predictor variables`\
`  number_of_outputs ``=`` ``1``, ``# output/response variable`\
`  n_layers ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``3``, ``4``)``, ``# possible number of layers`\
`  n_neurons ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``5``, ``14``, ``21``)``, ``# possible number of neurons on each layer`\
`  batch_norm ``=`` ``TRUE``, ``# batch normalization between layers`\
`  dropout ``=`` ``0`` ``# without training dropout`\
`)`\
\
`number_before`` ``<-`` ``archs``$``arch_list`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`length`](https://rdrr.io/r/base/length.html)`(``)`` ``# 117`\
\
`archs`` ``<-`` ``adm``::`[`select_arch_list`](https://sjevelazco.github.io/adm/reference/select_arch_list.md)`(`\
`  arch_list ``=`` ``archs``,`\
`  type ``=`` ``"dnn"``,`\
`  method ``=`` ``"percentile"``, ``# sample by number of parameters percentiles`\
`  n_samples ``=`` ``2``, ``# at least two with each number of layers`\
`  min_max ``=`` ``TRUE`` ``# keep the more simple and the more complex networks`\
`)`\
\
`archs``$``arch_list`` `[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)` `[`length`](https://rdrr.io/r/base/length.html)`(``)`\
`#> [1] 50`

However, for the sake of brevity in this tutorial, we will manually
reduce even more our architectures list to just a few:

\
`archs``$``arch_list`` ``<-`` ``archs``$``arch_list``[`[`seq`](https://rdrr.io/r/base/seq.html)`(``from ``=`` ``1``, to ``=`` `[`length`](https://rdrr.io/r/base/length.html)`(``archs``)``, by ``=`` ``10``)``]`\
\
[`length`](https://rdrr.io/r/base/length.html)`(``archs``$``arch_list``)`\
`#> [1] 1`

Now we can construct the grid with hyper-parameters combinations. Note
that hyper-parameters values will be tested for each architecture.
Therefore user needs to be careful with grid and architecture list
sizes.

\
`dnn_grid`` ``<-`` `[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(`\
`  batch_size ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``64``)``,`\
`  validation_patience ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``5``)``,`\
`  fitting_patience ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``5``)``,`\
`  learning_rate ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.005``, ``0.0005``)``,`\
`  n_epochs ``=`` ``200`\
`)`\
[`head`](https://rspatial.github.io/terra/reference/headtail.html)`(``dnn_grid``)`\
`#>   batch_size validation_patience fitting_patience learning_rate n_epochs`\
`#> 1         64                   5                5         5e-03      200`\
`#> 2         64                   5                5         5e-04      200`\
[`nrow`](https://rspatial.github.io/terra/reference/dimensions.html)`(``dnn_grid``)`\
`#> [1] 2`

Now we can use the architectures generated and the grid created within
the *tune_abund_dnn* function:

\
`mdnn`` ``<-`` `[`tune_abund_dnn`](https://sjevelazco.github.io/adm/reference/tune_abund_dnn.md)`(`\
`  data ``=`` ``species_data``,`\
`  response ``=`` ``"ind_ha_zscore"``, ``# using the transformed response`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``, ``"PC2"``, ``"PC3"``, ``"PC4"``, ``"PC5"``)``,`\
`  predictors_f ``=`` ``NULL``,`\
`  partition ``=`` ``".part"``,`\
`  predict_part ``=`` ``TRUE``,`\
`  grid ``=`` ``dnn_grid``,`\
`  architectures ``=`` ``archs``,`\
`  metrics ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"corr_pear"``, ``"mae"``)``,`\
`  n_cores ``=`` ``2``,`\
`  verbose ``=`` ``FALSE`\
`)`\
`#> Using provided architectures.`\
`#> Adding default hyperparameter for: weight_decay`\
`#> Testing 6 combinations.`\
`#> Searching for optimal hyperparameters...`\
`#> `\
`#> Fitting the best model...`\
`#> Warning: Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`#> Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`#> Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`#> Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`#> Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`#> Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`#> Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`#> The best model was achieved with: `\
`#>  learning_rate = 0.005`\
`#>  n_epochs = 200`\
`#>  patience = 5 and 5`\
`#>  batch_size = 64`\
`#>  arch = 2 layers with 5->5 neurons`

Again, the output is very similar as before, because they are
standardize for all *tune_abund\_* functions. Now, the “model” element
is a “luz_module_fitted” from *torch* and *luz* packages.

\
[`class`](https://rdrr.io/r/base/class.html)`(``mdnn``$``model``)`\
`#> [1] "luz_module_fitted"`\
`mdnn``$``model`\
`` #> A `luz_module_fitted` ``\
`#> ── Time ────────────────────────────────────────────────────────────────────────`\
`#> • Total time: 2s`\
`#> • Avg time per training epoch: 112ms`\
`#> `\
`#> ── Results ─────────────────────────────────────────────────────────────────────`\
`#> Metrics observed in the last epoch.`\
`#> `\
`#> ``ℹ`` Training:`\
`#> loss: 0.5333`\
`#> `\
`#> ── Model ───────────────────────────────────────────────────────────────────────`\
`` #> An `nn_module` containing 86 parameters. ``\
`#> `\
`#> ── Modules ─────────────────────────────────────────────────────────────────────`\
`#> • linear1: <nn_linear> #30 parameters`\
`#> • linear2: <nn_linear> #30 parameters`\
`#> • output: <nn_linear> #6 parameters`\
`#> • bn1: <nn_batch_norm1d> #10 parameters`\
`#> • bn2: <nn_batch_norm1d> #10 parameters`

### Summarizing results

In *adm* is possible to quick summarize several models evalutions in one
dataframe, using *adm_summarize* function:

\
`models_list`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``mraf``, ``mglm``)`\
`if`` ``(``has_torch``)`` ``models_list`` ``<-`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`list`](https://rdrr.io/r/base/list.html)`(``mdnn``)``, ``models_list``)`\
[`adm_summarize`](https://sjevelazco.github.io/adm/reference/adm_summarize.md)`(``models_list``)`\
`#> ``# A tibble: 3 × 14`\
`#>   model_ID model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean`\
`#>      ``<int>`` ``<chr>``    ``<dbl>``  ``<dbl>``           ``<dbl>``         ``<dbl>``          ``<dbl>`\
`#> ``1``        1 dnn      0.579  0.176           0.404         0.125          0.321`\
`#> ``2``        2 raf      7.84   1.83            0.384         0.200          0.316`\
`#> ``3``        3 glm      7.41   2.62            0.505         0.167          0.428`\
`#> ``# ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,`\
`#> ``#   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>`

## Predicting models

Once we fitted models, we can make spatial predictions of them. The
process is straightforward with the `adm_predict` function. It can make
predictions for more than one model at once. For illustration purposes,
we will make predictions for a calibration area delimited by a 100 km
buffered minimum convex polygon around species presence points, but this
process is optional.

\
`sp_train_a`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"external/cretusa_calib_area.gpkg"``, package ``=`` ``"adm"``)`\
`sp_train_a`` ``<-`` ``terra``::`[`vect`](https://rspatial.github.io/terra/reference/vect.html)`(``sp_train_a``)`\
\
`preds`` ``<-`` ``adm``::`[`adm_predict`](https://sjevelazco.github.io/adm/reference/adm_predict.md)`(`\
`  models ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``mraf``, ``mglm``)``,`\
`  pred ``=`` ``cretusa_predictors``,`\
`  predict_area ``=`` ``sp_train_a``,`\
`  training_data ``=`` ``species_data``,`\
`  transform_negative ``=`` ``TRUE`` ``# negative predictions will be considered 0`\
`)`\
`#> Predicting a list of models`

To predict DNN, we need to pay attention to some detail. As we trained
the DNN with tranformed response, if we want it to predict in the
original scale, we need to use the “invert_transform” argument in
*adm_predict*:

\
`pred_dnn`` ``<-`` `[`adm_predict`](https://sjevelazco.github.io/adm/reference/adm_predict.md)`(`\
`  models ``=`` ``mdnn``,`\
`  pred ``=`` ``cretusa_predictors``,`\
`  predict_area ``=`` ``sp_train_a``,`\
`  training_data ``=`` ``species_data``,`\
`  transform_negative ``=`` ``TRUE``,`\
`  invert_transform ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    method ``=`` ``"zscore"``,`\
`    a ``=`` `[`mean`](https://rspatial.github.io/terra/reference/summarize-generics.html)`(``species_data``$``ind_ha``)``,`\
`    b ``=`` `[`sd`](https://rdrr.io/r/stats/sd.html)`(``species_data``$``ind_ha``)`\
`  ``)`\
`)`\
`#> Predicting an individual model`\
`#> Warning: Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`

Note: transformation terms varies among tranformation methods. To learn
more about it, visit *adm_transform* documentation. Let’s visualize the
predictions:

\
[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``2`` ``+`` ``has_torch``)``)`\
`if`` ``(``has_torch``)`` `[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``pred_dnn``$``dnn``, main ``=`` ``"DNN"``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``preds``$``raf``, main ``=`` ``"RAF"``)`\
[`plot`](https://rspatial.github.io/terra/reference/plot.html)`(``preds``$``glm``, main ``=`` ``"GLM"``)`

![](v01_modelling_workflow_files/figure-html/unnamed-chunk-4-1.png)

\
[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``1``)``)`

## Exploring models

Finally, we can explore how a variable or a pair of variables affect the
predicted values. *adm* features univariate and bivariate Partial
Dependence Plots, PDP, and BPDP, respectively. To create PDP we use
`p_abund_pdp`. PDP illustrates the marginal response of one predictor.
It can provide very relevant information about residuals and model
extrapolation:

\
`# PDP for DNN model`\
`if`` ``(``has_torch``)`` ``{`\
`  ``pdp_dnn`` ``<-`` `[`p_abund_pdp`](https://sjevelazco.github.io/adm/reference/p_abund_pdp.md)`(`\
`    model ``=`` ``mdnn``, ``# the output of tune_abund_ or fit_abund_`\
`    predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``)``,`\
`    resolution ``=`` ``100``,`\
`    resid ``=`` ``TRUE``, ``# plot residuals`\
`    training_data ``=`` ``species_data``,`\
`    invert_transform ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      method ``=`` ``"zscore"``, ``# same as before`\
`      a ``=`` `[`mean`](https://rspatial.github.io/terra/reference/summarize-generics.html)`(``species_data``$``ind_ha``)``,`\
`      b ``=`` `[`sd`](https://rdrr.io/r/stats/sd.html)`(``species_data``$``ind_ha``)`\
`    ``)``,`\
`    response_name ``=`` ``"ind/ha"``, ``# this argument is for aesthetic only, and determines the name of y axis`\
`    projection_data ``=`` ``cretusa_predictors``, ``# to visualize extrapolation`\
`    rug ``=`` ``TRUE``, ``# rug plot of the predictor`\
`    colorl ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#462777"``, ``"#6DCC57"``)``, ``# projection and training values, respectively`\
`    colorp ``=`` ``"black"``, ``# residuals colors`\
`    alpha ``=`` ``0.2``,`\
`    theme ``=`` ``ggplot2``::`[`theme_classic`](https://ggplot2.tidyverse.org/reference/ggtheme.html)`(``)`` ``# a ggplot2 theme`\
`  ``)`\
`}`\
`#> Warning: Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`#> Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`# PDP for GLM model`\
`pdp_glm`` ``<-`` `[`p_abund_pdp`](https://sjevelazco.github.io/adm/reference/p_abund_pdp.md)`(`\
`  model ``=`` ``mglm``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``)``,`\
`  resolution ``=`` ``100``,`\
`  resid ``=`` ``TRUE``,`\
`  training_data ``=`` ``species_data``,`\
`  response_name ``=`` ``"ind/ha"``,`\
`  projection_data ``=`` ``cretusa_predictors``,`\
`  rug ``=`` ``TRUE``,`\
`  colorl ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#462777"``, ``"#6DCC57"``)``,`\
`  colorp ``=`` ``"black"``,`\
`  alpha ``=`` ``0.2``,`\
`  theme ``=`` ``ggplot2``::`[`theme_classic`](https://ggplot2.tidyverse.org/reference/ggtheme.html)`(``)`\
`)`\
\
`# PDP for RAF model`\
`pdp_raf`` ``<-`` `[`p_abund_pdp`](https://sjevelazco.github.io/adm/reference/p_abund_pdp.md)`(`\
`  model ``=`` ``mraf``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``)``,`\
`  resolution ``=`` ``100``,`\
`  resid ``=`` ``TRUE``,`\
`  training_data ``=`` ``species_data``,`\
`  response_name ``=`` ``"ind/ha"``,`\
`  projection_data ``=`` ``cretusa_predictors``,`\
`  rug ``=`` ``TRUE``,`\
`  colorl ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"#462777"``, ``"#6DCC57"``)``,`\
`  colorp ``=`` ``"black"``,`\
`  alpha ``=`` ``0.2``,`\
`  theme ``=`` ``ggplot2``::`[`theme_classic`](https://ggplot2.tidyverse.org/reference/ggtheme.html)`(``)`\
`)`

\
`if`` ``(``has_torch``)`` ``pdp_dnn`

![](v01_modelling_workflow_files/figure-html/all%20pdp-1.png)

\
`pdp_raf`

![](v01_modelling_workflow_files/figure-html/all%20pdp-2.png)

\
`pdp_glm`

![](v01_modelling_workflow_files/figure-html/all%20pdp-3.png)

BPDP are similar to PDP, but instead of one, it illustrates the marginal
response of a pair of variables. In this example we will use the first
and third PC. Note that for `p_abund_pdp` and `p_abund_bpdp`, any subset
of predictors can be used in the “predictors” argument. If “predictors”
argument is NULL, functions plot all variables or variables pair
combinations, respectively.

\
`# BPDP for DNN`\
`if`` ``(``has_torch``)`` ``{`\
`  ``bpdp_dnn`` ``<-`` `[`p_abund_bpdp`](https://sjevelazco.github.io/adm/reference/p_abund_bpdp.md)`(`\
`    model ``=`` ``mdnn``,`\
`    predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``, ``"PC3"``)``, ``# a pair of predictors`\
`    resolution ``=`` ``100``,`\
`    training_data ``=`` ``species_data``,`\
`    projection_data ``=`` ``cretusa_predictors``,`\
`    training_boundaries ``=`` ``"convexh"``, ``# the shape in which the training boundaries are drawn. Outside of it, its extrapolations`\
`    invert_transform ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      method ``=`` ``"zscore"``, ``# same as before`\
`      a ``=`` `[`mean`](https://rspatial.github.io/terra/reference/summarize-generics.html)`(``species_data``$``ind_ha``)``,`\
`      b ``=`` `[`sd`](https://rdrr.io/r/stats/sd.html)`(``species_data``$``ind_ha``)`\
`    ``)``,`\
`    response_name ``=`` ``"ind/ha"``,`\
`    color_gradient ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`      ``"#000004"``,`\
`      ``"#1B0A40"``,`\
`      ``"#4A0C69"``,`\
`      ``"#781B6C"``,`\
`      ``"#A42C5F"``,`\
`      ``"#CD4345"``,`\
`      ``"#EC6824"``,`\
`      ``"#FA990B"``,`\
`      ``"#F7CF3D"``,`\
`      ``"#FCFFA4"`\
`    ``)``, ``# gradient for response variable`\
`    color_training_boundaries ``=`` ``"white"``,`\
`    theme ``=`` ``ggplot2``::`[`theme_classic`](https://ggplot2.tidyverse.org/reference/ggtheme.html)`(``)`\
`  ``)`\
`}`\
`#> Warning: Some torch operators might not yet be implemented for the MPS device. A`\
`` #> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a ``\
`#> fall back for those operators:`\
`#> ``ℹ```  Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use ``\
`` #>   `usethis::edit_r_environ()`. ``\
`#> ``✖```  Using `Sys.setenv()` doesn't work because the env var must be set before R ``\
`#>   starts.`\
`# BPDP for GLM`\
`bpdp_glm`` ``<-`` `[`p_abund_bpdp`](https://sjevelazco.github.io/adm/reference/p_abund_bpdp.md)`(`\
`  model ``=`` ``mglm``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``, ``"PC3"``)``,`\
`  resolution ``=`` ``100``,`\
`  training_data ``=`` ``species_data``,`\
`  projection_data ``=`` ``cretusa_predictors``,`\
`  training_boundaries ``=`` ``"convexh"``,`\
`  response_name ``=`` ``"ind/ha"``,`\
`  color_gradient ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"#000004"``,`\
`    ``"#1B0A40"``,`\
`    ``"#4A0C69"``,`\
`    ``"#781B6C"``,`\
`    ``"#A42C5F"``,`\
`    ``"#CD4345"``,`\
`    ``"#EC6824"``,`\
`    ``"#FA990B"``,`\
`    ``"#F7CF3D"``,`\
`    ``"#FCFFA4"`\
`  ``)``,`\
`  color_training_boundaries ``=`` ``"white"``,`\
`  theme ``=`` ``ggplot2``::`[`theme_classic`](https://ggplot2.tidyverse.org/reference/ggtheme.html)`(``)`\
`)`\
\
`# BPDP for RAF`\
`bpdp_raf`` ``<-`` `[`p_abund_bpdp`](https://sjevelazco.github.io/adm/reference/p_abund_bpdp.md)`(`\
`  model ``=`` ``mraf``,`\
`  predictors ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"PC1"``, ``"PC3"``)``,`\
`  resolution ``=`` ``100``,`\
`  training_data ``=`` ``species_data``,`\
`  projection_data ``=`` ``cretusa_predictors``,`\
`  training_boundaries ``=`` ``"convexh"``,`\
`  response_name ``=`` ``"ind/ha"``,`\
`  color_gradient ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"#000004"``,`\
`    ``"#1B0A40"``,`\
`    ``"#4A0C69"``,`\
`    ``"#781B6C"``,`\
`    ``"#A42C5F"``,`\
`    ``"#CD4345"``,`\
`    ``"#EC6824"``,`\
`    ``"#FA990B"``,`\
`    ``"#F7CF3D"``,`\
`    ``"#FCFFA4"`\
`  ``)``,`\
`  color_training_boundaries ``=`` ``"white"``,`\
`  theme ``=`` ``ggplot2``::`[`theme_classic`](https://ggplot2.tidyverse.org/reference/ggtheme.html)`(``)`\
`)`

\
`if`` ``(``has_torch``)`` ``bpdp_dnn`

![](v01_modelling_workflow_files/figure-html/all%20bpdp-1.png)

\
`bpdp_raf`

![](v01_modelling_workflow_files/figure-html/all%20bpdp-2.png)

\
`bpdp_glm`

![](v01_modelling_workflow_files/figure-html/all%20bpdp-3.png)

## Conclusions

In this vignette, we explored *adm* tools to support complete ADM
workflows. The package features functions from model tuning, fitting,
and validation to model prediction and exploration.
