# Select probability distributions for GAM and GLM

Select probability distribution available in gamlss.dist suited for a
given response variables (e.g., count, zero-inflated) used to fit GAM
and GLM models. See
[gamlss.family](https://rdrr.io/pkg/gamlss.dist/man/gamlss.family.html)
for more details.

## Usage

``` r
family_selector(data, response)
```

## Arguments

- data:

  data.frame or tibble. Database with species abundance

- response:

  character. Column name with species abundance

## Value

tibble with family_name, family_call, range, and discrete columns (if
family distribution or not discrete)

## Examples

``` r
# \donttest{
data(sppabund)

family_selector(data = sppabund, response = "ind_ha")
#> Response variable is continuous and need to be transformed to integer in order to test for discrete families.
#> Selected 29 suitable families for the data.
#> # A tibble: 29 × 4
#>    family_name                         family_call range   discrete
#>    <chr>                               <chr>       <chr>      <int>
#>  1 Exponential Gaussian                exGAUS      (-∞, ∞)        0
#>  2 Exponential generalized Beta type 2 EGB2        (-∞, ∞)        0
#>  3 Generalized t                       GT          (-∞, ∞)        0
#>  4 Gumbel                              GU          (-∞, ∞)        0
#>  5 Johnson's SU                        JSU         (-∞, ∞)        0
#>  6 Logistic                            LO          (-∞, ∞)        0
#>  7 Normal                              NO          (-∞, ∞)        0
#>  8 Normal Exponential t                NET         (-∞, ∞)        0
#>  9 Normal Family                       NOF         (-∞, ∞)        0
#> 10 Normal Linear Quadratic             LQNO        (-∞, ∞)        0
#> # ℹ 19 more rows
# }
```
