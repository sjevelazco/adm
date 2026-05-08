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
if (FALSE) { # \dontrun{
data(sppabund)

family_selector(data = sppabund, response = "ind_ha")
} # }
```
