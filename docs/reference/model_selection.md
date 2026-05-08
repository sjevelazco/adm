# Best hyperparameter selection

Best hyperparameter selection

## Usage

``` r
model_selection(hyper_combinations, metrics)
```

## Arguments

- hyper_combinations:

  tibble or data.frame. with hyperparameter combinations and performance
  metrics

- metrics:

  character. with the performance metrics to be considered

## Value

A list containing:

- optimal_combination: tibble with the optimal hyperparameter
  combination.

- all_combinations: all the other combinations.
