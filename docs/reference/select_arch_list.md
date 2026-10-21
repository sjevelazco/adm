# Select architectures for Convolutional Neural Network or Deep Neural Network

Select architectures for Convolutional Neural Network or Deep Neural
Network

## Usage

``` r
select_arch_list(
  arch_list,
  type = c("dnn", "cnn"),
  method = "percentile",
  n_samples = 1,
  min_max = TRUE
)
```

## Arguments

- arch_list:

  list. Containing Convolutional Neural Network or Deep Neural Network
  architectures.

- type:

  character. Indicating the type of network. Options are "dnn" or "cnn".

- method:

  character. Indicating the method to select architectures. Default is
  "percentile".

- n_samples:

  integer. Specifying the number of samples to select per group. Default
  is 1.

- min_max:

  logical. If TRUE, include networks with minimal and maximal
  parameters.

## Value

A list with:

- arch_list: a list containing torch neural networks

- arch_dict: a list of matrices describing the structure of those
  networks

- changes: a tibble with information about neural networks name changes,
  number of parameters and deepness

## Examples

``` r
# \donttest{
require(dplyr)
# Generate some big list of architectures combining all argument values
big_arch_list <- generate_arch_list(
  type = "dnn",
  number_of_features = 4,
  number_of_outputs = 1,
  n_layers = seq(from = 2, to = 4, by = 1),
  n_neurons = c(8, 16, 32)
)

length(big_arch_list$arch_list) # number of architectures generated
#> [1] 117

# It can be reduced sampling network architectures by its parameters number

reduced_arch_list <- big_arch_list %>% select_arch_list(
  type = c("dnn"),
  method = "percentile",
  n_samples = 1, # Keep at least one of each deepness
  min_max = TRUE # Keep the network with the minimum and maximum number of parameters
)

length(reduced_arch_list$arch_list) # reduced number of architectures
#> [1] 28

# See architectures names, deepness and number of parameters
reduced_arch_list$changes
#>          class  old_name  new_name parameters deep
#> 1  2_layer_net  arch-2-1  arch-2-1        153    2
#> 2  2_layer_net  arch-2-4  arch-2-2        249    2
#> 3  2_layer_net  arch-2-3  arch-2-3        513    2
#> 4  2_layer_net  arch-2-8  arch-2-4        753    2
#> 5  2_layer_net  arch-2-6  arch-2-5        801    2
#> 6  2_layer_net  arch-2-9  arch-2-6       1377    2
#> 7  3_layer_net  arch-3-1  arch-3-1        241    3
#> 8  3_layer_net arch-3-19  arch-3-2        529    3
#> 9  3_layer_net arch-3-20  arch-3-3        649    3
#> 10 3_layer_net arch-3-22  arch-3-4        873    3
#> 11 3_layer_net  arch-3-6  arch-3-5        945    3
#> 12 3_layer_net arch-3-15  arch-3-6       1105    3
#> 13 3_layer_net arch-3-17  arch-3-7       1297    3
#> 14 3_layer_net arch-3-24  arch-3-8       1425    3
#> 15 3_layer_net arch-3-26  arch-3-9       1873    3
#> 16 3_layer_net arch-3-27 arch-3-10       2497    3
#> 17 4_layer_net  arch-4-1  arch-4-1        329    4
#> 18 4_layer_net  arch-4-4  arch-4-2        481    4
#> 19 4_layer_net arch-4-10  arch-4-3        481    4
#> 20 4_layer_net arch-4-37  arch-4-4        641    4
#> 21 4_layer_net arch-4-58  arch-4-5        769    4
#> 22 4_layer_net arch-4-20  arch-4-6        905    4
#> 23 4_layer_net arch-4-65  arch-4-7       1081    4
#> 24 4_layer_net arch-4-15  arch-4-8       1249    4
#> 25 4_layer_net arch-4-48  arch-4-9       1433    4
#> 26 4_layer_net arch-4-74 arch-4-10       1769    4
#> 27 4_layer_net arch-4-52 arch-4-11       2105    4
#> 28 4_layer_net arch-4-81 arch-4-12       3617    4
# }
```
