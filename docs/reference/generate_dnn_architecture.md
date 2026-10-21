# Generate architectures for Deep Neural Network

Generate architectures for Deep Neural Network

## Usage

``` r
generate_dnn_architecture(
  number_of_features = 7,
  number_of_outputs = 1,
  number_of_hidden_layers = 2,
  hidden_layers_size = c(14, 7),
  batch_norm = TRUE,
  dropout = 0,
  verbose = FALSE
)
```

## Arguments

- number_of_features:

  numeric. Value that specifies the number of features in the dataset.

- number_of_outputs:

  numeric. Value that specifies the number of outputs.

- number_of_hidden_layers:

  numeric. Number of hidden layers in the neural network. Default 2.

- hidden_layers_size:

  numeric vector. Size of each hidden layer in the neural network.
  Default c(14, 7).

- batch_norm:

  logical. Whether to include batch normalization layers. Default TRUE.

- dropout:

  logical. Specifies whether dropout is included in the architecture.
  Default FALSE.

- verbose:

  logical. Whether to print the architecture. Default FALSE.

## Value

A list containing:

- net: a instantiated torch neural net.

- arch: a string with a R expression to instantiate the neural network.

- arch_dict: a list with a matrix describing the architecture structure.

## See also

[`generate_arch_list`](https://sjevelazco.github.io/adm/reference/generate_arch_list.md),
[`generate_cnn_architecture`](https://sjevelazco.github.io/adm/reference/generate_cnn_architecture.md)

## Examples

``` r
# \donttest{
require(dplyr)
# Generate a Deep Neural Network with:
dnn_arch <- generate_dnn_architecture(
  number_of_features = 8, # eight input variables
  number_of_outputs = 1, # one output
  number_of_hidden_layers = 5, # five layers between input and output
  hidden_layers_size = c(8, 16, 32, 16, 8), # of this size, respectively
  batch_norm = TRUE, # with batch normalization
  dropout = 0, # without dropout
)

dnn_arch$net() # a torch net
#> An `nn_module` containing 1,593 parameters.
#> 
#> ── Modules ─────────────────────────────────────────────────────────────────────
#> • linear1: <nn_linear> #72 parameters
#> • linear2: <nn_linear> #144 parameters
#> • linear3: <nn_linear> #544 parameters
#> • linear4: <nn_linear> #528 parameters
#> • linear5: <nn_linear> #136 parameters
#> • output: <nn_linear> #9 parameters
#> • bn1: <nn_batch_norm1d> #16 parameters
#> • bn2: <nn_batch_norm1d> #32 parameters
#> • bn3: <nn_batch_norm1d> #64 parameters
#> • bn4: <nn_batch_norm1d> #32 parameters
#> • bn5: <nn_batch_norm1d> #16 parameters
dnn_arch$arch %>% cat() # the torch code to create it
#> net <- torch::nn_module(
#>   'neural_net',
#>   initialize = function() {
#>     self$linear1<-torch::nn_linear(8,8)
#>     self$linear2<-torch::nn_linear(8,16)
#>     self$linear3<-torch::nn_linear(16,32)
#>     self$linear4<-torch::nn_linear(32,16)
#>     self$linear5<-torch::nn_linear(16,8)
#>     self$output<-torch::nn_linear(8,1)
#>     self$bn1<-torch::nn_batch_norm1d(num_features=8)
#>     self$bn2<-torch::nn_batch_norm1d(num_features=16)
#>     self$bn3<-torch::nn_batch_norm1d(num_features=32)
#>     self$bn4<-torch::nn_batch_norm1d(num_features=16)
#>     self$bn5<-torch::nn_batch_norm1d(num_features=8)
#>   },
#>   forward = function(x) {
#>     x %>%
#>       self$linear1() %>%
#>       torch::nnf_relu() %>%
#>       self$bn1() %>%
#>       self$linear2() %>%
#>       torch::nnf_relu() %>%
#>       self$bn2() %>%
#>       self$linear3() %>%
#>       torch::nnf_relu() %>%
#>       self$bn3() %>%
#>       self$linear4() %>%
#>       torch::nnf_relu() %>%
#>       self$bn4() %>%
#>       self$linear5() %>%
#>       torch::nnf_relu() %>%
#>       self$bn5() %>%
#>       self$output()
#>   }
#> )
dnn_arch$arch_dict # and a quick description of its structure
#> $`5_layer_net`
#>         [,1]
#> layer_1    8
#> layer_2   16
#> layer_3   32
#> layer_4   16
#> layer_5    8
#> 
# }
```
