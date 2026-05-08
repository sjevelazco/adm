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
if (FALSE) { # \dontrun{
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
dnn_arch$arch %>% cat() # the torch code to create it
dnn_arch$arch_dict # and a quick description of its structure
} # }
```
