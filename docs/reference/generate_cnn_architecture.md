# Generate architectures for Convolutional Neural Network

Generate architectures for Convolutional Neural Network

## Usage

``` r
generate_cnn_architecture(
  number_of_features = 7,
  number_of_outputs = 1,
  sample_size = c(11, 11),
  number_of_conv_layers = 2,
  conv_layers_size = c(14, 28),
  conv_layers_kernel = 3,
  conv_layers_stride = 1,
  conv_layers_padding = 0,
  number_of_fc_layers = 1,
  fc_layers_size = c(28),
  pooling = NULL,
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

- sample_size:

  vector. Specifies the size. Default c(11, 11)

- number_of_conv_layers:

  numeric. Specifies the number of convolutional layers. Default 2.

- conv_layers_size:

  numeric. the size the convolutional layers. Default c(14, 28).

- conv_layers_kernel:

  numeric. Specifies the kernel size for layers. Default 3.

- conv_layers_stride:

  numeric. Specifies the stride for the convolutional layers. Default 1.

- conv_layers_padding:

  numeric. Specifies the padding for the convolutional layers. Default
  0.

- number_of_fc_layers:

  numeric. Specifies the number of fully connected layers. Default 1.

- fc_layers_size:

  vector. Specifies the size of the fully connected layers. Default 14.

- pooling:

  numeric. Specifies 2D average pooling kernel size. Default NULL

- batch_norm:

  logical. Specifies whether batch normalization is included in the
  architecture. Default TRUE.

- dropout:

  Numeric. The probability (p) of randomly zeroing elements of the input
  tensor during training to prevent overfitting. Must be between 0 (no
  dropout) and 1 (all inputs zeroed). Default is 0 (no dropout).

- verbose:

  logical. Specifies whether the architecture should be printed. Default
  FALSE.

## Value

A list containing:

- net: a instantiated torch neural net.

- arch: a string with a R expression to instantiate the neural network.

- arch_dict: a list with a matrix describing the architecture structure.

## See also

[`generate_arch_list`](https://sjevelazco.github.io/adm/reference/generate_arch_list.md),
[`generate_dnn_architecture`](https://sjevelazco.github.io/adm/reference/generate_dnn_architecture.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate a Conclutional Neural Network with:
cnn_arch <- generate_cnn_architecture(
  number_of_features = 7, # seven input variables
  number_of_outputs = 1, # one output
  sample_size = c(11, 11), # image dimensions
  number_of_conv_layers = 2, # two convolutional layers between input and output
  conv_layers_size = c(14, 28), # of this size, respectively
  conv_layers_kernel = 3, # with a 3 pixels kernel
  conv_layers_stride = 1, # walking 1 pixel at a time
  conv_layers_padding = 1, # with 1 pixel of padding
  number_of_fc_layers = 1, # followed by one fully connected layer
  fc_layers_size = c(28), # with 28 neurons
  pooling = NULL, # without average pooling
  batch_norm = TRUE, # with batch normalization
  dropout = 0 # and without dropout
)

cnn_arch$net() # a torch net
cnn_arch$arch %>% cat() # the torch code to create it
cnn_arch$arch_dict # and a quick description of its structure
} # }
```
