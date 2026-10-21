# Generate architecture list for Deep Neural Network and Convolutional Neural Network

This function generates a list of architectures for either a Deep Neural
Netwokd (DNN) or a Convolutional Neural Network (CNN).

## Usage

``` r
generate_arch_list(
  type,
  number_of_features,
  number_of_outputs,
  n_layers = c(1, 2),
  n_neurons = c(7),
  sample_size = c(11, 11),
  number_of_fc_layers = 1,
  fc_layers_size = c(14),
  conv_layers_kernel = 3,
  conv_layers_stride = 1,
  conv_layers_padding = 0,
  pooling = NULL,
  batch_norm = TRUE,
  dropout = 0
)
```

## Arguments

- type:

  string. Specifies the type the network. The valid inputs "dnn" and
  "cnn".

- number_of_features:

  numeric. Value that specifies the number of features in the dataset.

- number_of_outputs:

  numeric. Value that specifies the number of outputs.

- n_layers:

  numeric. Vector that specifies the number of layers in the networks.
  Default value are 1 and 2.

- n_neurons:

  vector. Specifies the number of neurons each layer. Default 7.

- sample_size:

  vector. Specifies the size. Default c(11, 11)

- number_of_fc_layers:

  numeric. Specifies the number of fully connected layers. Default 1.

- fc_layers_size:

  vector. Specifies the size of the fully connected layers. Default 14.

- conv_layers_kernel:

  numeric. Specifies the kernel size for layers. Default 3.

- conv_layers_stride:

  numeric. Specifies the stride for the convolutional layers. Default 1.

- conv_layers_padding:

  numeric. Specifies the padding for the convolutional layers. Default
  0.

- pooling:

  numeric. Specifies 2D average pooling kernel size. Default NULL

- batch_norm:

  logical. Specifies whether batch normalization is included in the
  architecture. Default TRUE.

- dropout:

  Numeric. The probability (p) of randomly zeroing elements of the input
  tensor during training to prevent overfitting. Must be between 0 (no
  dropout) and 1 (all inputs zeroed). Default is 0 (no dropout).

## Value

A list containing:

- arch_list: a list of generated architectures.

- arch_dict: a list of architecture dictionaries.

## See also

[`generate_dnn_architecture`](https://sjevelazco.github.io/adm/reference/generate_dnn_architecture.md),
[`generate_cnn_architecture`](https://sjevelazco.github.io/adm/reference/generate_cnn_architecture.md)

## Examples

``` r
# \donttest{
require(dplyr)
# Generating architectures for DNN, using batch normalization and dropout
dnn_archs <- generate_arch_list(
  type = "dnn",
  number_of_features = 4,
  number_of_outputs = 1,
  n_layers = c(2, 3),
  n_neurons = c(8, 16),
  batch_norm = TRUE,
  dropout = 0.2
)

dnn_archs$arch_dict # Matrices describing the networks
#> $`2_layer_net`
#>         [,1] [,2] [,3] [,4]
#> layer_1    8   16    8   16
#> layer_2    8    8   16   16
#> 
#> $`3_layer_net`
#>         [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> layer_1    8   16    8   16    8   16    8   16
#> layer_2    8    8   16   16    8    8   16   16
#> layer_3    8    8    8    8   16   16   16   16
#> 
length(dnn_archs$arch_list) # number of DNN architectures generated
#> [1] 12

# Generating architectures for CNN, using batch normalization, dropout and average pooling
# Note that arguments meaning change with the context

cnn_archs <- generate_arch_list(
  type = "cnn",
  number_of_features = 4,
  number_of_outputs = 1,
  n_layers = c(2, 3), # now convolutional layers
  n_neurons = c(8, 16),
  sample_size = c(11, 11),
  number_of_fc_layers = 2, # fully connected layers
  fc_layers_size = c(16, 8),
  conv_layers_kernel = 3,
  conv_layers_stride = 1,
  conv_layers_padding = 0,
  pooling = 1,
  batch_norm = TRUE,
  dropout = 0.2
)

cnn_archs$arch_dict # Matrices describing the networks
#> $`conv2-fc2-net`
#>        [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> conv_1    8   16    8   16    8   16    8   16    8    16     8    16     8
#> conv_2    8    8   16   16    8    8   16   16    8     8    16    16     8
#> fc_1     16   16   16   16    8    8    8    8   16    16    16    16     8
#> fc_2     16   16   16   16   16   16   16   16    8     8     8     8     8
#>        [,14] [,15] [,16]
#> conv_1    16     8    16
#> conv_2     8    16    16
#> fc_1       8     8     8
#> fc_2       8     8     8
#> 
#> $`conv3-fc2-net`
#>        [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> conv_1    8   16    8   16    8   16    8   16    8    16     8    16     8
#> conv_2    8    8   16   16    8    8   16   16    8     8    16    16     8
#> conv_3    8    8    8    8   16   16   16   16    8     8     8     8    16
#> fc_1     16   16   16   16   16   16   16   16    8     8     8     8     8
#> fc_2     16   16   16   16   16   16   16   16   16    16    16    16    16
#>        [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]
#> conv_1    16     8    16     8    16     8    16     8    16     8    16     8
#> conv_2     8    16    16     8     8    16    16     8     8    16    16     8
#> conv_3    16    16    16     8     8     8     8    16    16    16    16     8
#> fc_1       8     8     8    16    16    16    16    16    16    16    16     8
#> fc_2      16    16    16     8     8     8     8     8     8     8     8     8
#>        [,26] [,27] [,28] [,29] [,30] [,31] [,32]
#> conv_1    16     8    16     8    16     8    16
#> conv_2     8    16    16     8     8    16    16
#> conv_3     8     8     8    16    16    16    16
#> fc_1       8     8     8     8     8     8     8
#> fc_2       8     8     8     8     8     8     8
#> 
length(cnn_archs$arch_list) # number of CNN architectures generated
#> [1] 48

# The list size can be easily and greatly reduced with select_arch_list

dnn_archs_redux <- dnn_archs %>% select_arch_list(
  type = c("dnn"),
  method = "percentile",
  n_samples = 1,
  min_max = TRUE # Keep the network with the minimum and maximum number of parameters
)

length(dnn_archs_redux$arch_list) # reduced number of architectures
#> [1] 12

cnn_archs_redux <- cnn_archs %>% select_arch_list(
  type = c("cnn"),
  method = "percentile",
  n_samples = 1,
  min_max = TRUE # Keep the network with the minimum and maximum number of parameters
)

length(cnn_archs_redux$arch_list) # reduced number of architectures
#> [1] 19
# }
```
