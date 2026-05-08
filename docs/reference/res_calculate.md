# Calculate the output resolution of a layer

Calculate the output resolution of a layer or pooling operation in a
Convolutional Neural Network.

## Usage

``` r
res_calculate(
  type = c("layer", "pooling"),
  in_res,
  kernel_size,
  stride,
  padding
)
```

## Arguments

- type:

  string. Accepted values are "layer" and "pooling".

- in_res:

  integer. It represents the resolution of the input layer.

- kernel_size:

  integer. It refers to the size of the kernel used in the convolution
  or pooling operation.

- stride:

  integer. It is the stride length for the convolution or pooling
  operation. Only used when type is "layer".

- padding:

  integer. It is the amount padding added to the input layer. Only used
  when type is "layer"

## Value

The function returns integer which is the output resolution.

## Details

- When type is "layer" the output resolution is calculated as ((in_res -
  kernel_size + (2 \* padding)) / stride) + 1.

- When type is "pooling", the output resolution is calculated as the
  floor division of in_res by kernel_size.

## Examples

``` r
if (FALSE) { # \dontrun{

# Calculating output resolution for a convolution layer
res_calculate(type = "layer", in_res = 12, kernel_size = 2, stride = 2, padding = 0)

# Calculating output resolution for a pooling layer
res_calculate(type = "pooling", in_res = 12, kernel_size = 2)
} # }
```
