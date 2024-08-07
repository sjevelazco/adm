#' Generate Convolutional Neural Network architecture
#'
#' @param number_of_features
#' @param number_of_outputs
#' @param sample_size
#' @param number_of_conv_layers
#' @param conv_layers_size
#' @param conv_layers_kernel
#' @param conv_layers_stride
#' @param conv_layers_padding
#' @param number_of_fc_layers
#' @param fc_layers_size
#' @param verbose
#'
#' @importFrom torch nn_module nn_conv2d nn_linear nnf_relu torch_flatten torch_manual_seed
#'
#' @return
#' @export
#'
#' @examples
generate_cnn_architecture <-
  function(number_of_features = 7,
           number_of_outputs = 1,
           sample_size = c(11, 11),
           number_of_conv_layers = 2,
           conv_layers_size = c(14, 28),
           conv_layers_kernel = 3,
           conv_layers_stride = 1,
           conv_layers_padding = 0,
           number_of_fc_layers = 1,
           fc_layers_size = c(28),
           pooling = FALSE,
           batch_norm = TRUE,
           dropout = FALSE,
           verbose = FALSE) {
    ##
    if (any(sample_size < (conv_layers_kernel + conv_layers_padding))) {
      stop("Sample dimension is too small for the choosen configuration.")
    }

    if (pooling & !is.numeric(pooling)) {
      stop("Pooling should be a numeric value.")
    }

    ##
    arch <- "net <- torch::nn_module(
  'conv_neural_net',
  initialize = function() layer_definition,
  forward = function(x) foward_definition
)"
    conv_layer_names <- c()
    conv_layers_definition <- rep("\n    self$", number_of_conv_layers, sep = "")
    for (i in 1:length(conv_layers_definition)) {
      layer <- conv_layers_definition[[i]] %>% paste("conv", i, sep = "")
      stripped <- strsplit(layer, split = "")[[1]]
      layer_name <- stripped[6:length(stripped)] %>% paste0(collapse = "")
      conv_layer_names <- append(conv_layer_names, layer_name)

      if (i == 1) {
        layer <- paste(layer, "<-torch::nn_conv2d(", number_of_features, ",", conv_layers_size[[1]], ",", conv_layers_kernel, ",", conv_layers_stride, ",", conv_layers_padding, ")", sep = "")
      } else if (i == number_of_conv_layers) {
        layer <- paste(layer, "<-torch::nn_conv2d(", conv_layers_size[[i - 1]], ",", conv_layers_size[[length(conv_layers_size)]], ",", conv_layers_kernel, ",", conv_layers_stride, ",", conv_layers_padding, ")", sep = "")
      } else {
        layer <- paste(layer, "<-torch::nn_conv2d(", conv_layers_size[[i - 1]], ",", conv_layers_size[[i]], ",", conv_layers_kernel, ",", conv_layers_stride, ",", conv_layers_padding, ")", sep = "")
      }

      conv_layers_definition[[i]] <- layer
    }

    # simulates the resolution loss across the convolutional layers
    for (i in 1:number_of_conv_layers) {
      if (any(sample_size < (conv_layers_kernel + conv_layers_padding))) {
        stop("Sample dimension is too small for the choosen configuration.")
      }

      if (is.numeric(pooling)) {
        sample_size[[1]] <- res_calculate(
          "layer",
          sample_size[[1]],
          conv_layers_kernel,
          conv_layers_stride,
          conv_layers_padding
        )

        sample_size[[1]] <- res_calculate(
          "pooling",
          sample_size[[1]],
          pooling
        )

        sample_size[[2]] <- res_calculate(
          "layer",
          sample_size[[2]],
          conv_layers_kernel,
          conv_layers_stride,
          conv_layers_padding
        )

        sample_size[[2]] <- res_calculate(
          "pooling",
          sample_size[[2]],
          pooling
        )
      } else if (!pooling) {
        sample_size[[1]] <- res_calculate(
          "layer",
          sample_size[[1]],
          conv_layers_kernel,
          conv_layers_stride,
          conv_layers_padding
        )

        sample_size[[2]] <- res_calculate(
          "layer",
          sample_size[[2]],
          conv_layers_kernel,
          conv_layers_stride,
          conv_layers_padding
        )
      }
    }
    #

    fc_layer_names <- c()
    fc_layers_definition <- rep("\n    self$", number_of_fc_layers, sep = "")
    for (i in 1:length(fc_layers_definition)) {
      layer <- fc_layers_definition[[i]] %>% paste("linear", i, sep = "")
      stripped <- strsplit(layer, split = "")[[1]]
      layer_name <- stripped[6:length(stripped)] %>% paste0(collapse = "")
      fc_layer_names <- append(fc_layer_names, layer_name)

      if (i == 1) {
        layer <- paste(layer, "<-torch::nn_linear(", sample_size[[1]] * sample_size[[2]] * conv_layers_size[[length(conv_layers_size)]], ",", fc_layers_size[[1]], ")", sep = "")
      } else if (i == number_of_fc_layers) {
        layer <- paste(layer, "<-torch::nn_linear(", fc_layers_size[[i - 1]], ",", fc_layers_size[[length(fc_layers_size)]], ")", sep = "")
      } else {
        layer <- paste(layer, "<-torch::nn_linear(", fc_layers_size[[i - 1]], ",", fc_layers_size[[i]], ")", sep = "")
      }

      fc_layers_definition[[i]] <- layer
    }

    out_layer <- paste("\n    self$output<-torch::nn_linear(", fc_layers_size[[length(fc_layers_size)]], ",", number_of_outputs, ")", sep = "")

    layers_definition <- conv_layers_definition %>%
      append(fc_layers_definition) %>%
      append(out_layer) %>%
      paste0(collapse = "")

    #
    if (batch_norm) {
      # convolutional layers
      batch_normalization_conv <- rep("\n    self$", number_of_conv_layers, sep = "")
      for (i in 1:length(batch_normalization_conv)) {
        bn_layer <- batch_normalization_conv[[i]] %>% paste("bn_conv", i, sep = "")
        stripped <- strsplit(bn_layer, split = "")[[1]]

        if (i == 1) {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm2d(num_features=", conv_layers_size[[1]], ")", sep = "")
        } else if (i == number_of_conv_layers) {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm2d(num_features=", conv_layers_size[[length(conv_layers_size)]], ")", sep = "")
        } else {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm2d(num_features=", conv_layers_size[[i]], ")", sep = "")
        }

        batch_normalization_conv[[i]] <- bn_layer
      }
      # fully connected layers
      batch_normalization_fc <- rep("\n    self$", number_of_fc_layers, sep = "")
      for (i in 1:length(batch_normalization_fc)) {
        bn_layer <- batch_normalization_fc[[i]] %>% paste("bn_linear", i, sep = "")
        stripped <- strsplit(bn_layer, split = "")[[1]]

        if (i == 1) {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm1d(num_features=", fc_layers_size[[1]], ")", sep = "")
        } else if (i == number_of_fc_layers) {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm1d(num_features=", fc_layers_size[[length(fc_layers_size)]], ")", sep = "")
        } else {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm1d(num_features=", fc_layers_size[[i]], ")", sep = "")
        }

        batch_normalization_fc[[i]] <- bn_layer
      }

      batch_normalization_conv <- batch_normalization_conv %>%
        paste0(collapse = "")

      batch_normalization_fc <- batch_normalization_fc %>%
        paste0(collapse = "")

      batch_normalization <- paste0(batch_normalization_conv, batch_normalization_fc)
      # batch_normalization <- paste("{", batch_normalization, "\n  }", sep = "")
      layers_definition <- paste0(layers_definition, batch_normalization)
    }

    layers_definition <- paste("{", layers_definition, "\n  }", sep = "")
    #

    foward_definition <- c("{\n    x %>%")
    conv_foward <- c()
    for (i in 1:length(conv_layer_names)) {
      seq <- paste("\n      ", paste(conv_layer_names[i], "()", sep = ""), " %>%\n      torch::nnf_relu() %>%", sep = "")

      if (batch_norm) {
        seq <- paste0(seq, "\n      self$bn_conv", i, "() %>%")
      }

      if (is.numeric(pooling)) {
        seq <- paste0(seq, "\n      nnf_avg_pool2d(", pooling, ") %>%")
      }

      if (dropout & dropout > 0 & dropout < 1 & is.numeric(dropout)) {
        seq <- paste0(seq, "\n      torch::nnf_dropout(p=", dropout, ") %>%")
      }

      conv_foward <- append(conv_foward, seq)
    }

    fc_foward <- c()
    for (i in 1:length(fc_layer_names)) {
      seq <- paste("\n      ", paste(fc_layer_names[i], "()", sep = ""), " %>%\n      torch::nnf_relu() %>%", sep = "")

      if (batch_norm) {
        seq <- paste0(seq, "\n      self$bn_linear", i, "() %>%")
      }

      if (dropout & dropout > 0 & dropout < 1 & is.numeric(dropout) & i < length(fc_layer_names)) {
        seq <- paste0(seq, "\n      torch::nnf_dropout(p=", dropout, ") %>%")
      }

      fc_foward <- append(fc_foward, seq)
    }

    foward_definition <- foward_definition %>%
      append(conv_foward) %>%
      append("\n      torch::torch_flatten(start_dim = 2) %>%") %>%
      append(fc_foward) %>%
      append("\n      self$output()\n  }") %>%
      paste0(collapse = "")

    arch <- sub("layer_definition", layers_definition, arch)
    arch <- sub("foward_definition", foward_definition, arch)

    torch::torch_manual_seed(13)

    net <- eval(parse(text = arch))

    if (verbose == TRUE) {
      cat(arch)
    }

    return_list <- list(net = net, arch = arch)

    return(return_list)
  }
