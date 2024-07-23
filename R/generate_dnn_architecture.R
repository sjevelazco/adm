#' Generate Deep Neural Network architecture
#'
#' @param number_of_features
#' @param number_of_outputs
#' @param number_of_hidden_layers
#' @param hidden_layers_size
#' @param verbose
#'
#' @importFrom torch nn_module nn_linear nnf_relu torch_manual_seed
#'
#' @return
#' @export
#'
#' @examples
generate_dnn_architecture <-
  function(number_of_features = 7,
           number_of_outputs = 1,
           number_of_hidden_layers = 2,
           hidden_layers_size = c(14, 7),
           batch_norm = TRUE,
           dropout = 0.2,
           verbose = FALSE) {
    arch <- "net <- torch::nn_module(
  'neural_net',
  initialize = function() layer_definition,
  forward = function(x) foward_definition
)"

    layer_names <- c()
    layers_definition <- rep("\n    self$", number_of_hidden_layers, sep = "")
    for (i in 1:length(layers_definition)) {
      layer <- layers_definition[[i]] %>% paste("linear", i, sep = "")
      stripped <- strsplit(layer, split = "")[[1]]
      layer_name <- stripped[6:length(stripped)] %>% paste0(collapse = "")
      layer_names <- append(layer_names, layer_name)

      if (i == 1) {
        layer <- paste(layer, "<-torch::nn_linear(", number_of_features, ",", hidden_layers_size[[1]], ")", sep = "")
      } else if (i == number_of_hidden_layers) {
        layer <- paste(layer, "<-torch::nn_linear(", hidden_layers_size[[i - 1]], ",", hidden_layers_size[[length(hidden_layers_size)]], ")", sep = "")
      } else {
        layer <- paste(layer, "<-torch::nn_linear(", hidden_layers_size[[i - 1]], ",", hidden_layers_size[[i]], ")", sep = "")
      }

      layers_definition[[i]] <- layer
    }
    out_layer <- paste("\n    self$output<-torch::nn_linear(", hidden_layers_size[[length(hidden_layers_size)]], ",", number_of_outputs, ")", sep = "")

    layers_definition <- layers_definition %>%
      append(out_layer) %>%
      paste0(collapse = "")

    ##
    if (batch_norm) {
      batch_normalization <- rep("\n    self$", number_of_hidden_layers, sep = "")
      for (i in 1:length(batch_normalization)) {
        bn_layer <- batch_normalization[[i]] %>% paste("bn", i, sep = "")
        stripped <- strsplit(bn_layer, split = "")[[1]]
        # bn_layer_name <- stripped[6:length(stripped)] %>% paste0(collapse = "")
        # bn_layer_names <- append(layer_names, layer_name)

        if (i == 1) {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm1d(num_features=", hidden_layers_size[[1]], ")", sep = "")
        } else if (i == number_of_hidden_layers) {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm1d(num_features=", hidden_layers_size[[length(hidden_layers_size)]], ")", sep = "")
        } else {
          bn_layer <- paste(bn_layer, "<-torch::nn_batch_norm1d(num_features=", hidden_layers_size[[i]], ")", sep = "")
        }

        batch_normalization[[i]] <- bn_layer
      }
      batch_normalization <- batch_normalization %>%
        paste0(collapse = "")
      # batch_normalization <- paste("{", batch_normalization, "\n  }", sep = "")
      layers_definition <- paste0(layers_definition, batch_normalization)
    }

    layers_definition <- paste("{", layers_definition, "\n  }", sep = "")

    ##
    foward_definition <- c("{\n    x %>%")
    for (i in 1:length(layer_names)) {
      seq <- paste("\n      ", paste(layer_names[i], "()", sep = ""), " %>%\n      torch::nnf_relu() %>%", sep = "")

      if (batch_norm) {
        seq <- paste0(seq, "\n      self$bn", i, "() %>%")
      }

      if (dropout & dropout > 0 & dropout < 1 & is.numeric(dropout) & i < length(layer_names)) {
        seq <- paste0(seq, "\n      torch::nnf_dropout(p=", dropout, ") %>%")
      }

      foward_definition <- append(foward_definition, seq)
    }
    foward_definition <- foward_definition %>%
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
