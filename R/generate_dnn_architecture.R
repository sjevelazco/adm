#' Generate architectures for Deep Neural Network
#'
#' @param number_of_features numeric. Value that specifies the number of features in the dataset.
#' @param number_of_outputs numeric. Value that specifies the number of outputs.
#' @param number_of_hidden_layers  numeric. Number of hidden layers in the neural network. Default 2.
#' @param hidden_layers_size numeric vector. Size of each hidden layer in the neural network. Default c(14, 7).
#' @param verbose logical. Whether to print the architecture. Default FALSE.
#' @param batch_norm logical. Whether to include batch normalization layers. Default TRUE.
#' @param dropout logical. Specifies whether dropout is included in the architecture. Default FALSE.
#'
#' @importFrom torch nn_module nn_linear nnf_relu torch_manual_seed
#'
#' @return A list containing:
#' \itemize{
#' \item net: a instantiated torch neural net.
#' \item arch: a string with a R expression to instantiate the neural network.
#' \item arch_dict: a list with a matrix describing the architecture structure.
#' }
#'
#' @seealso \code{\link{generate_arch_list}}, \code{\link{generate_cnn_architecture}}
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate a Deep Neural Network with:
#' dnn_arch <- generate_dnn_architecture(
#'   number_of_features = 8, # eight input variables
#'   number_of_outputs = 1, # one output
#'   number_of_hidden_layers = 5, # five layers between input and output
#'   hidden_layers_size = c(8, 16, 32, 16, 8), # of this size, respectively
#'   batch_norm = TRUE, # with batch normalization
#'   dropout = 0, # without dropout
#' )
#'
#' dnn_arch$net() # a torch net
#' dnn_arch$arch %>% cat() # the torch code to create it
#' dnn_arch$arch_dict # and a quick description of its structure
#' }
generate_dnn_architecture <-
  function(number_of_features = 7,
           number_of_outputs = 1,
           number_of_hidden_layers = 2,
           hidden_layers_size = c(14, 7),
           batch_norm = TRUE,
           dropout = 0,
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

      if (dropout > 0 & dropout < 1 & is.numeric(dropout) & i < length(layer_names)) {
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

    # Creating arch_dict
    arch_dict <- as.matrix(hidden_layers_size) %>%
      list()
    names(arch_dict) <- paste0(number_of_hidden_layers, "_layer_net")
    row.names(arch_dict[[1]]) <- paste0("layer_", seq(1, number_of_hidden_layers))


    return_list <- list(net = net, arch = arch, arch_dict = arch_dict)

    return(return_list)
  }
