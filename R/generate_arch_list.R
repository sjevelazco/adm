#' Generate architecture list for Deep Neural Network and Convolutional Neural Network
#'
#' @description This function generates a list of architectures for either a Deep Neural Netwokd (DNN) or
#' a Convolutional Neural Network (CNN).
#'
#' @param type string. Specifies the type the network. The valid inputs "dnn" and "cnn".
#' @param number_of_features numeric. Value that specifies the number of features in the dataset.
#' @param number_of_outputs numeric. Value that specifies the number of outputs.
#' @param n_layers numeric. Vector that specifies the number of layers in the networks. Default value are 1 and 2.
#' @param n_neurons  vector. Specifies the number of neurons each layer. Default  7.
#' @param sample_size vector. Specifies the size. Default c(11, 11)
#' @param number_of_fc_layers numeric. Specifies the number of fully connected layers. Default 1.
#' @param fc_layers_size vector. Specifies the size of the fully connected layers. Default 14.
#' @param conv_layers_kernel numeric. Specifies the kernel size for layers. Default 3.
#' @param conv_layers_stride numeric. Specifies the stride for the convolutional layers. Default 1.
#' @param conv_layers_padding  numeric. Specifies the padding for the convolutional layers. Default 0.
#' @param pooling numeric. Specifies 2D average pooling kernel size. Default NULL
#' @param batch_norm logical. Specifies whether batch normalization is included in the architecture. Default TRUE.
#' @param dropout Numeric. The probability (p) of randomly zeroing elements of the input tensor during training to prevent overfitting. Must be between 0 (no dropout) and 1 (all inputs zeroed). Default is 0 (no dropout).
#'
#' @importFrom dplyr rename
#' @importFrom tidyr unnest_wider
#'
#' @return A list containing:
#' \itemize{
#' \item arch_list: a list of generated architectures.
#' \item arch_dict: a list of architecture dictionaries.
#' }
#'
#' @seealso \code{\link{generate_dnn_architecture}}, \code{\link{generate_cnn_architecture}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generating architectures for DNN, using batch normalization and dropout
#' dnn_archs <- generate_arch_list(
#'   type = "dnn",
#'   number_of_features = 4,
#'   number_of_outputs = 1,
#'   n_layers = c(2, 3, 4),
#'   n_neurons = c(8, 16, 32, 64),
#'   batch_norm = TRUE,
#'   dropout = 0.2
#' )
#'
#' dnn_archs$arch_dict # Matrices describing the networks
#' length(dnn_archs$arch_list) # Generated 336 DNN architectures
#'
#' # Generating architectures for CNN, using batch normalization, dropout and average pooling
#' # Note that arguments meaning change with the context
#'
#' cnn_archs <- generate_arch_list(
#'   type = "cnn",
#'   number_of_features = 4,
#'   number_of_outputs = 1,
#'   n_layers = c(2, 3, 4), # now convolutional layers
#'   n_neurons = c(8, 16, 32, 64),
#'   sample_size = c(11, 11),
#'   number_of_fc_layers = c(2, 4), # fully connected layers
#'   fc_layers_size = c(16, 8),
#'   conv_layers_kernel = 3,
#'   conv_layers_stride = 1,
#'   conv_layers_padding = 0,
#'   pooling = 1,
#'   batch_norm = TRUE,
#'   dropout = 0.2
#' )
#'
#' cnn_archs$arch_dict # Matrices describing the networks
#' length(cnn_archs$arch_list) # Generated 6720 CNN architectures
#'
#' # The list size can be easily and greatly reduced with select_arch_list
#'
#' dnn_archs_redux <- dnn_archs %>% select_arch_list(
#'   type = c("dnn"),
#'   method = "percentile",
#'   n_samples = 1,
#'   min_max = TRUE # Keep the network with the minimum and maximum number of parameters
#' )
#'
#' length(dnn_archs_redux$arch_list) # from 336 to 29 architectures
#'
#' cnn_archs_redux <- cnn_archs %>% select_arch_list(
#'   type = c("cnn"),
#'   method = "percentile",
#'   n_samples = 1,
#'   min_max = TRUE # Keep the network with the minimum and maximum number of parameters
#' )
#'
#' length(cnn_archs_redux$arch_list) # from 6720 to 77 architectures
#' }
generate_arch_list <-
  function(type,
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
           dropout = 0) {
    Var1 <- Var2 <- NULL

    if (type == "dnn") {
      arch_dict <- list()
      for (i in n_layers) {
        neuronal_list <- rep(list(n_neurons), i)
        org <- do.call(expand.grid, neuronal_list) %>%
          t() %>%
          as.matrix()

        row.names(org) <- paste0("layer_", seq(i))

        temp_list <- list(org)
        names(temp_list) <- paste0(i, "_layer_net")
        arch_dict <- append(arch_dict, temp_list)
      }

      arch_list <- list()
      for (arch in arch_dict) {
        number_of_hidden_layers <- nrow(arch)
        for (comb in 1:ncol(arch)) {
          hidden_layers_size <- arch[, comb]
          net <- generate_dnn_architecture(
            number_of_features = number_of_features,
            number_of_outputs = number_of_outputs,
            number_of_hidden_layers = number_of_hidden_layers,
            hidden_layers_size = hidden_layers_size,
            batch_norm = batch_norm,
            dropout = dropout
          )

          l <- list(net$net)
          names(l) <- paste0("arch-", number_of_hidden_layers, "-", comb)
          arch_list <- append(arch_list, l)
        }
      }

      return_list <- list(
        arch_list = arch_list,
        arch_dict = arch_dict
      )
    } else if (type == "cnn") {
      if (!is.vector(sample_size)) {
        stop("sample_size expected to be a vector.")
      } else {
        conv_dict <- list()
        for (i in n_layers) {
          neuronal_list <- rep(list(n_neurons), i)
          org <- do.call(expand.grid, neuronal_list) %>%
            t() %>%
            as.matrix()

          row.names(org) <- paste0("conv_", seq(i))

          temp_list <- list(org)
          names(temp_list) <- paste0(i, "_layer_net")
          conv_dict <- append(conv_dict, temp_list)
        }

        fc_dict <- list()
        for (i in number_of_fc_layers) {
          neuronal_list <- rep(list(fc_layers_size), i)
          org <- do.call(expand.grid, neuronal_list) %>%
            t() %>%
            as.matrix()

          row.names(org) <- paste0("fc_", seq(i))

          temp_list <- list(org)
          names(temp_list) <- paste0(i, "_layer_net")
          fc_dict <- append(fc_dict, temp_list)
        }

        arch_dict <- list()
        for (conv_set in conv_dict) {
          for (fc_set in fc_dict) {
            conv_list <- list()
            for (i in 1:ncol(conv_set)) {
              v <- as.vector(conv_set[, i])
              conv_list <- append(conv_list, list(v))
            }

            fc_list <- list()
            for (i in 1:ncol(fc_set)) {
              v <- as.vector(fc_set[, i])
              fc_list <- append(fc_list, list(v))
            }

            connected_set <- expand.grid(conv_list, fc_list) %>%
              dplyr::rename(conv = Var1, fc = Var2) %>%
              tidyr::unnest_wider(c("conv", "fc"), names_sep = "_") %>%
              t() %>%
              as.matrix()

            l <- list(connected_set)
            names(l) <- paste0("conv", nrow(conv_set), "-fc", nrow(fc_set), "-net")
            arch_dict <- append(arch_dict, l)
          }
        }

        arch_list <- list()
        for (arch in arch_dict) {
          layers <- rownames(arch)
          conv_layers_names <- layers[grep("^conv_", layers)]
          fc_layers_names <- layers[grep("^fc_", layers)]
          number_of_conv_layers <- length(conv_layers_names)
          number_of_fc_layers <- length(fc_layers_names)
          for (comb in 1:ncol(arch)) {
            conv_layers_size <- arch[, comb][conv_layers_names]
            fc_layers_size <- arch[, comb][fc_layers_names]
            net <- generate_cnn_architecture(
              number_of_features = number_of_features,
              number_of_outputs = number_of_outputs,
              sample_size = sample_size,
              number_of_conv_layers = number_of_conv_layers,
              conv_layers_size = conv_layers_size,
              conv_layers_kernel = conv_layers_kernel,
              conv_layers_stride = conv_layers_stride,
              conv_layers_padding = conv_layers_padding,
              number_of_fc_layers = number_of_fc_layers,
              fc_layers_size = fc_layers_size,
              pooling = pooling,
              batch_norm = batch_norm,
              dropout = dropout,
              verbose = FALSE
            )

            l <- list(net$net)
            names(l) <- paste0("arch-conv", number_of_conv_layers, "-fc", number_of_fc_layers, "-", comb)
            arch_list <- append(arch_list, l)
          }
        }

        return_list <- list(
          arch_list = arch_list,
          arch_dict = arch_dict
        )
      }
    } else {
      stop('type is expected to be one of "cnn" or "dnn".')
    }

    return(return_list)
  }
