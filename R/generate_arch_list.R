#' Title
#'
#' @param type
#' @param number_of_features
#' @param number_of_outputs
#' @param n_layers
#' @param n_neurons
#' @param sample_size
#' @param number_of_fc_layers
#' @param fc_layers_size
#'
#' @importFrom dplyr rename
#' @importFrom tidyr unnest_wider
#'
#' @return
#' @export
#'
#' @examples
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
           pooling = FALSE,
           batch_norm = TRUE,
           dropout = FALSE) {
    
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
