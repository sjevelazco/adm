#' Select arch list
#'
#' @param arch_list
#' @param type
#' @param method
#' @param n_samples
#' @param min_max
#'
#' @importFrom dplyr bind_rows group_by slice_min tibble select rename left_join join_by
#' @importFrom stats quantile
#' @importFrom stringr str_extract_all
#'
#' @return
#' @export
#'
#' @examples
select_arch_list <-
  function(arch_list,
           type = c("dnn", "cnn"),
           method = "percentile", # TODO
           n_samples = 1, # TODO
           min_max = TRUE # TODO
  ) {
    parameters <- name <- combination <- old_name <- NULL

    architectures <- arch_list

    net_size_df <- list()
    for (arch in architectures$arch_list %>% names()) {
      net <- architectures$arch_list[[arch]]
      net_test <- net()

      net_size <- c()
      for (param in net_test$parameters) {
        net_size <- append(net_size, param %>% length())
      }

      net_size <- sum(net_size)
      net_df <- data.frame(arch = arch, parameters = net_size) %>%
        list()
      net_size_df <- append(net_size_df, net_df)
    }

    net_size_df <- dplyr::bind_rows(net_size_df)

    arch_indexes <- stringr::str_extract_all(net_size_df[["arch"]], "\\d+")
    deep <- c()
    if (type == "dnn") {
      for (x in arch_indexes) {
        deep <- append(deep, x[[1]] %>% as.numeric())
      }
    } else if (type == "cnn") {
      for (x in arch_indexes) {
        deep <- append(deep, x[[1]] %>% as.numeric() + x[[2]] %>% as.numeric())
      }
    }

    net_size_df$deep <- deep

    if (min_max) {
      minimal_net <- min(net_size_df[["parameters"]])
      maximal_net <- max(net_size_df[["parameters"]])

      min_row <- net_size_df[which(net_size_df[["parameters"]] == minimal_net), ]
      max_row <- net_size_df[which(net_size_df[["parameters"]] == maximal_net), ]

      net_size_df <- net_size_df[
        -c(
          which(net_size_df[["parameters"]] == minimal_net),
          which(net_size_df[["parameters"]] == maximal_net)
        ),
      ]
    }

    # TODO - FAZ SÓ PERCENTIL
    percentile <- stats::quantile(net_size_df$parameters, probs = seq(0, 1, by = 0.1))
    net_size_df$percentile <-
      cut(net_size_df$parameters,
        breaks = percentile,
        include.lowest = TRUE,
        labels = FALSE
      )

    net_size_df <- net_size_df %>%
      dplyr::group_by(percentile, deep) %>%
      dplyr::slice_min(dplyr::tibble(parameters), n = n_samples)

    if (min_max) {
      net_size_df <-
        dplyr::bind_rows(
          min_row,
          net_size_df %>%
            as.data.frame() %>%
            dplyr::select(-percentile),
          max_row
        )
    } else {
      net_size_df <-
        dplyr::bind_rows(
          net_size_df %>%
            as.data.frame() %>%
            dplyr::select(-percentile)
        )
    }

    new_arch_list <- architectures$arch_list[net_size_df$arch]

    layer_classes <- c()
    combination_numbers <- c()
    if (type == "dnn") {
      for (arch in names(new_arch_list)) {
        arch_indexes <- stringr::str_extract_all(arch, "\\d+")[[1]] %>%
          as.numeric()
        layer_class <- paste0(arch_indexes[[1]], "_layer_net")
        combination_number <- arch_indexes[[2]]

        layer_classes <- append(layer_classes, layer_class)
        combination_numbers <- append(combination_numbers, combination_number)
      }
    } else if (type == "cnn") {
      for (arch in names(new_arch_list)) {
        arch_indexes <- stringr::str_extract_all(arch, "\\d+")[[1]] %>%
          as.numeric()
        layer_class <- paste0("conv", arch_indexes[[1]], "-fc", arch_indexes[[2]], "-net")
        combination_number <- arch_indexes[[3]]

        layer_classes <- append(layer_classes, layer_class)
        combination_numbers <- append(combination_numbers, combination_number)
      }
    }

    class_comb_df <- data.frame(
      class = layer_classes,
      combination = combination_numbers,
      name = names(new_arch_list)
    )

    classes_df <- list()
    new_arch_dict <- list()
    for (layer_class in unique(class_comb_df$class)) {
      if (type == "dnn") {
        class_df <- class_comb_df[class_comb_df$class == layer_class, ]
        layer_class_dict <- architectures$arch_dict[[layer_class]][, class_df[["combination"]]] %>%
          as.matrix() %>%
          list()

        class_number <- stringr::str_extract_all(layer_class, "\\d+")[[1]]
        class_df$new_name <- paste0(paste0("arch-", class_number, "-"), 1:ncol(layer_class_dict[[1]]))
        classes_df <- append(classes_df, list(class_df))

        names(layer_class_dict) <- layer_class
        new_arch_dict <- append(new_arch_dict, layer_class_dict)
      } else if (type == "cnn") {
        class_df <- class_comb_df[class_comb_df$class == layer_class, ]
        layer_class_dict <- architectures$arch_dict[[layer_class]][, class_df[["combination"]]] %>%
          as.matrix() %>%
          list()

        class_number <- stringr::str_extract_all(layer_class, "\\d+")[[1]]
        class_df$new_name <- paste0(paste0("arch-conv", class_number[[1]], "-fc", class_number[[2]], "-"), 1:ncol(layer_class_dict[[1]]))
        classes_df <- append(classes_df, list(class_df))

        names(layer_class_dict) <- layer_class
        new_arch_dict <- append(new_arch_dict, layer_class_dict)
      }
    }
    class_comb_df <- dplyr::bind_rows(classes_df) %>%
      dplyr::rename(old_name = name)

    new_arch_list <- new_arch_list[match(class_comb_df[["old_name"]], names(new_arch_list))]
    names(new_arch_list) <- class_comb_df[["new_name"]]

    class_comb_df <- class_comb_df %>% dplyr::select(-combination)
    net_size_df <- net_size_df %>%
      dplyr::rename(old_name = arch)

    final_df <- dplyr::left_join(class_comb_df, net_size_df, by = dplyr::join_by(old_name))

    return_list <- list(
      arch_list = new_arch_list,
      arch_dict = new_arch_dict,
      changes = final_df
    )

    return(return_list)
  }
