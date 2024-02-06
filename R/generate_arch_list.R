generate_arch_list <-
  function(number_of_features,
           number_of_outputs,
           n_layers = seq(1,3),
           n_neurons = seq(2,14,2)){
    
    arch_dict <- list()
    for (i in n_layers) {
      neuronal_list <- rep(list(n_neurons),i)
      org <- do.call(expand.grid, neuronal_list) %>% 
        t() %>%
        as.matrix()
      
      row.names(org) <- paste0("layer_",seq(i))

      temp_list <- list(org)
      names(temp_list) <- paste0(i,"_layer_net")
      arch_dict <- append(arch_dict, temp_list)
    }
    
    arch_list <- list()
    for (arch in arch_dict) {
      number_of_hidden_layers <- nrow(arch)
      for (comb in 1:ncol(arch)) {
        hidden_layers_size <- arch[,comb]
        net <- generate_dnn_architecture(
          number_of_features = number_of_features,
          number_of_outputs = number_of_outputs,
          number_of_hidden_layers = number_of_hidden_layers,
          hidden_layers_size = hidden_layers_size
        )
        
        l <- list(net$net)
        names(l) <- paste0("arch-",number_of_hidden_layers,"-",comb)
        arch_list <- append(arch_list,l)
      }
    }
    
    return_list <- list(arch_list = arch_list,
                        arch_dict = arch_dict)
    
    return(return_list)
  }
