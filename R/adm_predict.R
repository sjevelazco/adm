#' Spatial predictions from individual and ensemble models
#'
#' @description This function allows the geographical prediction of one or more models constructed
#' with the fit_ or tune_ function set, models fitted with esm_ function set (i.e., ensemble of
#' small models approach), or models constructed with fit_ensemble function. It can return
#' continuous or continuous and binary predictions for one or more thresholds
#'
#' @param models list of one or more models fitted with fit_ or tune_ functions. In case use models fitted with fit_ensemble or esm_ family function only one model could be used. Usage models = mglm or models = list(mglm, mraf, mgbm)
#' @param pred SpatRaster. Raster layer with predictor variables. Names of layers must exactly
#' match those used in model fitting.
#' @param nchunk integer. Number of chunks to split data used to predict models (i.e., SpatRaster
#' used in pred argument). Predicting models in chunks helps reduce memory requirements in cases
#' where models are predicted for large scales and high resolution. Default = 1
#' @param thr character. Threshold used to get binary suitability values (i.e., 0,1). It is possible
#' to use more than one threshold type. It is mandatory to use the same threshold/s used to fit the
#' models. The following threshold types are available:
#' \itemize{
#'   \item lpt: The highest threshold at which there is no omission.
#'   \item equal_sens_spec: Threshold at which the sensitivity and specificity are equal.
#'   \item max_sens_spec: Threshold at which the sum of the sensitivity and specificity is the
#'   highest (aka threshold that maximizes the TSS).
#'   \item max_jaccard: The threshold at which the Jaccard index is the highest.
#'   \item max_sorensen: The threshold at which the Sorensen index is highest.
#'   \item max_fpb: The threshold at which FPB is highest.
#'   \item sensitivity: Threshold based on a specified sensitivity value used to fit the models.
#'   \item all: All the threshold used in the model outputs used in 'models' argument will be used.
#'   }
#' Usage thr = c('lpt', 'max_sens_spec', 'max_jaccard'), thr=c('lpt', 'max_sens_spec',
#' 'sensitivity'), or thr='all'. If no threshold is specified (i.e., thr = NULL) function
#' will return continuous prediction only. Default NULL
#' @param con_thr logical. If true predictions with suitability values above threshold/s will be
#' returned. Default = FALSE
#' @param predict_area SpatVector, SpatialPolygon, or SpatialPolygonDataFrame. Spatial polygon
#' used for restring prediction into only a given region. Default = NULL
#' @param clamp logical. It is set with TRUE, predictors and features are restricted to the range
#' seen during model training. Only valid for Maxent model (see tune_mx and fit_mx). Default TRUE.
#' @param pred_type character. Type of response required available "link", "exponential",
#' "cloglog" and "logistic". Only valid for Maxent model (see tune_mx and fit_mx).
#' Default "cloglog".
#' @param training_data 
#' @param transform_negative 
#'
#' @return A list of SpatRaster with continuous and/or binary predictions
#'
#' @export
#'
#' @importFrom dplyr mutate across left_join pull bind_rows filter select
#' @importFrom gamlss predictAll
#' @importFrom kernlab predict
#' @importFrom stats predict median
#' @importFrom stringr str_detect
#' @importFrom terra vect crop mask as.data.frame is.factor rast app weighted.mean lapp crs
#' @importFrom torch dataset torch_tensor
#' 
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#' 
#' data("sppabund")
#' envar <- system.file("external/envar.tif", package = "adm")
#' envar <- terra::rast(envar)
#' 
#' # Extract data
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species one") %>%
#'   dplyr::select(species, ind_ha, x,y)
#' 
#' some_sp
#' 
#' some_sp <-
#'   adm_extract(
#'     data = some_sp,
#'     x = "x",
#'     y = "y",
#'     env_layer = envar
#'   )
#' 
#' # Partition
#' some_sp <- flexsdm::part_random(
#'   data = some_sp,
#'   pr_ab = "ind_ha",
#'   method = c(method = "rep_kfold", folds = 3, replicates = 3)
#' )
#' 
#' 
#' ## %######################################################%##
#' #                                                          #
#' ####          Create different type of models           ####
#' #                                                          #
#' ## %######################################################%##
#' # Fit some models
#' # require(gamlss)
#' # m1 <- gamlss::fitDist(some_sp$ind_ha, type="realline")
#' # m1$fits
#' # m1$failed
#' # 
#' # m1 <- gamlss(ind_ha ~ pb(elevation) + pb(sand) + pb(bio3) + pb(bio12), family=NO, data=some_sp)
#' # choosen_dist <- gamlss::chooseDist(m1, parallel="snow", ncpus=4, type="realAll")
#' 
#' a_gam <- fit_abund_gam(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("elevation", "sand", "bio3", "bio12"),
#'   predictors_f = "eco",
#'   partition = ".part",
#'   distribution = gamlss.dist::NO())
#' 
#' mraf <- fit_abund_raf(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("elevation", "sand", "bio3", "bio12"),
#'   partition = ".part",
#' )
#' 
#' mgbm <- fit_abund_gbm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("elevation", "sand", "bio3", "bio12"),
#'   partition = ".part",
#'   distribution = 
#' )
#'
#' 
#' ##%######################################################%##
#' #                                                          #
#' ####            ' ####      Predict models              ####
#' #                                                          #
#' ##%######################################################%##
#' 
#' # adm_predict can be used for predict one or more models fitted with fit_ or tune_ functions
#'
#' # a single model
#' ind_p <- sdm_predict(
#'   models = mglm,
#'   pred = somevar,
#'   thr = "max_fpb",
#'   con_thr = FALSE,
#'   predict_area = NULL
#' )
#'
#' # a list of models
#' list_p <- sdm_predict(
#'   models = list(mglm, mraf, mgbm),
#'   pred = somevar,
#'   thr = "max_fpb",
#'   con_thr = FALSE,
#'   predict_area = NULL
#' )
#'
#' # Predict an ensemble model
#' # (only is possilbe use one fit_ensemble)
#' ensemble_p <- sdm_predict(
#'   models = mensemble,
#'   pred = somevar,
#'   thr = "max_fpb",
#'   con_thr = FALSE,
#'   predict_area = NULL
#' )
#'
#' # Predict an ensemble of small models
#' # (only is possible to use one ensemble of small models)
#' small_p <- sdm_predict(
#'   models = msmall,
#'   pred = somevar,
#'   thr = "max_fpb",
#'   con_thr = FALSE,
#'   predict_area = NULL
#' )
#'
#' ##%######################################################%##
#' #                                                          #
#' ####              Predict model using chunks            ####
#' #                                                          #
#' ##%######################################################%##
#' # Predicting models in chunks helps reduce memory requirements in
#' # cases where models are predicted for large scales and high resolution
#'
#' ind_p <- sdm_predict(
#'   models = mglm,
#'   pred = somevar,
#'   thr = "max_fpb",
#'   con_thr = FALSE,
#'   predict_area = NULL,
#'   nchunk = 4
#' )
#'
#' }
#'
adm_predict <-
  function(models,
           pred,
           training_data = NULL,
           nchunk = 1,
           thr = NULL,
           con_thr = FALSE,
           predict_area = NULL,
           clamp = TRUE,
           pred_type = "cloglog",
           transform_negative = TRUE) {
    . <- model <- threshold <- thr_value <- self <- NULL
    
    # TODO write codes to predict CNN ANN adapt GLM and GAM to use gamlss
    
    if (is.null(names(models))) {
      message("Predicting list of individual models")
      ensembles <- NULL
      esm <- NULL
    } else if (all(names(models) %in% c("models", "thr_metric", "predictors", "performance"))) {
      message("Predicting ensembles")
      ensembles <- models
      models <- NULL
      esm <- NULL
    } else if (all(names(models) %in% c("esm_model", "predictors", "performance"))) {
      message("Predicting ensemble of small models")
      esm <- models
      models <- NULL
      ensembles <- NULL
    } else {
      message("Predicting individual models")
      models <- list(models)
      ensembles <- NULL
      esm <- NULL
    }
    
    #### Prepare datasets ####
    # Crop and mask projection area
    if (!is.null(predict_area)) {
      if (class(predict_area) %in% c("SpatialPolygons", "SpatialPolygonsDataFrame")) {
        predict_area <- terra::vect(predict_area)
      }
      pred <-
        pred %>%
        terra::crop(., predict_area) %>%
        terra::mask(., predict_area)
    }
    
    #### Model predictions
    if (!is.null(models)) {
      # Prepare model list
      m <- lapply(models, function(x) {
        x[[1]]
      })
      
      names(m) <- paste0("m_", 1:length(m))
      
      # Extract model names object
      clss <- sapply(m, function(x) {
        class(x)[1]
      }) %>%
        tolower() %>%
        gsub(".formula", "", .)
      
      if (any(lapply(models, function(x){class(x[[1]])[1]}) %>% unlist == "gamlss")){
        indx <- which(lapply(models, function(x){class(x[[1]])[1]}) %>% unlist == "gamlss")
        gamlss_classes <- lapply(models[indx], function(x){
          x$predictors$model
        }) %>% unlist()
        
        clss[indx] <- paste0(clss[indx],"_",gamlss_classes)
      }
    }
    
    # Transform raster to data.frame
    
    # if(chunk){
    cell <- terra::as.data.frame(pred, cells = TRUE, na.rm = TRUE)[,"cell"]
    set <-
      c(seq(1, length(cell), length.out = nchunk) # length.out = nchunk
        %>% round(),
        length(cell) + 1
      )
    # } #else {
    #   pred_df <-
    #   terra::as.data.frame(pred, cells = FALSE, na.rm = TRUE)
    # }
    
    r <- pred[[!terra::is.factor(pred)]][[1]]
    r[!is.na(r)] <- NA
    
    # Create list for storing raster for current condition
    model_c <- as.list(names(m))
    names(model_c) <- names(m)
    for(i in seq_along(model_c)){
      model_c[[i]] <- r
    }
    
    
    
    # Write here the loop
    for(CH in seq_len((length(set) - 1))){
      rowset <- set[CH]:(set[CH + 1] - 1)
      rowset <- cell[rowset]
      pred_df <- pred[rowset]
      rownames(pred_df) <- rowset
      
      
      ## %######################################################%##
      #                                                          #
      ####          Prediction for different models           ####
      #                                                          #
      ## %######################################################%##
      
      #### xgboost models ####
      # wm <- which(clss == "xgb.booster")
      # if (length(clss) > 0){
      #   wm <- names(wm)
      #   for (i in wm) {
      #     r <- pred[[!terra::is.factor(pred)]][[1]]
      #     r[!is.na(r)] <- NA
      #     r[as.numeric(rownames(pred_df))] <-
      #     
      #     model_c[[i]][rowset] <- r[rowset]
      #   }
      # }
      
      #### dnn models ####
      wm <- which(clss == "luz_module_fitted")
      if (length(wm) > 0){
        wm <- names(wm)
        
        # create_dataset definition
        create_dataset <- torch::dataset(
          "dataset",
          initialize = function(df, response_variable=0) {
            self$df <- df
          },
          .getitem = function(index) {
            x <- torch::torch_tensor(as.numeric(self$df[index, ]))
            list(x = x)
          },
          .length = function() {
            nrow(self$df)
          }
        )
        
        pred_dataset <- create_dataset(pred_df)
        for (i in wm) {
          r <- pred[[!terra::is.factor(pred)]][[1]]
          r[!is.na(r)] <- NA
          r[as.numeric(rownames(pred_df))] <-
            suppressMessages(stats::predict(m[[i]], pred_dataset) %>% as.numeric())
          
          model_c[[i]][rowset] <- r[rowset]
        }
      }
      
      #### gam and glm ####
      wm <- clss[stringr::str_detect(clss,"gamlss")]
      if (length(wm) > 0){
        wm <- names(wm)
        for (i in wm) {
          r <- pred[[!terra::is.factor(pred)]][[1]]
          r[!is.na(r)] <- NA
          r[as.numeric(rownames(pred_df))] <-
            suppressMessages(predict(m[[i]], newdata = pred_df, data = training_data, type = "response"))
          
          model_c[[i]][rowset] <- r[rowset]
        }
      }
      
      #### gbm models ####
      wm <- which(clss == "gbm")
      if (length(wm) > 0) {
        wm <- names(wm)
        for (i in wm) {
          r <- pred[[!terra::is.factor(pred)]][[1]]
          r[!is.na(r)] <- NA
          r[as.numeric(rownames(pred_df))] <-
            suppressMessages(gamlss::predictAll(m[[i]], pred_df, type = "response"))
          
          model_c[[i]][rowset] <- r[rowset]
        }
      }
      
      
      #### randomforest class ####
      wm <- which(clss == "randomforest")
      if (length(wm) > 0) {
        wm <- names(wm)
        for (i in wm) {
          r <- pred[[!terra::is.factor(pred)]][[1]]
          r[!is.na(r)] <- NA
          
          # Test factor levels
          f <-
            m[[i]]$forest$xlevels[sapply(m[[i]]$forest$xlevels, function(x) {
              class(x) == "character"
            })]
          if (length(f) > 0) {
            for (ii in 1:length(f)) {
              vf <- f[[ii]] %>%
                unique()
              vf2 <- pred_df[, names(f[ii])] %>% unique()
              vfilter <- list()
              if (sum(!vf2 %in% vf) > 0) {
                vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
              }
            }
            if (length(vfilter) > 0) {
              if (length(vfilter) > 1) {
                vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
              } else {
                vfilter <- vfilter[[1]]
              }
            } else {
              vfilter <- 0
            }
          } else {
            vfilter <- 0
          }
          
          
          if (sum(vfilter) > 0) {
            v <- rep(0, nrow(pred_df))
            v[!vfilter] <-
              suppressMessages(stats::predict(m[[i]], pred_df[!vfilter, ] %>%
                               dplyr::mutate(dplyr::across(
                                 .cols = names(f),
                                 .fns = ~ droplevels(.)
                               )),
                             type = "response")
              )
            r[as.numeric(rownames(pred_df))] <- v
            rm(v)
          } else {
            r[as.numeric(rownames(pred_df))] <-
              suppressMessages(stats::predict(m[[i]], pred_df, type = "response"))
          }
          
          model_c[[i]][rowset] <- r[rowset]
        }
      }
      
      #### ksvmj class ####
      wm <- which(clss == "ksvm")
      if (length(wm) > 0) {
        wm <- names(wm)
        for (i in wm) {
          r <- pred[[!terra::is.factor(pred)]][[1]]
          r[!is.na(r)] <- NA
          
          # Test factor levels
          f_n <- which(sapply(pred_df, class) == "factor") %>% names()
          f_n2 <- m[[i]]@xmatrix[[1]] %>% colnames()
          f <- lapply(f_n, function(x) {
            gsub(x, "", grep(x, f_n2, value = TRUE))
          })
          names(f) <- f_n
          
          if (length(f) > 0) {
            for (ii in 1:length(f)) {
              vf <- f[[ii]] %>%
                unique()
              vf2 <- pred_df[, names(f[ii])] %>% unique()
              vfilter <- list()
              if (sum(!vf2 %in% vf) > 0) {
                vfilter[[ii]] <- !pred_df[, names(f[ii])] %in% vf
              }
            }
            if (length(vfilter) > 0) {
              if (length(vfilter) > 1) {
                vfilter <- vapply(do.call("rbind", vfilter), any, logical(1))
              } else {
                vfilter <- vfilter[[1]]
              }
            } else {
              vfilter <- 0
            }
          } else {
            vfilter <- 0
          }
          
          if (sum(vfilter) > 0) {
            v <- rep(0, nrow(pred_df))
            v[!vfilter] <-
              kernlab::predict(m[[i]], pred_df[!vfilter, ] %>%
                                 dplyr::mutate(dplyr::across(
                                   .cols = names(f),
                                   .fns = ~ droplevels(.)
                                 )), type = "response")[, 2]
            r[as.numeric(rownames(pred_df))] <- v
            rm(v)
          } else {
            r[as.numeric(rownames(pred_df))] <-
              kernlab::predict(m[[i]], pred_df, type = "response")[, 2]
          }
          model_c[[i]][rowset] <- r[rowset]
        }
      }
    }
    
    rm(pred_df)
    
    df <- data.frame(
      alg = c(
        "luz_module_fitted",
        "gamlss_gam",
        "gamlss_glm",
        "gbm",
        "nnet",
        "randomforest",
        "ksvm"
      ),
      names = c("dnn", "gam", "glm", "gbm", "net", "raf", "svm")
    )
    
    names(model_c) <-
      dplyr::left_join(data.frame(alg = clss), df, by = "alg")[, 2]
    model_c <- mapply(function(x, n) {
      names(x) <- n
      x
    }, model_c, names(model_c))
    
    
    ## %######################################################%##
    #                                                          #
    ####                  Predict ensemble                  ####
    #                                                          #
    ## %######################################################%##
    if (!is.null(ensembles)) {
      model_c <- terra::rast(model_c) # stack individual models
      ens_perf <- ensembles[["performance"]] # get performance of ensembles
      ens_method <- ens_perf %>%
        dplyr::pull(model) %>%
        unique() # get ensemble methods
      single_model_perf <- lapply(ensembles[[1]], function(x) { # Get performance of individual model used for performing ensemble
        x[["performance"]]
      }) %>% dplyr::bind_rows()
      
      # Threshold and metric values for performing some ensembles
      if (any(ens_method %in% c("meanw", "meansup", "meanthr"))) {
        weight_data <- single_model_perf %>%
          dplyr::filter(threshold == ensembles$thr_metric[1]) %>%
          dplyr::select(model, thr_value, ensembles$thr_metric[2])
      }
      
      ensemble_c <- as.list(ens_method)
      names(ensemble_c) <- ensemble_c
      
      if (any("mean" == ens_method)) {
        ensemble_c[["mean"]] <- terra::app(model_c, fun = mean, cores = 1)
      }
      
      if (any("meanw" == ens_method)) {
        ensemble_c[["meanw"]] <- terra::weighted.mean(model_c, weight_data[[3]])
      }
      
      if (any("meansup" == ens_method)) {
        ensemble_c[["meansup"]] <-
          terra::app(model_c[[which(weight_data[[3]] >= mean(weight_data[[3]]))]],
                     fun = mean, cores = 1
          )
      }
      
      if (any("meanthr" == ens_method)) {
        mf1 <- function(x, y) {
          terra::lapp(x, function(x) ifelse(x >= y, x, 0))
        }
        model_c2 <- mapply(mf1, model_c, weight_data[[2]], SIMPLIFY = FALSE) %>% terra::rast()
        ensemble_c[["meanthr"]] <- terra::app(model_c2, fun = mean, cores = 1)
        rm(model_c2)
      }
      
      if (any("median" == ens_method)) {
        ensemble_c[["median"]] <- terra::app(model_c, fun = stats::median, cores = 1)
      }
      
      ensemble_c <- mapply(function(x, y) {
        names(x) <- y
        x
      }, ensemble_c, ens_method, SIMPLIFY = FALSE)
      
      model_c <- ensemble_c
      models <- split(ensembles[["performance"]], ensembles[["performance"]]$model)
      models <- lapply(models, function(x) {
        x <- list(x)
        names(x) <- "performance"
        x
      })
      rm(ensembles)
      rm(ensemble_c)
    }
    
    ## %######################################################%##
    #                                                          #
    ####       Predict ensemble of small models (esm)       ####
    #                                                          #
    ## %######################################################%##
    
    if (!is.null(esm)) {
      model_c <- terra::rast(model_c) # stack individual models
      weight_data <- esm$esm_model %>%
        names() %>%
        as.numeric() # get performance of esm
      
      ensemble_c <-
        terra::app(model_c * weight_data, fun = sum, cores = 1)
      ensemble_c <- ensemble_c/sum(weight_data)
      names(ensemble_c) <- paste0("esm_", unique(names(model_c)))
      
      model_c <- list(ensemble_c)
      models <- list("performance" = esm["performance"])
      
      rm(esm)
      rm(ensemble_c)
    }
    
    
    ## %######################################################%##
    #                                                          #
    ####              Get binary predictions                ####
    #                                                          #
    ## %######################################################%##
    # TODO
    if (is.null(thr)) {
      
      if (transform_negative) {
        for (i in 1:length(model_c)) {
          x <- model_c[[i]]
          x[x<0] <- 0
          model_c[[i]] <- x
        }
        rm(x)
      }
      
      return(model_c)
    } else {
      if (any("all" == thr)) {
        thr_df <- lapply(models, function(x) {
          x[["performance"]]
        })
      } else {
        thr_df <- lapply(models, function(x) {
          x[["performance"]] %>%
            dplyr::filter(threshold %in% thr)
        })
      }
      
      model_b <- list()
      
      for (i in 1:length(model_c)) {
        model_b[[i]] <-
          lapply(thr_df[[i]]$thr_value, function(x) {
            model_c[[i]] >= x
          }) %>% terra::rast()
        names(model_b[[i]]) <- names(thr_df[[i]]$thr_value)
      }
      
      names(model_b) <- names(model_c)
      
      # Return suitability values above thresholds
      if (con_thr) {
        for (i in 1:length(model_b)) {
          nms <- names(model_b[[i]])
          model_b[[i]] <- model_b[[i]] * model_c[[i]]
          names(model_b[[i]]) <- nms
        }
      }
      
      mf2 <- function(x, x2) {
        terra::rast(list(x, x2))
      }
      
      result <- mapply(mf2, model_c, model_b, SIMPLIFY = FALSE)
      if(grepl("esm_", names(model_c[[1]]))) {
        names(result) <- names(model_c[[1]])
      }
      for(f in 1:length(result)){
        terra::crs(result[[f]]) <- terra::crs(pred)
      }
      
      return(result)
    }
  }
