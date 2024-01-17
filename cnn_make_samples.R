# cropping.hood
croppin.hood = function(occ, longitude, latitude, raster, size) {
  require('terra')
  
  long  = as.numeric(occ[, longitude])
  lat   = as.numeric(occ[, latitude])
  
  rst.col  =  terra::colFromX(raster, long)
  rst.row  =  terra::rowFromY(raster, lat)
  
  x.max = terra::xFromCol(raster, rst.col + size)
  x.min = terra::xFromCol(raster, rst.col - size)
  y.max = terra::yFromRow(raster, rst.row - size)
  y.min = terra::yFromRow(raster, rst.row + size)
  
  r = terra::rast()
  ext(r) = c(x.min, x.max, y.min, y.max)
  
  cropped = terra::crop(raster, r, snap = "out")
  
  #plot(cropped)
  #points(x=longitude,y=latitude)
  
  return(cropped)
}

# make samples
cnn_make_samples <-
  function(df,
           longitude,
           latitude,
           response,
           raster,
           size = 5) {
    
  
  data <- df
  
  predictors <- list()
  responses <- list()
  for (i in 1:nrow(data)) {
    x <- croppin.hood(occ = data[i, ],
                      longitude = longitude,
                      latitude = latitude,
                      raster = raster,
                      size = size
    ) %>%
      as.array()
    
    for (j in 1:dim(x)[3]) {
      ary <- x[,,j]
      ary <- ifelse(is.na(ary), mean(ary, na.rm = TRUE), ary) 
      x[,,j] <- ary
    }
    
    y <- data[[i, response]]
    
    predictors <- append(predictors, list(x))
    responses <- append(responses, y)
  }
  
  data_list <- list(
    "predictors" = predictors,
    "response" = responses
  )
  
  return(data_list)
}
