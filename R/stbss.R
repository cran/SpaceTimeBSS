#-------------------------------------------#
# stbss 
#-------------------------------------------#
stbss <- function(x, ...) UseMethod("stbss")

stbss.default <- function(x, coords, time, kernel_type, 
                          kernel_parameters, lags, ordered = TRUE, kernel_list = NULL, ...) {
  
  # white data
  x_w <- white_data(x)  
  
  # spatial covariance matrices
  fcall <- match.call()
  fcall$x <- x_w$x_w
  fcall$center <- FALSE  
  m <- match(methods::formalArgs(lacov), 
             names(fcall), 0)
  fcall[[1]] <- as.symbol("lacov") 
  fcall <- fcall[c(1, m)]
  lacov_list <- eval(fcall, parent.frame())
  
  # diagonalization
  lacov_d <- diag_scatters(cov_list = lacov_list, ordered = ordered, ...)
  
  # unmixing matrix
  w <- crossprod(lacov_d$u, x_w$s_inv_sqrt)
  w_inv <- crossprod(x_w$s_sqrt, lacov_d$u)
  s <- tcrossprod(x_w$x_0, w)
  colnames(s) <- paste0('IC.', 1:ncol(s))
  
  # results
  res <- list(s = s, coords = if ("coords" %in% names(fcall)) coords else NULL, 
              time = if ("time" %in% names(fcall)) time else NULL, 
              w = w, w_inv = w_inv, pevals = lacov_d$pevals,
              d = lacov_d$d, x_mu = x_w$mu, cov_inv_sqrt = x_w$s_inv_sqrt)
  return(structure(res, class = "stbss"))
}

# methods for classes of spacetime package
stbss.STFDF <- function(x, ...) {
  if (!requireNamespace('spacetime', quietly = TRUE)) {
    stop('Please install the package spacetime to use this function.')
  }
  
  if (!requireNamespace('xts', quietly = TRUE)) {
    stop('Please install the package xts to use this function.')
  }
  
  if (!requireNamespace('zoo', quietly = TRUE)) {
    stop('Please install the package zoo to use this function.')
  }
  
  x <- methods::as(x, "STIDF")
  ts_class <- xts::tclass(x@time)
  xts::tclass(x@time) <- "POSIXct"
  result <- stbss.default(x = as.matrix(x@data), coords = x@sp@coords, time = as.vector(zoo::index(x@time)), ...)
  result$s <- methods::as(spacetime::STIDF(sp = x@sp, time = x@time, data = result$s), "STFDF")
  xts::tclass(results$s@time) <- ts_class
  result$coords <- NULL
  result$time <- NULL
  return(result)
}

stbss.STSDF <- function(x, ...) {
  if (!requireNamespace('spacetime', quietly = TRUE)) {
    stop('Please install the package spacetime to use this function.')
  }
  
  if (!requireNamespace('xts', quietly = TRUE)) {
    stop('Please install the package xts to use this function.')
  }
  
  if (!requireNamespace('zoo', quietly = TRUE)) {
    stop('Please install the package xts to use this function.')
  }
  
  x <- methods::as(x, "STIDF")
  ts_class <- xts::tclass(x@time)
  xts::tclass(x@time) <- "POSIXct"
  result <- stbss.default(x = as.matrix(x@data), coords = x@sp@coords, time = as.vector(zoo::index(x@time)), ...)
  result$s <- methods::as(spacetime::STIDF(sp = x@sp, time = x@time, data = result$s), "STSDF")
  xts::tclass(results$s@time) <- ts_class
  result$coords <- NULL
  result$time <- NULL
  return(result)
}

# methods for classes of sftime package
stbss.sftime <- function(x, ...) {
  if (!requireNamespace('sftime', quietly = TRUE)) {
    stop('Please install the package sftime to use this function.')
  }
  
  if (!requireNamespace('sf', quietly = TRUE)) {
    stop('Please install the package sf to use this function.')
  }
  
  result <- stbss.default(x = as.matrix(sf::st_drop_geometry(sftime::st_drop_time(x))), 
                          coords = sf::st_coordinates(x), 
                          time = as.vector(sftime::st_time(x)), ...)
  result$s <- sftime::st_set_time(x = sf::st_set_geometry(x = data.frame(result$s), 
                                                          value = sf::st_geometry(x)), 
                                  value = sftime::st_time(x))
  result$coords <- NULL
  result$time <- NULL
  return(result)    

}

# methods for stbss class
coef.stbss <- function(object, ...) {
  object$w
}

print.stbss <- function(x, ...) {
  print.listof(x[!(names(x) %in% c('s', 'coords', 'time', 'w_inv', 'x_mu', 'cov_inv_sqrt'))], ...)
}

