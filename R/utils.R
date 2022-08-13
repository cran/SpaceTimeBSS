#------------------------------------------------#
#   data whitening
#------------------------------------------------#
white_data <- function(x) {
  mu <- colMeans(x)
  x_0 <- sweep(x, MARGIN = 2, STATS = mu, FUN = '-')
  
  s <- crossprod(x_0) / (nrow(x) - 1)
  
  s_evd <- eigen(s, symmetric = TRUE)
  s_inv_sqrt <- s_evd$vectors %*% tcrossprod(diag(1 / sqrt(s_evd$values)), s_evd$vectors)
  s_sqrt <- s_evd$vectors %*% tcrossprod(diag(sqrt(s_evd$values)), s_evd$vectors)
  
  x_w <- tcrossprod(x_0, s_inv_sqrt)
  colnames(x_w) <- colnames(x_0)
  
  return(list(mu = mu, x_0 = x_0, x_w = x_w, s = s, s_inv_sqrt = s_inv_sqrt, s_sqrt = s_sqrt))
}

#------------------------------------------------#
#   space time kernel function computation
#------------------------------------------------#
stkmat <- function(coords, time, kernel_type, kernel_parameters, lags) {
  
  kernel_list <- list()  
  
  k <- max(c(length(kernel_type), length(kernel_parameters), length(lags)))
  
  if (length(kernel_type) == 1) {
    kernel_type = rep(kernel_type, k)
  }
  
  if (length(kernel_parameters) == 1) {
    kernel_parameters = rep(kernel_parameters, k)
  }
  
  if (length(lags) == 1) {
    lags = rep(lags, k)
  }
  
  for (idx in 1:k) {
    if (kernel_type[idx] == 'ball') {
      k_mat <- if (kernel_parameters[[idx]] >= 0) stkmat_ball(coords = coords, time = time,
                                                              tau = lags[idx], h = kernel_parameters[[idx]])
      else stop('Radius must be zero or positive.')  
    } else if (kernel_type[idx] == 'ring') {
      
      k_mat <- if (kernel_parameters[[idx]][1] >= kernel_parameters[[idx]][2]) stop('Inner radius must be smaller than outer radius.') 
      else stkmat_ring(coords = coords, time = time, tau = lags[idx],
                       h1 = kernel_parameters[[idx]][1], h2 = kernel_parameters[[idx]][2]) 
    } else if (kernel_type[idx] == 'gauss') {
      k_mat <- if (kernel_parameters[[idx]] >= 0) stkmat_gauss(coords = coords, time = time,
                                                               tau = lags[idx], h = kernel_parameters[[idx]])
      else stop('Radius must be zero or positive.')
    } else {
      stop('Invalid input!')
    }
    
    kernel_list <- c(kernel_list, list(k_mat))
  }
  
  return(kernel_list)
}

#------------------------------------------------#
#  local autocovariance matrix computation
#------------------------------------------------#
lacov <- function(x, coords, time, kernel_type, kernel_parameters, 
                  lags, kernel_list = NULL, center = TRUE) {
  
  if (center) {
    x <- scale(x, center = TRUE, scale = FALSE)
  }
  
  if (!missing(coords) && !missing(time) && !missing(lags) &&
      !missing(kernel_parameters) && is.vector(kernel_parameters)) {
    k <- max(c(length(kernel_type), length(kernel_parameters), length(lags)))
    
    if (length(kernel_type) == 1) {
      kernel_type = rep(kernel_type, k)
    }
    
    if (length(kernel_parameters) == 1) {
      kernel_parameters = rep(kernel_parameters, k)
    }
    
    if (length(lags) == 1) {
      lags = rep(lags, k)
    }
    lacov_list <- list()     
    
    for (idx in 1:k) {
      if (kernel_type[idx] == 'ball') {
        lacov_mat <- if (kernel_parameters[[idx]] >= 0) lacov_ball(coords = coords, time = time, x = x,
                                                                   h = kernel_parameters[[idx]], tau = lags[idx])
        else stop('Radius must be zero or positive.')  
      } else if (kernel_type[idx] == 'ring') {
        lacov_mat <- if (kernel_parameters[[idx]][1] >= kernel_parameters[[idx]][2]) stop('Inner radius must be smaller than outer radius.') 
        else lacov_ring(coords = coords, time = time, x = x, tau = lags[idx],
                        h1 = kernel_parameters[[idx]][1], h2 = kernel_parameters[[idx]][2]) 
      } else if (kernel_type[idx] == 'gauss') {
        lacov_mat <- if (kernel_parameters[[idx]] >= 0) lacov_gauss(coords = coords, time = time, x = x,
                                                                    h = kernel_parameters[[idx]], tau = lags[idx]) 
        else stop('Radius must be zero or positive.')  
      } else {
        stop('Invalid input!')
      }
      
      lacov_list <- c(lacov_list, list(lacov_mat))
    }      
  } else if (!is.null(kernel_list) && is.list(kernel_list)) {
    if (missing(coords)) {
      lacov_list <- lapply(kernel_list, function(k_mat) lacov_kmat(x = x, k = k_mat))
    }
  } else {
    stop('Invalid input for kernels.')
  }
  
  lacov_list <- lapply(lacov_list, function(x) (x + t(x)) / 2)
  
  return(lacov_list)
}

#------------------------------------------------#
#  function for scatter diagonalization
#------------------------------------------------#
diag_scatters <- function(cov_list, ordered, ...) {
  k <- length(cov_list)
  if (k == 1) {
    cov_evd <- eigen(cov_list[[1]], symmetric = TRUE)
    u <- cov_evd$vectors
    d <- diag(cov_evd$values)
  } else {
    jade <- JADE::frjd(do.call(rbind, cov_list), ...)
    u <- jade$V
    d <- jade$D
  }
  
  p <- ncol(d)    
  diags_mat <- matrix(0, nrow = k, ncol = p)
  for (idx in 1:k) {
    diags_mat[idx, ] <- diag(d[(1:p) + (idx - 1) * p, ])
  }
  pevals <- colSums(diags_mat ^ 2)
  if (ordered) {
    diag_order <- order(pevals, decreasing = TRUE)
    u <- u[, diag_order]
    for (idx in 1:k) {
      d[(1:p) + (idx - 1) * p, ] <- d[(1:p) + (idx - 1) * p, ][diag_order, diag_order]
    }
  } 
  
  return(list(u = u, d = d, pevals = pevals))
}

