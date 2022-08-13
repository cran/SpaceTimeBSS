\name{stbss}
\alias{stbss}
\alias{stbss.default}
\alias{stbss.STFDF}
\alias{stbss.STSDF}
\alias{stbss.sftime}

\title{
Space Time Blind Source Separation
}

\description{
For a given multivariate space-time dataset, \code{stbss} estimates the realization of spatio-temporally uncorrelated random fields through a linear transformation which is defined by a so-called mixing matrix and a location vector. This is done assuming a spatio-temporal blind source separation model and simultaneously/jointly diagonalizing the sample covariance matrix and one/many local autocovariance matrices. 
}

\usage{
stbss(x, ...)

\method{stbss}{default}(x, coords, time, kernel_type, 
      kernel_parameters, lags, ordered = TRUE, kernel_list = NULL, \dots)
\method{stbss}{STFDF}(x, \dots)
\method{stbss}{STSDF}(x, \dots)
\method{stbss}{sftime}(x, \dots)
}

\arguments{
  \item{x}{
either a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the space-time random field and the \code{n} rows are the observations, an object of class \code{\link[spacetime]{STFDF}}, an object of class \code{\link[spacetime]{STSDF}} or an object of class \code{\link[sftime]{st_sftime}}. 
}
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} where each row represents the coordinates of a point in the spatial domain over a 2D spatial domain. Only needed if \code{x} is a matrix and the argument \code{kernel_list} is \code{NULL}.
}
  \item{time}{
a numeric vector of length \code{n} where each entry represents the time of a point in the temporal domain. Only needed if \code{x} is a matrix and the argument \code{kernel_list} is \code{NULL}.
}
  \item{kernel_type}{
either a string or a string vector of length \code{K} (or \code{1}) indicating which spatio-temporal kernel function to use. Implemented choices are \code{'ring'}, \code{'ball'} or \code{'gauss'}.
}
  \item{kernel_parameters}{
a numeric vector of length \code{K} (or \code{1}) for the \code{'ball'} and \code{'gauss'} kernel function or a list of length \code{K} (or \code{1}) for the \code{'ring'} kernel, see details.
}
  \item{lags}{
an integer vector of length \code{K} (or \code{1}) that provides the temporal lags for the spatio-temporal kernel functions, see details.
}
  \item{ordered}{
logical. If \code{TRUE} the entries of the latent field are ordered by the sum of squared (pseudo-)eigenvalues of the diagonalized local covariance matrix/matrices. Default is \code{TRUE}.
}
  \item{kernel_list}{
a list of spatio-temporal kernel matrices with dimension \code{c(n,n)}, see details. Usually computed by the function \code{\link{stkmat}}.
}
  \item{\dots}{
further arguments for the fast real joint diagonalization algorithm that jointly diagonalizes the local covariance matrices. See details and \code{\link[JADE]{frjd}}.
}
}

\details{
It is assumed that the p-variate space-time random field \eqn{x(s,t)} is formed by \deqn{x(s,t) = A z(s,t) + b,} where \eqn{z(s,t)} is the latent p-variate space-time random field, \eqn{A} and \eqn{b} are the mixing matrix and a location vector and \eqn{s} and \eqn{t} are the space and time coordinates. Furthermore, it is assumed that \eqn{z(s,t)} is white and consists of space-time uncorrelated components. The goal is to reverse the linear form by estimating an unmixing matrix and the location vector. This is done by simultaneously/jointly diagonalizing local autocovariance matrices which are defined by \deqn{ LACov(f) = 1/(n F_{f,n}) \sum_{i,j} f(s_i-s_j,t_i-t_j) (x(s_i,t_i)-\bar{x}) (x(s_j,t_j)-\bar{x})',} with \deqn{ F^2_{f,n} = 1 / n \sum_{i,j} f^2(s_i-s_j,t_i-t_j).} Here, \eqn{x(s_i,t_i)} are the \code{p} random field values at location \eqn{s_i,t_i}, \eqn{\bar{x}} is the sample mean vector, and the space-time kernel function \eqn{f} determines the locality. The following kernel functions are implemented and chosen with the argument \code{kernel_type}:
\itemize{
  \item \code{'ring'}: the spatial parameters are inner radius \eqn{r_i} and outer radius \eqn{r_o}, with \eqn{r_i < r_o}, and \eqn{r_i, r_o \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = I(r_i < d_s \le r_o)I(d_t = u)}
  \item \code{'ball'}: the spatial parameter is the radius \eqn{r}, with \eqn{r \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = I(d_s \le r)I(d_t = u)}
  \item \code{'gauss'}: Gaussian function where 95\% of the mass is inside the spatial parameter \eqn{r}, with \eqn{r \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = exp(-0.5 (\Phi^{-1}(0.95) d_s/r)^2)I(d_t = u)}
}

Above, \eqn{I()} represents the indicator function. The argument \code{kernel_type} determines the used kernel function as presented above, the argument \code{lags} provides the used temporal lags for the kernel functions (\eqn{u} in the above formulas) and the argument \code{kernel_parameters} gives the spatial parameters for the kernel function. Each of the arguments \code{kernel_type}, \code{lags} and \code{kernel_parameters} can be of length \code{K} or \code{1}. Specifically, \code{kernel_type} can be either one kernel, then each local autocovariance matrix use the same kernel type, or of length \code{K} which leads to different kernel functions for the provided kernel parameters. \code{lags} can be either one integer, then for each kernel the same temporal lag is used, or an integer vector of length  \code{K} which leads to different temporal lags. In the same fashion \code{kernel_parameters} is a vector of length \code{K} or \code{1}. If \code{kernel_type} equals \code{'ball'} or \code{'gauss'} then the corresponding entry of \code{kernel_parameters} gives the single spatial radius parameter. In contrast, if (at least one entry of) \code{kernel_type} equals \code{'ring'}, then \code{kernel_parameters} must be a list of length \code{K} (or \code{1}) where each entry is a numeric vector of length \code{2} defining the inner and outer spatial radius.

Internally, \code{stbss} calls \code{\link{stkmat}} to compute a list of \code{c(n,n)} kernel matrices based on the parameters given, where each entry of those matrices corresponds to \eqn{f(s_i-s_j,t_i-t_j)}. Alternatively, such a list of kernel matrices can be given directly to the function \code{stbss} through the \code{kernel_list} argument. This is useful when \code{stbss} is called numerous times with the same coordinates/kernel functions as the computation of the kernel matrices is then only done once prior the actual \code{stbss} calls. For details see also \code{\link{lacov}}.

If more than one local autocovariance matrix is used \code{stbss} jointly diagonalizes these matrices with the function \code{\link[JADE]{frjd}}. \code{\dots} provides arguments for \code{frjd}, useful arguments might be:
\itemize{
  \item \code{eps}: tolerance for convergence.
  \item \code{maxiter}: maximum number of iterations.
}

}

\value{
\code{stbss} returns a list of class \code{'stbss'} with the following entries: 
  \item{s}{
  object of \code{class(x)} containing the estimated source space-time values.
}  
  \item{coords}{
  coordinates of the observations. Is \code{NULL} if \code{class(x)} is not a matrix or if \code{kernel_list} is provided at the \code{stbss} call.
}
  \item{time}{
  time of the observations. Is \code{NULL} if \code{kernel_list} is provided or if \code{class(x)} is not a matrix at the \code{stbss} call.
}
  \item{w}{
  estimated unmixing matrix.
}
  \item{w_inv}{
  inverse of the estimated unmixing matrix.
}
  \item{pevals}{
  (pseudo-)eigenvalues for each latent field entry.
}
  \item{d}{
  matrix of stacked (jointly) diagonalized local autocovariance matrices with dimension \code{c(length(kernel_parameters)*p,p)}.
}
  \item{x_mu}{
  columnmeans of \code{x}.
}
  \item{cov_inv_sqrt}{
  square root of the inverse sample covariance matrix of \code{x}.
}

}

\seealso{
\code{\link{stkmat}}, \code{\link[JADE]{frjd}}
}

\examples{
# space and time coordinates
n_t <- 50
n_sp <- 10
st_coords <- as.matrix(expand.grid(1:n_sp, 1:n_sp, 1:n_t))

# simulate three latent white noise fields
field_1 <- rnorm(nrow(st_coords))
field_2 <- rnorm(nrow(st_coords))
field_3 <- rnorm(nrow(st_coords))
  
# compute the observed field
latent_field <- cbind(field_1, field_2, field_3)
mixing_matrix <- matrix(rnorm(9), 3, 3)
observed_field <- latent_field %*% t(mixing_matrix)

# apply stbss with lag 1 and a ring kernel
stbss_res <- stbss(observed_field, coords = st_coords[, 1:2], time = st_coords[, 3], 
                   kernel_type = 'ring', kernel_parameters = list(c(0, 1)), lags = 1)

# print object
print(stbss_res)

# unmixing matrix
w_unmix <- coef(stbss_res)

# apply the same stbss with a kernel list
kernel_list <- stkmat(coords = st_coords[, 1:2], time = st_coords[, 3], 
                      kernel_type = 'ring', kernel_parameters = list(c(0, 1)), lags = 1)
stbss_res_k <- stbss(observed_field, kernel_list = kernel_list)

# apply stbss with three ball kernels 
stbss_res_b <- stbss(observed_field, coords = st_coords[, 1:2], time = st_coords[, 3], 
                     kernel_type = 'ball', kernel_parameters = 1:3, lags = 1:3)


}

\keyword{ multivariate }
\keyword{ spatial }
