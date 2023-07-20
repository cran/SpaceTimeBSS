\name{lacov}
\alias{lacov}

\title{
Local Autocovariance Matrices
}

\description{
Computation of local autocovariance matrices for a multivariate space-time dataset based on a given set of spatio-temporal kernel functions.
}

\usage{
lacov(x, coords, time, kernel_type, kernel_parameters, 
      lags, kernel_list = NULL, center = TRUE)
}
\arguments{
  \item{x}{
either a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the space-time random field and the \code{n} rows are the observations. 
}
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} where each row represents the spatial coordinates of the corresponding observation over a 2D spatial domain.
}
  \item{time}{
a numeric vector of length \code{n} where each entry represents the temporal coordinate of the corresponding observation.
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
  \item{kernel_list}{
a list of spatio-temporal kernel matrices with dimension \code{c(n,n)}, see details. Usually computed by the function \code{\link{stkmat}}.
}
  \item{center}{
logical. If \code{TRUE} the data \code{x} is centered prior computing the local covariance matrices. Default is \code{TRUE}.
}
}

\details{
Local autocovariance matrices are defined by \deqn{ LACov(f) = 1/(n F_{f,n}) \sum_{i,j} f(s_i-s_j,t_i-t_j) (x(s_i,t_i)-\bar{x}) (x(s_j,t_j)-\bar{x})',} with \deqn{ F^2_{f,n} = 1 / n \sum_{i,j} f^2(s_i-s_j,t_i-t_j).} Here, \eqn{x(s_i,t_i)} are the \code{p} random field values at location \eqn{s_i,t_i}, \eqn{\bar{x}} is the sample mean vector, and the space-time kernel function \eqn{f} determines the locality. The following kernel functions are implemented and chosen with the argument \code{kernel_type}:
\itemize{
  \item \code{'ring'}: the spatial parameters are inner radius \eqn{r_i} and outer radius \eqn{r_o}, with \eqn{r_i < r_o}, and \eqn{r_i, r_o \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = I(r_i < d_s \le r_o)I(d_t = u)}
  \item \code{'ball'}: the spatial parameter is the radius \eqn{r}, with \eqn{r \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = I(d_s \le r)I(d_t = u)}
  \item \code{'gauss'}: Gaussian function where 95\% of the mass is inside the spatial parameter \eqn{r}, with \eqn{r \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = exp(-0.5 (\Phi^{-1}(0.95) d_s/r)^2)I(d_t = u)}
}

Above, \eqn{I()} represents the indicator function. The argument \code{kernel_type} determines the used kernel function as presented above, the argument \code{lags} provides the used temporal lags for the kernel functions (\eqn{u} in the above formulas) and the argument \code{kernel_parameters} gives the spatial parameters for the kernel function. Each of the arguments \code{kernel_type}, \code{lags} and \code{kernel_parameters} can be of length \code{K} or \code{1}. Specifically, \code{kernel_type} can be either one kernel, then each local autocovariance matrix use the same kernel type, or of length \code{K} which leads to different kernel functions for the provided kernel parameters. \code{lags} can be either one integer, then for each kernel the same temporal lag is used, or an integer vector of length  \code{K} which leads to different temporal lags. In the same fashion \code{kernel_parameters} is a vector of length \code{K} or \code{1}. If \code{kernel_type} equals \code{'ball'} or \code{'gauss'} then the corresponding entry of \code{kernel_parameters} gives the single spatial radius parameter. In contrast, if (at least one entry of) \code{kernel_type} equals \code{'ring'}, then \code{kernel_parameters} must be a list of length \code{K} (or \code{1}) where each entry is a numeric vector of length \code{2} defining the inner and outer spatial radius. See examples below.

Alternatively, a list of kernel matrices can be given directly to the function \code{lacov} through the \code{kernel_list} argument. A list with kernel matrices can be computed with the function \code{\link{stkmat}}.

}

\value{
\code{lacov} returns a list of length \code{K} where each entry is a numeric matrix of dimension \code{c(p, p)} corresponding to a local autocovariance matrix. 
}

\references{
\cite{Muehlmann, C., De Iaco, S. and Nordhausen, K. (2023), Blind Recovery of Sources for Multivariate Space-Time Environmental Data. \emph{Stochastic and Environmental Research and Risk Assessment},  \bold{37}, 1593--1613, <doi:10.1007/s00477-022-02348-2>.}
}

\seealso{
\code{\link{stkmat}}, \code{\link{stbss}}
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

# lacov with different ring kernels and same lags
lacov_r <- lacov(observed_field, coords = st_coords[, 1:2], time = st_coords[, 3], 
                 kernel_type = 'ring', 
                 kernel_parameters = list(c(0, 1), c(1, 2)), lags = 1)

# lacov with same ball kernels and different lags
lacov_b <- lacov(observed_field, coords = st_coords[, 1:2], time = st_coords[, 3], 
                 kernel_type = 'ball', kernel_parameters = 1, lags = c(1, 2, 3))

# lacov with different gauss kernels and different lags
lacov_g <- lacov(observed_field, coords = st_coords[, 1:2], time = st_coords[, 3], 
                 kernel_type = 'gauss', kernel_parameters = 1, lags = 1:3)

# lacov mixed kernels
lacov_m <- lacov(observed_field, coords = st_coords[, 1:2], time = st_coords[, 3], 
                 kernel_type = c('ball', 'ring', 'gauss'), 
                 kernel_parameters = list(1, c(1:2), 3), lags = 1:3)

# lacov with a kernel list
kernel_list <- stkmat(coords = st_coords[, 1:2], time = st_coords[, 3], 
                      kernel_type = 'ring', 
                      kernel_parameters = list(c(0, 1)), lags = 1)
lacov_k <- lacov(observed_field, kernel_list = kernel_list)
}

\keyword{ multivariate }
\keyword{ spatial }