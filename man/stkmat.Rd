\name{stkmat}
\alias{stkmat}

\title{
Spatio-Temporal Kernel Matrices
}

\description{
Computation of spatio-temporal kernel matrices for given kernel functions.
}

\usage{
stkmat(coords, time, kernel_type, kernel_parameters, lags)
}

\arguments{
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} where each row represents the spatial coordinates of the corresponding observation over a 2D spatial domain.
}
  \item{time}{
an integer vector of length \code{n} where each entry represents the temporal coordinate of the corresponding observation.
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
}

\details{
The following kernel functions are implemented and chosen with the argument \code{kernel_type}:
\itemize{
  \item \code{'ring'}: the spatial parameters are inner radius \eqn{r_i} and outer radius \eqn{r_o}, with \eqn{r_i < r_o}, and \eqn{r_i, r_o \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = I(r_i < d_s \le r_o)I(d_t = u)}
  \item \code{'ball'}: the spatial parameter is the radius \eqn{r}, with \eqn{r \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = I(d_s \le r)I(d_t = u)}
  \item \code{'gauss'}: Gaussian function where 95\% of the mass is inside the spatial parameter \eqn{r}, with \eqn{r \ge 0}, the temporal parameter is the temporal lag \eqn{u}: \deqn{f(d_s,d_t) = exp(-0.5 (\Phi^{-1}(0.95) d_s/r)^2)I(d_t = u)}
}

Above, \eqn{I()} represents the indicator function. The argument \code{kernel_type} determines the used kernel function as presented above, the argument \code{lags} provides the used temporal lags for the kernel functions (\eqn{u} in the above formulas) and the argument \code{kernel_parameters} gives the spatial parameters for the kernel function. Each of the arguments \code{kernel_type}, \code{lags} and \code{kernel_parameters} can be of length \code{K} or \code{1}. Specifically, \code{kernel_type} can be either one kernel, then each local autocovariance matrix use the same kernel type, or of length \code{K} which leads to different kernel functions for the provided kernel parameters. \code{lags} can be either one integer, then for each kernel the same temporal lag is used, or an integer vector of length  \code{K} which leads to different temporal lags. In the same fashion \code{kernel_parameters} is a vector of length \code{K} or \code{1}. If \code{kernel_type} equals \code{'ball'} or \code{'gauss'} then the corresponding entry of \code{kernel_parameters} gives the single spatial radius parameter. In contrast, if (at least one entry of) \code{kernel_type} equals \code{'ring'}, then \code{kernel_parameters} must be a list of length \code{K} (or \code{1}) where each entry is a numeric vector of length \code{2} defining the inner and outer spatial radius. See examples below.

The output of this function can be used with the function \code{\link{stbss}} to avoid unnecessary computation of kernel matrices when \code{\link{stbss}} is called multiple times with the same coordinate/kernel function setting. Additionally, the output can be used with the function \code{\link{lacov}}.
}

\value{
\code{stkmat} returns a list of length \code{K} containing numeric matrices of dimension \code{c(n,n)} corresponding to the spatio-temporal kernel matrices.
}

\references{
\cite{Muehlmann, C., De Iaco, S. and Nordhausen, K. (2023), Blind Recovery of Sources for Multivariate Space-Time Environmental Data. \emph{Stochastic and Environmental Research and Risk Assessment},  \bold{37}, 1593--1613, <doi:10.1007/s00477-022-02348-2>.}
}

\seealso{
\code{\link{stbss}}, \code{\link{lacov}}
}

\examples{

# space and time coordinates
n_t <- 50
n_sp <- 10
coords <- runif(n_sp ^ 2 * 2) * n_sp
dim(coords) <- c(n_sp ^ 2, 2)
time <- 1:n_t

st_coords <- as.matrix(expand.grid(1:nrow(coords), 1:length(time)))
st_coords <- cbind(coords[st_coords[, 1], ], time[st_coords[, 2]])

# different ring kernels and same lags
stkmat_r <- stkmat(coords = st_coords[, 1:2], time = st_coords[, 3], 
                   kernel_type = 'ring', 
                   kernel_parameters = list(c(0, 1), c(1, 2)), lags = c(1, 1))

# same ball kernels and different lags
stkmat_b <- stkmat(coords = st_coords[, 1:2], time = st_coords[, 3], 
                   kernel_type = 'ball', kernel_parameters = 1:3, lags = c(1, 2, 3))

# different gauss kernels and different lags
stkmat_g <- stkmat(coords = st_coords[, 1:2], time = st_coords[, 3], 
                   kernel_type = 'gauss', kernel_parameters = 1:3, lags = 1:3)

# mixed kernels
stkmat_m <- stkmat(coords = st_coords[, 1:2], time = st_coords[, 3], 
                   kernel_type = c('ball', 'ring', 'gauss'), 
                   kernel_parameters = list(1, c(1:2), 3), lags = 1:3)
}

\keyword{ array }

