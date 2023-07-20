\name{SpaceTimeBSS-package}
\alias{SpaceTimeBSS-package}
\docType{package}
\title{
Blind Source Separation for Multivariate Spatio-Temporal Data
}

\description{
Simultaneous/joint diagonalization of local autocovariance matrices to estimate spatio-temporally uncorrelated random fields. 
}

\details{
\tabular{ll}{
Package: \tab SpaceTimeBSS\cr
Type: \tab Package\cr
Version: \tab 0.3-0 \cr
Date: \tab 2023-07-20\cr
License: \tab GPL (>= 2)\cr
}

Solving the second order blind source separation problem for multivariate space-time random fields. The random fields can be irregular in space but must be regular in time. The main function of this package is:
\itemize{
\item{\code{\link{stbss}}} {This function computes local autocovariance matrices. The considered temporal lags are integer numbers and the spatial lags are defined by spatial kernel functions. Then, these local autocovariance matrices and the sample covariance are simultaneously/jointly diagonalized.}
}

The package also contains a 9-variate dataset of deseasonalized weekly climate and meteorological measurements from the Italian Veneto region between 2000 and 2022 \code{\link{meteo_veneto}}.

Joint diagonalization is computed with the \code{\link[JADE]{frjd}} (fast real joint diagonalization) algorithm from the package \code{\link[JADE:JADE-package]{JADE}}.

The available finite realizations of the space time random fields can be defined by matrices or an object of classes \code{\link[spacetime]{STFDF}}, \code{\link[spacetime]{STSDF}} or \code{\link[sftime]{st_sftime}}.

}

\references{
\cite{Muehlmann, C., De Iaco, S. and Nordhausen, K. (2023), Blind Recovery of Sources for Multivariate Space-Time Environmental Data. \emph{Stochastic and Environmental Research and Risk Assessment},  \bold{37}, 1593--1613, <doi:10.1007/s00477-022-02348-2>.}
}

\author{
Christoph Muehlmann, Nikolaus Piccolotto, Claudia Cappello, Sandra De Iaco, Klaus Nordhausen

Maintainer: 
Klaus Nordhausen \email{klausnordhausenR@gmail.com}
}


\keyword{ package }
