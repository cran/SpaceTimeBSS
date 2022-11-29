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
Version: \tab 0.2-0 \cr
Date: \tab 2022-11-29\cr
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

\author{
Christoph Muehlmann, Nikolaus Piccolotto, Claudia Cappello, Sandra De Iaco, Klaus Nordhausen

Maintainer: 
Christoph Muehlmann \email{christoph.muehlmann@tuwien.ac.at}
}


\keyword{ package }
