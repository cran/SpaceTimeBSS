\name{meteo_veneto}
\alias{meteo_veneto}
\docType{data}
\title{
Climate and Meteorological Deseasonalized Data in Veneto
}

\description{
Weekly aggregated climate and meteorological deseasonalized data in Veneto region (Italy) for a 23-year span (2000-2022).
}

\usage{data("meteo_veneto")}

\format{
  Object of class data.frame with 85248 rows and 13 variables, where 85248 consists of 1184 weekly observations times 72 spatial locations. The variables are as follows:
  \describe{
    \item{\code{x}}{x coordinates in meters (Gauss Boaga - EPSG:3003)}
    \item{\code{y}}{y coordinates in meters (Gauss Boaga - EPSG:3003)}
    \item{\code{sp.ID}}{code for the spatial location}
    \item{\code{timeIndex}}{code for the temporal observation, from 1 to 1184}
    \item{\code{deseas_ET0}}{evapotranspiration levels (mm)}
    \item{\code{deseas_rad}}{solar radiation (MJ/m^2)}
    \item{\code{deseas_tmax}}{maximum temperature (degrees C)}
    \item{\code{deseas_taver}}{average temperature (degrees C)}
    \item{\code{deseas_tmin}}{minimum temperature (degrees C)}
    \item{\code{deseas_hmax}}{maximum humidity (\%)}
    \item{\code{deseas_hmin}}{minimum humidity (\%)}
    \item{\code{deseas_wind}}{wind velocity (m/s)}
    \item{\code{deseas_log_prec}}{log of precipitation (mm)}
  }
}

\details{
The evapotranspiration levels were estimated by ARPA Veneto according to the Hargreaves equation. The data have been obtained by removing the annual periodicity from the raw data and then computing weekly averages.
}

\source{
The raw data can be downloaded from the Environmental Protection Agency of Veneto Region (ARPA Veneto) website.
}

\keyword{datasets}
