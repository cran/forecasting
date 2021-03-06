\name{fitted.Arima}
\alias{fitted.Arima}
\title{One-step in-sample forecasts using ARIMA models}
\usage{fitted.Arima(object,...)
}

\arguments{
\item{object}{An object of class "\code{Arima}". Usually the result of a call to \code{\link[stats]{arima}}.}
\item{...}{Other arguments.}
}

\description{Returns one-step forecasts for the data used in fitting the ARIMA model.}

\value{An time series of the one-step forecasts.}

\seealso{\code{\link{forecast.Arima}}.}

\author{Rob J Hyndman}

\examples{fit <- Arima(WWWusage,c(3,1,0))
plot(WWWusage)
lines(fitted(fit),col=2)
}
\keyword{ts}
