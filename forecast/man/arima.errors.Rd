\name{arima.errors}
\alias{arima.errors}
\title{ARIMA errors}
\usage{arima.errors(z)}

\arguments{
\item{z}{Fitted ARIMA model from \code{\link[stats]{arima}}}
}
\description{Returns original time series after adjusting for regression variables.}
}
\value{A time series containing the errors.
}

\author{Rob J Hyndman}

\seealso{\code{\link[stats]{arima}}, \code{\link[base]{residuals}}}

\examples{ukdeaths.fit <- arima(UKDriverDeaths,c(1,0,1),c(0,1,1),xreg=Seatbelts[,"law"])
ukdeaths.errors <- arima.errors(ukdeaths.fit)
par(mfrow=c(2,1))
plot(UKDriverDeaths)
plot(ukdeaths.errors)
}
\keyword{ts}
