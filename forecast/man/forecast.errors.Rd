\name{forecasterrors}
\alias{forecasterrors}
\title{Forecast errors}
\usage{forecasterrors(f, x, test=1:length(x))
}

\arguments{
\item{f}{Vector of univariate forecasts or object of class code{"forecast"}.}
\item{x}{Vector of actual values of the same length as f.}
\item{test}{Indicator of which elements of x and f to test.}
}
\description{Returns range of summary measures of the error x-f.}

\value{Vector giving forecast error measures.
}

\author{Rob J Hyndman}
\examples{
fit1 <- rwf(EuStockMarkets[1:200,1],h=100)
fit2 <- meanf(EuStockMarkets[1:200,1],h=100)
forecasterrors(fit1,EuStockMarkets[201:300,1])
forecasterrors(fit2,EuStockMarkets[201:300,1])
plot(fit1)
lines(EuStockMarkets[1:300,1])
}
\keyword{ts}
