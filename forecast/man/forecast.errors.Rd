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
\description{Returns range of summary measures of the error x-f. These are
(usually) out-of-sample forecast errors. Compare \code{\link{gof}} which provides a 
summary of in-sample forecast errors.
All measures are defined and discussed in Hyndman and Koehler (2006).}
}

\value{Vector giving forecast error measures.}
}

\seealso{\code{\link{gof}}}

\references{Hyndman, R.J. and Koehler, A.B. (2006) "Another look at measures of forecast accuracy". \emph{International Journal of Forecasting},
\bold{22}(4).}


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
