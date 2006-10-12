\name{best.arima}
\alias{best.arima}
\title{Fit best ARIMA model to univariate time series}
\usage{
best.arima(x, d = ndiffs(x), D = 0, max.p = 5,
             max.q = 5, max.P = 2, max.Q = 2, max.order = 5,
             stationary = FALSE, drift = TRUE, trend = FALSE,
             ic = c("aic", "bic", "aicc"))
}

\arguments{
\item{x}{a univariate time series}
\item{d}{Order of first-differencing. Will choose something sensible based on KPSS test.}
\item{D}{Order of seasonal-differencing. Default is 0. Must be manually set.}
\item{max.p}{Maximum value of p}
\item{max.q}{Maximum value of q}
\item{max.P}{Maximum value of P}
\item{max.Q}{Maximum value of Q}
\item{max.order}{Maximum value of p+q+P+Q}
\item{stationary}{If TRUE, restricts search to stationary models.}
\item{drift}{If TRUE, allows (but doesn't force) the inclusion of drift in models with d>0.}
\item{trend}{If TRUE, allows (but doesn't force) the inclusion of a linear trend in models with d=0.}
\item{ic}{Information criterion to be used in model selection.}
}
\description{Returns best ARIMA model according to either AIC or BIC value. The function conducts a search over
possible model within the order constraints provided.}
}
\details{Can be slow, especially for seasonal data. Non-seasonal differences chosen using \code{\link[tseries]{kpss.test}}.
Seasonal differences must be chosen manually.}
}
\value{Same as for \code{\link{arima}}}
}


\seealso{\code{\link{arima}}}

\author{Rob J Hyndman}
\examples{fit <- best.arima(WWWusage)
plot(forecast(fit,h=20))
}
\keyword{ts}
