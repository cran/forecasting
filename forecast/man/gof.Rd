\name{gof}
\alias{gof}
\title{Goodness-of-fit for forecast model}
\usage{gof(object)
}

\arguments{
\item{object}{An object of class \code{"forecast"}}
}

\description{Returns range of summary measures of the in-sample errors. All measures are defined and discussed in Hyndman and Koehler (2006).}

\value{Vector giving in-sample error measures.
}

\references{Hyndman, R.J. and Koehler, A.B. (2006) "Another look at measures of forecast accuracy". \emph{International Journal of Forecasting},
\bold{22}(4).}

\author{Rob J Hyndman}
\examples{
fit1 <- rwf(Nile)
fit2 <- meanf(Nile)
gof(fit1)
gof(fit2)
}
\keyword{ts}
