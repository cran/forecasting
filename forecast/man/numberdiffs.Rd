\name{ndiffs}
\alias{ndiffs}
\title{Number of differences}
\usage{ndiffs(x,alpha=0.05)
}
\arguments{
\item{x}{A univariate time series}
\item{alpha}{Level of the test}
}

\description{Uses a KPSS test for the null hypothesis that \code{x} has a stationary root against a unit-root alternative.
    Returns the least number of differences required to pass the test at the level \code{alpha}.}

\value{An integer.}


\seealso{\code{\link[tseries]{kpss.test}}}

\author{Rob J Hyndman}

\examples{ndiffs(WWWusage)}

\keyword{ts}
