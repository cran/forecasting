\name{na.interp}
\alias{na.interp}
\title{Interpolate missing values in a time series}
\usage{na.interp(x)}

\arguments{  
\item{x}{time series}
}
\description{Uses linear interpolation to replace missing values.}

\details{A more general and flexible approach is available using \code{\link[zoo]{na.approx}}.}

\value{Time series}

\seealso{\code{\link[zoo]{na.approx}}}


\author{Rob J Hyndman}
\examples{
data(gold)
plot(na.interp(gold))
}
\keyword{ts}
