\name{tsdisplay}
\alias{tsdisplay}
\title{Time series display}
\usage{tsdisplay(x, plot.type = "partial", points = TRUE, ci.type =
                 "white", lag.max = round(10 * log10(length(x))),
                 na.action = na.interp, main = NULL, ylab = "", xlab =
                 "", pch = 1, cex = 0.5, ...)
}

\arguments{  \item{x}{a numeric vector or time series.}
\item{plot.type}{type of plot to include in lower right corner. Possible values are "partial", "scatter" or "spectrum".}
\item{points}{logical flag indicating whether to show the individual points or not in the time plot.}
\item{ci.type}{type of confidence limits for ACF. Possible values are as for \code{\link[stats]{acf}}.}
\item{lag.max}{the maximum lag to plot for the acf and pacf.}
\item{na.action}{how to handle missing values. Default is to use linear interpolation.}
\item{main}{Main title.}
\item{ylab}{Y-axis label}
\item{xlab}{X-axis label}
\item{pch}{Plotting character}
\item{cex}{Character size}
\item{\dots}{additional arguments to \code{\link[stats]{acf}}.}
}
\description{Plots a time series along with its acf and either its pacf, lagged scatterplot or spectrum.
}
\value{  None.
}
\references{  Makridakis, Wheelwright and Hyndman (1998) \emph{Forecasting:
methods and applications}, Wiley: New York. \url{http://www.robhyndman.info/forecasting/}
}
\author{Rob J Hyndman}
\seealso{\code{\link[stats]{plot.ts}}, \code{\link[stats]{acf}}}
\examples{tsdisplay(diff(WWWusage))
}
\keyword{ts}
