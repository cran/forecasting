\name{dm.test}
\alias{dm.test}
\title{Diebold-Mariano test for predictive accuracy}
\description{The Diebold-Mariano test 
compares the forecast accuracy of two forecast methods. The null hypothesis is that they have the same
forecast accuracy.}

\usage{
dm.test(e1, e2, alternative = c("two.sided", "less", "greater"), h = 1, power = 2)
}
\arguments{
  \item{e1}{Forecast errors from method 1.}
  \item{e2}{Forecast errors from method 2.}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}. 
    You can specify just the initial letter.}
  \item{h}{The forecast horizon used in calculating \code{e1} and \code{e2}.}
  \item{power}{The power used in the loss function. Usually 1 or 2.}
}

\value{
A list with class \code{"htest"} containing the following components: 
\item{statistic}{the value of the DM-statistic.}
\item{parameter}{the forecast horizon and loss function power used in the test.} 
\item{alternative}{a character string describing the alternative hypothesis.}
\item{p.value}{the p-value for the test.}
\item{method}{a character string with the value "Diebold-Mariano Test".} 
\item{data.name}{a character vector giving the names of the two error series.}
}

\references{
Diebold, F.X. and Mariano, R.S. (1995) Comparing predictive accuracy. \emph{Journal of Business
and Economic Statistics}, \bold{13}, 253-263.}

\author{George Athanasopoulos and Rob Hyndman}

\examples{
# Test on in-sample one-step forecasts
f1 <- ets(WWWusage)
f2 <- auto.arima(WWWusage)
accuracy(f1)
accuracy(f2)
dm.test(residuals(f1),residuals(f2),h=1)

# Test on out-of-sample one-step forecasts
f1 <- ets(WWWusage[1:80])
f2 <- auto.arima(WWWusage[1:80])
f1.out <- ets(WWWusage[81:100],model=f1)
f2.out <- Arima(WWWusage[81:100],model=f2)
accuracy(f1.out)
accuracy(f2.out)
dm.test(residuals(f1.out),residuals(f2.out),h=1)
}

\keyword{htest}
\keyword{ts}
