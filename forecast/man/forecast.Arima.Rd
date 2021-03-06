\name{forecast.Arima}
\alias{forecast.Arima}
\alias{forecast.ar}
\title{Forecasting using ARIMA models}
\usage{\method{forecast}{Arima}(object, h=ifelse(object$arma[5]>1,2*object$arma[5],10),
    level=c(80,95), fan=FALSE, xreg=NULL,...)
\method{forecast}{ar}(object, h=10, level=c(80,95), fan=FALSE, ...)
}

\arguments{
\item{object}{An object of class "\code{Arima}" or "\code{ar}". Usually the result of a call to \code{\link[stats]{arima}} or \code{\link[stats]{ar}}.}
\item{h}{Number of periods for forecasting. If \code{xreg} is used, h is ignored and the number of forecast periods is set to the number of rows of xreg.}
\item{level}{Confidence level for prediction intervals.}
\item{fan}{If TRUE, level is set to seq(50,99,by=1). This is suitable for fan plots.}
\item{xreg}{Future values of an regression variables (for class \code{Arima} objects only).}
\item{...}{Other arguments.}
}

\description{Returns forecasts and other information for univariate ARIMA models.}

\details{This function calls \code{predict.arima} or \code{predict.ar} and constructs an object of class "\code{forecast}" from the results.
}


\value{An object of class "\code{forecast}".

The function \code{summary} is used to obtain and print a summary of the
results, while the function \code{plot} produces a plot of the forecasts and prediction intervals.

The generic accessor functions \code{fitted.values} and \code{residuals} extract useful features of
the value returned by \code{forecast.Arima}.

An object of class "\code{forecast}" is a list containing at least the following elements:
\item{model}{A list containing information about the fitted model}
\item{method}{The name of the forecasting method as a character string}
\item{mean}{Point forecasts as a time series}
\item{lower}{Lower limits for prediction intervals}
\item{upper}{Upper limits for prediction intervals}
\item{level}{The confidence values associated with the prediction intervals}
\item{x}{The original time series (either \code{object} itself or the time series used to create the model stored as \code{object}).}
\item{residuals}{Residuals from the fitted model. That is x minus fitted values.}
\item{fitted}{Fitted values (one-step forecasts)}
}

\seealso{\code{\link[stats]{predict.arima}}, \code{\link[stats]{ar}}, \code{\link[stats]{arima}}, \code{\link{rwf}}.}

\author{Rob J Hyndman}

\examples{fit <- Arima(WWWusage,c(3,1,0))
plot(forecast(fit))
}
\keyword{ts}
