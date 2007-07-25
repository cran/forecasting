## Measures of forecast accuracy
## Forecasts in f. This may be a numerical vector or the output from arima or pegels or derivatives.
## Actual values in x
## test enables a subset of x and f to be tested.
forecasterrors <- function(f,x,test=1:length(x))
{
    if(is.list(f))
    {
        if(is.element("mean",names(f)))
            f = f$mean
        else
            stop("Unknown list structure")
    }
    n <- length(x)
    error <- (x-f[1:n])[test]
    me <- mean(error)
    mse <- mean(error^2)
    mae <- mean(abs(error))
    mape <- mean(100*abs(error/x[test]))
    mpe <-  mean(100*error/x[test])
    fpe <- (c(f[2:n])/c(x[1:(n-1)]) - 1)[test-1]
    ape <- (c(x[2:n])/c(x[1:(n-1)]) - 1)[test-1]
    theil <- sqrt(sum((fpe - ape)^2)/sum(ape^2))
    r1 <- acf(error,plot=FALSE,lag.max=2)$acf[2,1,1]
    junk <- c(me,mae,mse,mpe,mape,r1,theil)
    names(junk) <- c("ME","MAE","MSE","MPE","MAPE","ACF1","Theil's U")
    junk
}


gof <- function(object)
{
    res <- object$x-object$fitted
    # Don't use object$resid as this may contain multiplicative errors.
    pe <- res/object$x * 100 # Percentage error
    scale <- mean(abs(diff(object$x)),na.rm=TRUE)
    out <- c(mean(res,na.rm=TRUE), sqrt(mean(res^2,na.rm=TRUE)), mean(abs(res),na.rm=TRUE), mean(pe,na.rm=TRUE), mean(abs(pe),na.rm=TRUE),
        mean(abs(res/scale),na.rm=TRUE))
    names(out) <- c("ME","RMSE","MAE","MPE","MAPE","MASE")
    return(out)
}
