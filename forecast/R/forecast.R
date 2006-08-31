## Generic forecast functions
## Part of forecast and demography packages

forecast <- function(...) UseMethod("forecast")

forecast.default <- function(...) forecast.ts(...)

forecast.ts <- function(x, h=ifelse(frequency(x)>1, 2*frequency(x), 10), conf=c(80,95), fan=FALSE, ...)
{
    forecast(ets(x,...),h=h,conf=conf,fan=fan)
}

print.forecast <- function(x,...)
{
#    cat(paste("Call:\n",deparse(x$call),"\n\n"))
    nconf <- length(x$conf)
    out <- ts(matrix(x$mean,ncol=1))
    attributes(out)$tsp <- attributes(x$mean)$tsp
    names <- c("Point Forecast")
    if(!is.null(x$lower) & !is.null(x$upper) & !is.null(x$conf))
    {
        x$upper <- as.matrix(x$upper)
        x$lower <- as.matrix(x$lower)
        for(i in 1:nconf)
        {
            out <- cbind(out,x$lower[,i],x$upper[,i])
            names <- c(names,paste("Lo",x$conf[i]),paste("Hi",x$conf[i]))
        }
    }
    colnames(out) <- names
    rownames(out) <- time(x$mean)
    if(any(frequency(x$x) == c(4,12)))
        print(out)
    else
        print(as.data.frame(out))
}


summary.forecast <- function(object,...)
{
    cat(paste("\nForecast method:",object$method))
#    cat(paste("\n\nCall:\n",deparse(object$call)))
    cat(paste("\n\nModel Information:\n"))
    print(object$model)
    cat("\nIn-sample error measures:\n")
    print(gof(object))
    if(is.null(object$mean))
        cat("\n No forecasts\n")
    else
    {
        cat("\nForecasts:\n")
        print(object)
    }
}

plot.forecast <- function(x, include, plot.conf=TRUE, shaded=TRUE,
        shadecols=switch(1+(length(x$conf)>1),7,length(x$conf):1),
        shadepalette=heat.colors(length(x$conf)),
        lambda=NULL, col=1,fcol=4, ylim=NULL, main=NULL, ylab="",xlab="",...)
{
    if(is.element("x",names(x))) # Assume stored as x
        data=as.ts(x$x)
    else
    {
        data=NULL
        include=0
    }

    if(missing(include))
        include <- length(data)

    if(is.null(x$lower) | is.null(x$upper) | is.null(x$conf))
        plot.conf=FALSE

    # Extract components of predict if it exists
    # This occurs for the pegels functions.
    if(is.null(main))
        main=paste("Forecasts from ",x$method,sep="")

    freq <- frequency(data)
    strt <- start(data)
    n <- length(data)
    if(plot.conf)
    {
        upper <- as.matrix(x$upper)
        lower <- as.matrix(x$lower)
    }
    pred.mean <- x$mean
    if(!is.null(lambda))  # undo Box-Cox transformation
    {
        pred.mean <- InvBoxCox(x$mean,lambda)
        if(plot.conf)
        {
            lower <- InvBoxCox(lower,lambda)
            upper <- InvBoxCox(upper,lambda)
        }
        xx <- InvBoxCox(data,lambda)
        if(lambda<0 & plot.conf)
        {
            junk <- upper
            upper <- lower
            lower <- junk
        }
    }
    else
        xx <- data
    if(is.null(ylim))
    {
        ylim <- range(c(xx[(n-include+1):n],pred.mean),na.rm=TRUE)
        if(plot.conf)
            ylim <- range(ylim,lower,upper,na.rm=TRUE)
    }
    npred <- length(pred.mean)
    plot(ts(c(xx[(n-include+1):n], rep(NA, npred)), end=tsp(xx)[2]+npred/freq, f = freq),
        xlab=xlab,ylim=ylim,ylab=ylab,main=main,col=col,...)
    xxx <- tsp(xx)[2] + (1:npred)/freq
    if(plot.conf)
    {
        idx <- rev(order(x$conf))
        nint <- length(x$conf)
        if(shaded)
        {
            if(nint>1)
                palette(shadepalette)
        }
        for(i in 1:nint)
        {
            if(shaded)
                polygon(c(xxx, rev(xxx)), c(lower[,idx[i]], rev(upper[,idx[i]])),
                            col = shadecols[i], border=FALSE)
            else
            {
                lines(xxx,lower[,idx[i]],lty=2)
                lines(xxx,upper[,idx[i]],lty=2)
            }
        }
        if(shaded)
            palette("default")
    }
    lines(ts(pred.mean, start = tsp(xx)[2]+1/freq, f = freq), lty = 1,col=fcol)
    if(plot.conf)
        invisible(list(mean=pred.mean,lower=lower,upper=upper))
    else
        invisible(list(mean=pred.mean))
}


