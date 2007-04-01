# Best.arima kept for backwards compatibility
best.arima <- function(...)
{
    warning("best.arima() will be deprecated in a later version. Please use auto.arima() instead.\n")
    return(auto.arima(...))
}

auto.arima <- function(x, d=NA, D=NA, max.p=5, max.q=5,
    max.P=2, max.Q=2, max.order=5,
    start.p=2, start.q=2, start.P=1, start.Q=1,
    stationary=FALSE, ic=c("aic","aicc","bic"),
    stepwise=TRUE, trace=FALSE)
{
    if(!stepwise)
        return(search.arima(x,d,D,max.p,max.q,max.P,max.Q,max.order,stationary,ic,trace))

    ic <- match.arg(ic)
    oldwarn <- options()$warn
    options(warn=-1)
    on.exit(options(warn=oldwarn))

    # Choose order of differencing
    m <- frequency(x)
    if(stationary)
        d <- D <- 0
    if(m == 1)
        D <- max.P <- max.Q <- 0
    else if(is.na(D))
        D <- forecast:::nsdiffs(x)
    if(D > 0)
        dx <- diff(x,D)
    else
        dx <- x
    if(is.na(d))
        d <- forecast:::ndiffs(dx)

    if(m > 1)
    {
        max.p <- min(max.p, m-1)
        max.q <- min(max.q, m-1)
    }

    # Starting model
    p <- start.p <- min(start.p,max.p)
    q <- start.q <- min(start.q,max.q)
    P <- start.P <- min(start.P,max.P)
    Q <- start.Q <- min(start.Q,max.Q)
    constant <- (d+D <= 1)
    results <- matrix(NA,nrow=100,ncol=8)

    bestfit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace)
    results[1,] <- c(p,d,q,P,D,Q,constant,bestfit$ic)
    # Null model
    fit <- myarima(x,order=c(0,d,0),seasonal=c(0,D,0),constant=constant,ic,trace)
    results[2,] <- c(0,d,0,0,D,0,constant,fit$ic)
    if(fit$ic < bestfit$ic)
    {
        bestfit <- fit
        p <- q <- P <- Q <- 0
    }
    # Basic AR model
    fit <- myarima(x,order=c(1,d,0),seasonal=c(m>1,D,0),constant=constant,ic,trace)
    results[3,] <- c(1,d,0,m>1,D,0,constant,fit$ic)
    if(fit$ic < bestfit$ic)
    {
        bestfit <- fit
        p <- 1
        P <- m>1
        q <- Q <- 0
    }
    # Basic MA model
    fit <- myarima(x,order=c(0,d,1),seasonal=c(0,D,m>1),constant=constant,ic,trace)
    results[4,] <- c(0,d,1,0,D,m>1,constant,fit$ic)
    if(fit$ic < bestfit$ic)
    {
        bestfit <- fit
        p <- P <- 0
        Q <- m>1
        q <- 1
    }

    startk <- 0
    k <- 4
    while(startk < k & k < 94)
    {
        startk <- k
        if(P > 0 & newmodel(p,d,q,P-1,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P-1,D,Q),constant=constant,ic,trace)
            results[k,] <- c(p,d,q,P-1,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                P <- (P-1)
            }
        }
        if(P < max.P & newmodel(p,d,q,P+1,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P+1,D,Q),constant=constant,ic,trace)
            results[k,] <- c(p,d,q,P+1,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                P <- (P+1)
            }
        }
        if(Q > 0 & newmodel(p,d,q,P,D,Q-1,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q-1),constant=constant,ic,trace)
            results[k,] <- c(p,d,q,P,D,Q-1,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                Q <- (Q-1)
            }
        }
        if(Q < max.Q & newmodel(p,d,q,P,D,Q+1,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q+1),constant=constant,ic,trace)
            results[k,] <- c(p,d,q,P,D,Q+1,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                Q <- (Q+1)
            }
        }
        if(Q > 0 & P > 0 & newmodel(p,d,q,P-1,D,Q-1,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P-1,D,Q-1),constant=constant,ic,trace)
            results[k,] <- c(p,d,q,P-1,D,Q-1,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                Q <- (Q-1)
                P <- (P-1)
            }
        }
        if(Q < max.Q & P < max.P & newmodel(p,d,q,P+1,D,Q+1,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P+1,D,Q+1),constant=constant,ic,trace)
            results[k,] <- c(p,d,q,P+1,D,Q+1,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                Q <- (Q+1)
                P <- (P+1)
            }
        }

        if(p > 0 & newmodel(p-1,d,q,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p-1,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace)
            results[k,] <- c(p-1,d,q,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                p <- (p-1)
            }
        }
        if(p < max.p & newmodel(p+1,d,q,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p+1,d,q),seasonal=c(P,D,Q),constant=constant,ic,trace)
            results[k,] <- c(p+1,d,q,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                p <- (p+1)
            }
        }
        if(q > 0 & newmodel(p,d,q-1,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q-1),seasonal=c(P,D,Q),constant=constant,ic,trace)
            results[k,] <- c(p,d,q-1,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                q <- (q-1)
            }
        }
        if(q < max.q & newmodel(p,d,q+1,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q+1),seasonal=c(P,D,Q),constant=constant,ic,trace)
            results[k,] <- c(p,d,q+1,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                q <- (q+1)
            }
        }
        if(q > 0 & p > 0 & newmodel(p-1,d,q-1,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p-1,d,q-1),seasonal=c(P,D,Q),constant=constant,ic,trace)
            results[k,] <- c(p-1,d,q-1,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                q <- (q-1)
                p <- (p-1)
            }
        }
        if(q < max.q & p < max.p & newmodel(p+1,d,q+1,P,D,Q,constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p+1,d,q+1),seasonal=c(P,D,Q),constant=constant,ic,trace)
            results[k,] <- c(p+1,d,q+1,P,D,Q,constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                q <- (q+1)
                p <- (p+1)
            }
        }
        if(newmodel(p,d,q,P,D,Q,!constant,results[1:k,]))
        {
            k <- k + 1
            fit <- myarima(x,order=c(p,d,q),seasonal=c(P,D,Q),constant=!constant,ic,trace)
            results[k,] <- c(p,d,q,P,D,Q,!constant,fit$ic)
            if(fit$ic < bestfit$ic)
            {
                bestfit <- fit
                constant <- !constant
            }
        }
    }

    # Return best fit
    bestfit$x <- x
    bestfit$series <- deparse(substitute(x))
    bestfit$ic=NULL

    if(trace)
        cat("\n\n Best model:",arima.string(bestfit),"\n\n")

    return(bestfit)
}


# Calls arima from stats package and adds data to the returned object
# Also allows refitting to new data
# and drift terms to be included.
myarima <- function(x, order = c(0, 0, 0), seasonal = c(0, 0, 0), constant=TRUE, ic="aic", trace=FALSE)
{
    n <- length(x)
    m <- frequency(x)
    use.season <- (sum(seasonal)>0) & m>0
    diffs <- order[2]+seasonal[2]
    if(diffs==1 & constant)
    {
        xreg <- 1:length(x)
        if(use.season)
            fit <- try(stats:::arima(x=x,order=order,seasonal=list(order=seasonal,period=m),xreg=xreg),silent=TRUE)
        else
            fit <- try(stats:::arima(x=x,order=order,xreg=xreg),silent=TRUE)
    }
    else
    {
        if(use.season)
            fit <- try(stats:::arima(x=x,order=order,seasonal=list(order=seasonal,period=m),include.mean=constant),silent=TRUE)
        else
            fit <- try(stats:::arima(x=x,order=order,include.mean=constant),silent=TRUE)
    }
    if(class(fit) != "try-error")
    {
        if(diffs==1 & constant)
        {
            fitnames <- names(fit$coef)
            fitnames[length(fitnames)] <- "drift"
            names(fit$coef) <- fitnames
            fit$call$xreg <- xreg
        }
        if(!is.na(fit$aic))
        {
            npar <- order[1]+order[3]+seasonal[1]+seasonal[3]+constant
            fit$bic <- fit$aic + npar*(log(n) - 2)
            fit$aicc <- fit$aic + 2*npar*(n/(n-npar-1) - 1)
            fit$ic <- switch(ic,bic=fit$bic,aic=fit$aic,aicc=fit$aicc)
        }
        else
            fit$aic <- fit$bic <- fit$aicc <- fit$ic <- 1e20
        # Check for unit roots
        minroot <- 2
        if(order[1] + seasonal[1] > 0)
            minroot <- min(minroot,abs(polyroot(c(1,-fit$model$phi))))
        if(order[3] + seasonal[3] > 0)
            minroot <- min(minroot,abs(polyroot(c(1,fit$model$theta))))
        if(minroot < 1 + 1e-3)
            fit$ic <- 1e20 # Don't like this model
        if(trace)
            cat("\n",arima.string(fit),":",fit$ic)
        return(fit)
    }
    else
    {
        if(trace)
        {
            cat("\n ARIMA(",order[1],",",order[2],",",order[3],")",sep="")
            if(use.season)
                cat("(",seasonal[1],",",seasonal[2],",",seasonal[3],")[",m,"]",sep="")
            if(constant & (order[2]+seasonal[2] == 0))
                cat(" with non-zero mean")
            else if(constant & (order[2]+seasonal[2] == 1))
                cat(" with drift        ")
            else if(!constant & (order[2]+seasonal[2] == 0))
                cat(" with zero mean    ")
            else
                cat("         ")
            cat(" :",1e20,"*")
        }
        return(list(ic=1e20))
    }
}

newmodel <- function(p,d,q,P,D,Q,constant,results)
{
    n <- nrow(results)
    for(i in 1:n)
    {
        if(identical(c(p,d,q,P,D,Q,constant),results[i,1:7]))
            return(FALSE)
    }
    return(TRUE)
}

arima.string <- function(object)
{
    order <- object$arma[c(1,6,2,3,7,4,5)]
    result <- paste("ARIMA(",order[1],",",order[2],",",order[3],")",sep="")
    if(order[7]>1 & sum(order[4:6]) > 0)
        result <- paste(result,"(",order[4],",",order[5],",",order[6],")[",order[7],"]",sep="")
    if(is.element("constant",names(object$coef)))
        result <- paste(result,"with non-zero mean")
    else if(is.element("drift",names(object$coef)))
        result <- paste(result,"with drift        ")
    else if(order[2]==0 & order[5]==0)
        result <- paste(result,"with zero mean    ")
    else
        result <- paste(result,"                  ")
    return(result)
}
