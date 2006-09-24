best.arima <- function(x,d=ndiffs(x),D=0,max.p=5,max.q=5,max.P=2,max.Q=2,max.order=5,stationary=FALSE,
    drift=TRUE, trend=FALSE, ic=c("aic","bic"))
{
    ic <- match.arg(ic)
    m <- frequency(x)
    n <- length(x)
    seasonal <- (m > 1)
    best.ic <- 1e9
    oldwarn <- options()$warn
    options(warn=-1)
    on.exit(options(warn=oldwarn))
    if(!seasonal)
        max.P <- max.Q <- 0
    if(stationary)
    {
        d <- 0
        drift <- trend <- FALSE
    }
    for(i in 0:max.p)
    {
        for(j in 0:max.q)
        {
            for(I in 0:max.P)
            {
                for(J in 0:max.Q)
                {
                    for(KK in 0:(2-d))
                    {
                        if(i+j+I+J <= max.order)
                        {
                            if(KK==0) # Zero mean
                            {
                                if(seasonal)
                                    fit <- try(arima(x,order=c(i,d,j),seasonal=list(order=c(I,D,J),period=m),include.mean=FALSE),silent=TRUE)
                                else
                                    fit <- try(arima(x,order=c(i,d,j),include.mean=FALSE),silent=TRUE)
                            }
                            else if(d==0 & KK==1) # Allow non-zero mean when d=0
                            {
                                if(seasonal)
                                    fit <- try(arima(x,order=c(i,d,j),seasonal=list(order=c(I,D,J),period=m),include.mean=TRUE),silent=TRUE)
                                else
                                    fit <- try(arima(x,order=c(i,d,j),include.mean=TRUE),silent=TRUE)
                            }
                            else if(trend & d==0 & KK==2 & !stationary) # Allow linear trend when d=0
                            {
                                xdrift <- seq(1,n)
                                if(seasonal)
                                    fit <- try(arima(x,order=c(i,d,j),seasonal=list(order=c(I,D,J),period=m),xreg=xdrift,include.mean=TRUE),silent=TRUE)
                                else
                                    fit <- try(arima(x,order=c(i,d,j),xreg=xdrift,include.mean=TRUE),silent=TRUE)
                                if(class(fit) != "try-error")
                                {
                                    fitnames <- names(fit$coef)
                                    fitnames[length(fitnames)] <- "drift"
                                    names(fit$coef) <- fitnames
                                }
                            }
                            else if(drift & d==1 & KK==1)# allow drift when d=1
                            {
                                xdrift <- seq(1,n)
                                if(seasonal)
                                    fit <- try(arima(x,order=c(i,d,j),seasonal=list(order=c(I,D,J),period=m),xreg=xdrift),silent=TRUE)
                                else
                                    fit <- try(arima(x,order=c(i,d,j),xreg=xdrift),silent=TRUE)
                                if(class(fit) != "try-error")
                                {
                                    fitnames <- names(fit$coef)
                                    fitnames[length(fitnames)] <- "drift"
                                    names(fit$coef) <- fitnames
                                }
                            }
                            else # Do nothing
                                class(fit) <- "try-error"
                            if(class(fit) != "try-error")
                            {
                                if(!is.na(fit$aic))
                                {
                                    npar <- i+j+I+J+KK
                                    fit$bic <- fit$aic + npar*(log(n) - 2)
                                    fit$ic <- if(ic=="bic")
                                                fit$bic
                                              else if(ic=="aic")
                                                fit$aic
                                              else
                                                stop("This shouldn't happen")
                                    if(fit$code==0 & fit$ic < best.ic)
                                    {
                                        order <- c(i,d,j)
                                        Order <- c(I,D,J)
                                        best.ic <- fit$ic
                                        bestfit <- fit
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if(exists("bestfit"))
        bestfit$x <- x
    else
    {
        options(show.error.messages=olderror,warn=oldwarn)
        stop("No ARIMA model able to be estimated")
    }
    bestfit$series <- deparse(substitute(x))
    bestfit$ic=NULL
    return(bestfit)
}

ndiffs <- function(x,alpha=0.05)
{
    require(tseries)
    x <- c(na.omit(c(x)))
    d <- 0
    oldwarn <- options(warn=-1)
    p.v <- kpss.test(x)$p.value
    if(is.na(p.v))
    {
        options(warn=oldwarn$warn)
        return(d)
    }
    while(p.v < alpha & d<2)
    {
        x <- diff(x)
        d <- d+1
        p.v <- kpss.test(x)$p.value
        if(is.na(p.v))
            return(d-1)
    }
    options(warn=oldwarn$warn)
    return(d)
}


forecast.Arima <- function (object, h = ifelse(object$arma[5] > 1, 2 * object$arma[5], 10),
    conf = c(80, 95), fan=FALSE, xreg = NULL,...)
{
#    use.constant <- is.element("constant",names(object$coef))
    use.drift <- is.element("drift", names(object$coef))
    if (is.element("x", names(object)))
        x <- object$x
    else
        x <- eval(parse(text = object$series))
    usexreg <- (!is.null(xreg) | use.drift)# | use.constant)
#    if(use.constant)
#        xreg <- as.matrix(rep(1,h))
    if (use.drift)
    {
        n <- length(x)
        if(!is.null(xreg))
            xreg <- cbind(xreg,(1:h)+n)
        else
            xreg <- as.matrix((1:h)+n)
    }
    if(usexreg)
        pred <- predict(object, n.ahead = h, newxreg = xreg)
    else
        pred <- predict(object, n.ahead = h)

    if(fan)
        conf <- seq(50,99,by=1)
    else
    {
        if(min(conf) > 0 & max(conf) < 1)
            conf <- 100*conf
        else if(min(conf) < 0 | max(conf) > 99.99)
            stop("Confidence limit out of range")
    }

    nint <- length(conf)
    lower <- matrix(NA, ncol = nint, nrow = length(pred$pred))
    upper <- lower
    for (i in 1:nint)
    {
        qq <- qnorm(0.5 * (1 + conf[i]/100))
        lower[, i] <- pred$pred - qq * pred$se
        upper[, i] <- pred$pred + qq * pred$se
    }
    colnames(lower) = colnames(upper) = paste(conf, "%", sep = "")
    method <- paste("ARIMA(", object$arma[1], ",", object$arma[6],
        ",", object$arma[2], ")", sep = "")
    if (object$arma[5] > 1)
        method <- paste(method, "(", object$arma[3], ",", object$arma[7],
            ",", object$arma[4], ")", object$arma[5], sep = "")
    return(structure(list(method = method, model = object, conf = conf,
        mean = pred$pred, lower = lower, upper = upper, x = x,
        xname = deparse(substitute(x)), fitted = fitted(object), residuals = residuals(object)),
        class = "forecast"))
}


forecast.ar <- function(object,h=10,conf=c(80,95),fan=FALSE,...)
{
    pred <- predict(object,n.ahead=h)
    if(fan)
        conf <- seq(50,99,by=1)
    else
    {
        if(min(conf) > 0 & max(conf) < 1)
            conf <- 100*conf
        else if(min(conf) < 0 | max(conf) > 99.99)
            stop("Confidence limit out of range")
    }
    nint <- length(conf)
    lower <- matrix(NA,ncol=nint,nrow=length(pred$pred))
    upper <- lower
    for(i in 1:nint)
    {
        qq <- qnorm(0.5*(1+conf[i]/100))
        lower[,i] <- pred$pred - qq*pred$se
        upper[,i] <- pred$pred + qq*pred$se
    }
    colnames(lower) = colnames(upper) = paste(conf,"%",sep="")
    method <- paste("AR(",object$order,")",sep="")
    x <- as.ts(eval(parse(text=object$series)))
    f=frequency(x)
    res <- ts(object$resid[-(1:object$order)],start=tsp(x)[1]+object$order/f,f=f)
    fits <- x-res

    return(structure(list(method=method,model=object,conf=conf,mean=pred$pred,lower=lower,upper=upper,
        x=x, xname = deparse(substitute(x)), fitted=fits,residuals=res)
        ,class="forecast"))
}

# Extract errors from ARIMA model (as distinct from residuals)
arima.errors <- function(z)
{
    if(!is.list(z))
        stop("z must be a list")
    if(is.element("x",names(z)))
        x <- z$x
    else
    {
        series.name <- z$series
        if(is.null(series.name))
            stop("missing component series in argument z\n")
        x <- eval(parse(text = series.name))
    }
    if(!is.element("xreg",names(z$call)))
        return(x)
    else
        xreg <- z$call$xreg
    norder <- sum(z$arma[1:4])
    if(is.element("intercept",names(z$coef)))
        xreg <- cbind(rep(1,length(x)),xreg)
    return(ts(x - xreg %*% as.matrix(z$coef[(norder+1):length(z$coef)]),f=frequency(x),s=start(x)))
}

# Return one-step fits
fitted.Arima <- function(object,...)
{
    if(is.element("x",names(object)))
        x <- object$x
    else
        x <- eval(parse(text=object$series))

    return(x - object$residuals)
}

# Calls arima from stats package and adds data to the returned object
arima <- function(x, order = c(0, 0, 0),
      seasonal = list(order = c(0, 0, 0), period = NA),
      xreg = NULL, include.mean = TRUE, transform.pars = TRUE,
      fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"),
      n.cond, optim.control = list(), kappa = 1e6, model=NULL)
{
    if(!is.null(model))
        return(arima2(x,model))
    if(is.null(xreg))
        tmp <- stats:::arima(x=x,order=order,seasonal=seasonal,include.mean=include.mean,
            transform.pars=transform.pars,fixed=fixed,init=init,method=method,n.cond=n.cond,optim.control=optim.control,kappa=kappa)
    else
        tmp <- stats:::arima(x=x,order=order,seasonal=seasonal,xreg=xreg,include.mean=include.mean,
               transform.pars=transform.pars,fixed=fixed,init=init,method=method,n.cond=n.cond,optim.control=optim.control,kappa=kappa)
    tmp$x <- x
    tmp$series <- deparse(substitute(x))
    if(!is.null(xreg))
        tmp$call$xreg <- xreg
    return(tmp)
}

# Refits the model to new data x
arima2 <- function (x, model)
{
    use.drift <- is.element("drift",names(model$coef))
    use.intercept <- is.element("intercept",names(model$coef))
    if(use.drift)
        xreg <- length(model$x)+(1:length(x))

    if(model$arma[5]>1 & sum(abs(model$arma[c(3,4,7)]))>0) # Seasonal model
    {
        if(use.drift)
            refit <- arima(x,order=model$arma[c(1,6,2)],seasonal=list(order=model$arma[c(3,7,4)],period=model$arma[5]),
                fixed=model$coef,include.mean=use.intercept,xreg=xreg)
        else
            refit <- arima(x,order=model$arma[c(1,6,2)],seasonal=list(order=model$arma[c(3,7,4)],period=model$arma[5]),
                fixed=model$coef,include.mean=use.intercept)
    }
    else if(length(model$coef)>0) # Nonseasonal model with some parameters
    {
        if(use.drift)
           refit <- arima(x,order=model$arma[c(1,6,2)],fixed=model$coef,xreg=xreg,include.mean=use.intercept)
        else
            refit <- arima(x,order=model$arma[c(1,6,2)],fixed=model$coef,include.mean=use.intercept)
    }
    else # No parameters
        refit <- arima(x,order=model$arma[c(1,6,2)],include.mean=FALSE)
    refit$var.coef <- matrix(0,length(refit$coef),length(refit$coef))
    return(refit)
}

print.Arima <- function (x, digits = max(3, getOption("digits") - 3), se = TRUE,
    ...)
{
    cat("Series:",x$series,"\n")
    order <- x$arma[c(1,6,2,3,7,4,5)]
    cat("ARIMA(",order[1],",",order[2],",",order[3],")",sep="")
    if(order[7]>1)
        cat("(",order[4],",",order[5],",",order[6],")[",order[7],"]",sep="")
    cat(" model\n")
    if(!is.null(x$call$xreg))
    {
        cat("\nRegression variables fitted:\n")
        x$call$xreg <- as.matrix(x$call$xreg)
        for(i in 1:3)
            cat("  ",x$call$xreg[i,],"\n")
        cat("   . . .\n")
        for(i in 1:3)
            cat("  ",x$call$xreg[nrow(x$call$xreg)-3+i,],"\n")
    }
    if (length(x$coef) > 0) {
        cat("\nCoefficients:\n")
        coef <- round(x$coef, digits = digits)
        if (se && nrow(x$var.coef)) {
            ses <- rep(0, length(coef))
            ses[x$mask] <- round(sqrt(diag(x$var.coef)), digits = digits)
            coef <- matrix(coef, 1, dimnames = list(NULL, names(coef)))
            coef <- rbind(coef, s.e. = ses)
        }
        print.default(coef, print.gap = 2)
    }
    cm <- x$call$method
    if (is.null(cm) || cm != "CSS")
        cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
            ":  log likelihood = ", format(round(x$loglik, 2)),
            ",  aic = ", format(round(x$aic, 2)), "\n", sep = "")
    else cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
        ":  part log likelihood = ", format(round(x$loglik, 2)),
        "\n", sep = "")
    invisible(x)
}
