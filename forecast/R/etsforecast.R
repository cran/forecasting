predict.ets <- function(...){forecast.ets(...)}

forecast.ets <- function(obj, h=ifelse(obj$m>1, 2*obj$m, 10),
    conf=c(80,95), fan=FALSE, simulate=FALSE, bootstrap=FALSE, npaths=5000,...)
{
    # Check inputs
    if(h>2000 | h<=0)
        stop("Forecast horizon out of bounds")
    if(fan)
        conf <- seq(50,99,by=1)
    else
    {
        if(min(conf) > 0 & max(conf) < 1)
            conf <- 100*conf
        else if(min(conf) < 0 | max(conf) > 99.99)
            stop("Confidence limit out of range")
    }

    n <- length(obj$x)
    damped <- as.logical(obj$components[4])

    if(simulate)
        f <- pegelsfcast.C(h,obj,conf=conf,bootstrap=bootstrap,npaths=npaths)
    else if(obj$components[1]=="A" & is.element(obj$components[2],c("A","N")) & is.element(obj$components[3],c("N","A")))
        f <- class1(h,obj$states[n+1,],obj$components[2],obj$components[3],damped,obj$m,obj$sigma2,obj$par)
    else if(obj$components[1]=="M" & is.element(obj$components[2],c("A","N")) & is.element(obj$components[3],c("N","A")))
        f <- class2(h,obj$states[n+1,],obj$components[2],obj$components[3],damped,obj$m,obj$sigma2,obj$par)
    else if(obj$components[1]=="M" & obj$components[3]=="M" & obj$components[2]!="M")
        f <- class3(h,obj$states[n+1,],obj$components[2],obj$components[3],damped,obj$m,obj$sigma2,obj$par)
    else
        f <- pegelsfcast.C(h,obj,conf=conf,bootstrap=bootstrap,npaths=npaths)

    tsp.x <- tsp(obj$x)
    if(!is.null(tsp.x))
        start.f <- tsp(obj$x)[2] + 1/obj$m
    else
        start.f <- length(obj$x)+1
    out <- list(model=obj,mean=ts(f$mu,f=obj$m,s=start.f),conf=conf,x=obj$x)
    if(!is.null(f$var))
    {
        out$lower <- out$upper <- matrix(NA,ncol=length(conf),nrow=h)
        for(i in 1:length(conf))
        {
            marg.error <- sqrt(f$var) * abs(qnorm((100-conf[i])/200))
            out$lower[,i] <- out$mean - marg.error
            out$upper[,i] <- out$mean + marg.error
        }
    }
    else if(!is.null(f$lower))
    {
        out$lower <- ts(f$lower,f=obj$m,s=start.f)
        out$upper <- ts(f$upper,f=obj$m,s=start.f)
    }
    else
        warning("No prediction intervals for this model")
    out$method <- obj$method
    out$residuals <- residuals(obj)
    out$fitted <- fitted(obj)
    return(structure(out,class="forecast"))
}

pegelsfcast.C <- function(h,obj,npaths,conf,bootstrap)
{
    y.paths <- matrix(NA,nrow=npaths,ncol=h)
    for(i in 1:npaths)
        y.paths[i,] <- simulate.ets(obj, h, initstate=obj$state[length(obj$x)+1,], bootstrap=bootstrap)
    y.f <- .C("etsforecast",
            as.double(obj$state[length(obj$x)+1,]),
            as.integer(obj$m),
            as.integer(switch(obj$components[2],"N"=0,"A"=1,"M"=2)),
            as.integer(switch(obj$components[3],"N"=0,"A"=1,"M"=2)),
            as.double(ifelse(obj$components[4]=="FALSE",1,obj$par["phi"])),
            as.integer(h),
            as.double(numeric(h)),
        PACKAGE="forecast")[[7]]
    if(abs(y.f[1]+99999) < 1e-7)
        stop("Problem with multiplicative damped trend")

    lower <- apply(y.paths,2,quantile,0.5 - conf/200,type=8)
    upper <- apply(y.paths,2,quantile,0.5 + conf/200,type=8)
    if(length(conf)>1)
    {
        lower <- t(lower)
        upper <- t(upper)
    }
    return(list(mu=y.f,lower=lower,upper=upper))
}

class1 <- function(h,last.state,trendtype,seasontype,damped,m,sigma2,par)
{
    p <- length(last.state)
    H <- matrix(c(1,rep(0,p-1)),nrow=1)
    if(seasontype=="A")
        H[1,p] <- 1
    if(trendtype=="A")
        H[1,2] <- 1
    F <- matrix(0,p,p)
    F[1,1] <- 1
    if(trendtype=="A")
    {
        F[1,2] <- 1
        if(damped)
            F[2,2] <- par["phi"]
        else
            F[2,2] <- 1
    }
    if(seasontype=="A")
    {
        F[p-m+1,p] <- 1
        F[(p-m+2):p,(p-m+1):(p-1)] <- diag(m-1)
    }
    G <- matrix(0,nrow=p,ncol=1)
    G[1,1] <- par["alpha"]
    if(trendtype=="A")
        G[2,1] <- par["beta"]
    if(seasontype=="A")
        G[3,1] <- par["gamma"]
    mu <- numeric(h)
    Fj <- diag(p)
    if(h>1)
    {
        cj <- numeric(h-1)
        for(i in 1:(h-1))
        {
            mu[i] <- H %*% Fj %*% last.state
            Fj <- Fj %*% F
            cj[i] <- H %*% Fj %*% G
        }
        cj2 <- cumsum(cj^2)
        var <- sigma2 * c(1,1+cj2)
    }
    else
        var <- sigma2
    mu[h] <- H %*% Fj %*% last.state

    return(list(mu=mu,var=var,cj=cj))
}

class2 <- function(h,last.state,trendtype,seasontype,damped,m,sigma2,par)
{
    tmp <- class1(h,last.state,trendtype,seasontype,damped,m,sigma2,par)
    theta <- numeric(h)
    theta[1] <- tmp$mu[1]^2
    if(h>1)
    {
        for(j in 2:h)
            theta[j] <- tmp$mu[j]^2 + sigma2 * sum(tmp$cj[1:(j-1)]^2*theta[(j-1):1])
    }
    var <- (1+sigma2)*theta - tmp$mu^2
    return(list(mu=tmp$mu,var=var))
}

class3 <- function(h,last.state,trendtype,seasontype,damped,m,sigma2,par)
{
    p <- length(last.state)
    H1 <- matrix(rep(1,1+(trendtype!="N")),nrow=1)
    H2 <- matrix(c(rep(0,m-1),1),nrow=1)
    if(trendtype=="N")
    {
        F1 <- 1
        G1 <- par["alpha"]
    }
    else
    {
        F1 <- rbind(c(1,1),c(0,ifelse(damped,par["phi"],1)))
        G1 <- rbind(c(par["alpha"],par["alpha"]),c(par["beta"],par["beta"]))
    }
    F2 <- rbind(c(rep(0,m-1),1),cbind(diag(m-1),rep(0,m-1)))

    G2 <- matrix(0,m,m)
    G2[1,m] <- par["gamma"]
    Mh <- matrix(last.state[1:(p-m)]) %*% matrix(last.state[(p-m+1):p],nrow=1)
    Vh <- matrix(0,length(Mh),length(Mh))
    H21 <- H2 %x% H1
    F21 <- F2 %x% F1
    G21 <- G2 %x% G1
    K <- (G2 %x% F1) + (F2 %x% G1)
    mu <- var <- numeric(h)
    for(i in 1:h)
    {
        mu[i] <- H1 %*% Mh %*% t(H2)
        var[i] <- (1+sigma2) * H21 %*% Vh %*% t(H21) + sigma2*mu[i]^2
        vecMh <- c(Mh)
        Vh <- F21 %*% Vh %*% t(F21) + sigma2 * (F21 %*% Vh %*% t(G21) + G21 %*% Vh %*% t(F21) +
            K %*% (Vh + vecMh %*% t(vecMh)) %*% t(K) + sigma2 * G21 %*% (3*Vh + 2*vecMh%*%t(vecMh))%*%t(G21))
        Mh <- F1 %*% Mh %*% t(F2) + G1 %*% Mh %*% t(G2) * sigma2
    }
    return(list(mu=mu,var=var))
}

ses <- function(x,h=10,conf=c(80,95),fan=FALSE,...)
{
    fcast <- forecast(ets(x,"ANN"),h,conf=conf,fan=fan,...)
    fcast$method <- "Simple exponential smoothing"
    fcast$model$call <- match.call()
    return(fcast)
}

holt <- function(x,h=10, damped=FALSE, conf=c(80,95), fan=FALSE, ...)
{
    junk <- forecast(ets(x,"AAN",damped=damped),h,conf=conf,fan=fan,...)
    if(damped)
        junk$method <- "Damped Holt's method"
    else
        junk$method <- "Holt's method"
    junk$model$call <- match.call()
    return(junk)
}

hw <- function(x,h=2*frequency(x),seasonal="additive",damped=FALSE,conf=c(80,95), fan=FALSE, ...)
{
    if(seasonal=="additive")
    {
        junk <- forecast(ets(x,"AAA",damped=damped),h,conf=conf,fan=fan,...)
        junk$method <- "Holt-Winters' additive method"
    }
    else
    {
        junk <- forecast(ets(x,"MAM",damped=damped),h,conf=conf,fan=fan,...)
        junk$method <- "Holt-Winters' multiplicative method"
    }
    junk$model$call <- match.call()
    return(junk)
}

