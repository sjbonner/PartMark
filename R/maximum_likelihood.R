#' Fit N-Mix/MR Model via Maximum Likelihood Inference
#'
#' @param model List of model components
#' @param data List of data components
#' @param upper_limit Vector of upper limits on summation over abundance for each site
#' @param method Optimization method.
#' @param trace Trace level for optim (defaults to none).
#' @param maxit Maximum iterations in optim (defaults to optim default).
#' @param submodel One of "pm" (Partial Marking), "mr" (Mark Recapture), "um" (Unmarked)
#'
#' @return
#' @export
#'
#' @examples
ml_fit <-
  function(model,
           data,
           upper_limit,
           submodel = "pm",
           inits=NULL,
           method = "BFGS",
           trace = NULL,
           maxit = NULL) {

   ## Select appropriate likelihood 
    lhd_wrap <- switch(submodel,
                       pm = lhd_wrap,
                       mr = lhd_wrap,
                       um = lhd_um,
                       stop("Unknown submodel", submodel, ".\n\n"))
    
    
    ## Set initial values
    if(is.null(inits)){
      if (model$mixture == "Poisson") {
        ## Lambda
        if (submodel == "pm") {
          minN <-
            max(sapply(data, function(list)
              max(list$u + list$Mstar))) # Max of min number of individuals known to be alive at each site
        }
        else if (submodel == "mr") {
          minN <- sapply(data, function(X)
            X$n)
        }
        else if (submodel == "um") {
          minN <- sapply(data, function(X)
            max(X$y))
        }
        
        lambda <- mean(minN / (1 - (1 - .4) ^ model$T[1]))
        
        ## p
        p <- mean(sapply(data,function(X) mean(X$y)/lambda))
        
        ## Transform initial values
        inits <- c(model$lambda$link$linkfun(lambda),rep(0,ncol(model$lambda$X)-1),
                   model$p$link$linkfun(p),rep(0,ncol(model$p$X)-1))
      }
      else if (model$mixture == "Negative Binomial") {
        ## Lambda
        if (submodel == "pm") {
          minN <- sapply(data, function(list)
            max(list$u + list$Mstar)) # Max of min number of individuals known to be alive at each site
        }
        else if (submodel == "mr") {
          minN <- sapply(data, function(X)
            X$n)
        }
        else if (submodel == "um") {
          minN <- sapply(data, function(X)
            max(X$y))
        }
        
        Ninit <- minN / (1 - (1 - .4) ^ model$T[1])
        
        lambda <- mean(Ninit)
        
        alpha <- lambda ^ 2 / (max(c(1.1 * lambda, var(Ninit))) - lambda)
        
        ## p
        p <- mean(sapply(data,function(X) mean(X$y)/lambda))
        
        ## Transform initial values
        inits <- c(model$lambda$link$linkfun(lambda),rep(0,ncol(model$lambda$X)-1),
                   log(alpha),
                   model$p$link$linkfun(p),rep(0,ncol(model$p$X)-1))
      }
    }
    else{
      ## Transform initial values
      if(model$mixture=="Poisson"){
        ## Transform initial values
        inits <- c(model$lambda$link$linkfun(inits$lambda),
                   model$p$link$linkfun(inits$p))
      }
      else if(model$mixture=="Negative Binomial"){
        ## Transform initial values
        inits <- c(model$lambda$link$linkfun(inits$lambda),
                   log(inits$alpha),
                   model$p$link$linkfun(inits$p))
      }
    }
    
    
    ## Set control list for optim
    control <- list(fnscale = -1) ## Maximize!
    
    if (!is.null(trace))
      control$trace <- trace
    if (!is.null(maxit))
      control$maxit <- maxit
    
    ## Maximize log-likelihood
    opt_out <-
      optim(
        inits,
        lhd_wrap,
        control = control,
        method = method,
        hessian = TRUE,
        model = model,
        data = data,
        upper_limit = upper_limit
      )
    
    ## Add names to parameters
    if(model$mixture=="Poisson")
      names(opt_out$par) <- 
      colnames(opt_out$hessian) <-
      rownames(opt_out$hessian) <- c(paste0("lambda:",colnames(model$lambda$X)),
                                    paste0("p:",colnames(model$p$X)))

    if(model$mixture=="Negative Binomial")
      names(opt_out$par) <- 
      colnames(opt_out$hessian) <-
      rownames(opt_out$hessian) <- c(paste0("lambda:",colnames(model$lambda$X)),
                                     "alpha",
                                     paste0("p:",colnames(model$p$X)))
    
    ## Compute standard errors and CIs on link scale
    opt_out$var <- -solve(opt_out$hessian)
    ses <- sqrt(diag(opt_out$var))
    ci <- opt_out$par + 1.96 * outer(ses, c(-1, 1))
    
    estimates <- data.frame(
      Estimate = opt_out$par,
      SE = ses,
      Lower95 = ci[, 1],
      Upper95 = ci[, 2]
    )
    
    ## Return output
    return(list(
      optim = opt_out,
      inits = inits,
      estimates = estimates
    ))
  }