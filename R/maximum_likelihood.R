#' Fit N-Mix/MR Model via Maximum Likelihood Inference
#'
#' @param model List of model components
#' @param data List of data components
#' @param upper_limit Vector of upper limits on summation over abundance for each site
#' @param method Optimization method.
#' @param trace Trace level for optim (defaults to none).
#' @param maxit Maximum iterations in optim (defaults to optim default).
#' @param Xlambda Matrix of covariates of abundance model
#' @param Xp Matri of covariats of detection model
#' @param submodel One of "pm" (Partial Marking), "mr" (Mark Recapture), "um" (Unmarked)
#'
#' @return
#' @export
#'
#' @examples
ml_fit <- function(model, data, Xlambda=NULL, Xp=NULL, upper_limit,submodel="pm",method="BFGS",trace=NULL,maxit=NULL) {

  ## Identify likelihood wrapper
  lhd_wrap <- switch(submodel,pm=lhd_pm_wrap,
                     mt=lhd_mt_wrap,
                     mr=lhd_mr_wrap,
                     um=lhd_um_wrap,
                     stop("Unknown submodel",submodel,".\n\n"))
  
  ## Build design matrices
  ## Abundance
  if(is.null(Xlambda))
    Xlambda <- data.frame(Intercept=rep(1,model$K))

  model$Xlambda <- model.matrix(model$p$formula,Xlambda)

  ## Detection
  if(is.null(Xp))
    Xp <- data.frame(Intercept=rep(1,sum(model$T)))
  
  model$Xp <- model.matrix(model$p$formula,Xp)
    
  ## Set initial values
  if (model$mixture == "Poisson") {
    ## Lambda
    if(submodel=="pm"){
    minN <-
      max(sapply(data, function(list)
        max(list$u + list$Mstar))) # Max of min number of individuals known to be alive at each site
    }
    else if(submodel=="mr")
      minN <- sapply(data,function(X) X$n)
    
    lambda <- mean(minN / (1 - (1 - .4) ^ model$T[1]))
    
    ## p
    p <- .4
    
    ## Transform
    inits <- c(model$lambda$link$linkfun(lambda), 
               model$p$link$linkfun(p))
  }
  else if (model$mixture == "Negative Binomial") {
    ## Lambda
    if(submodel=="pm"){
      minN <- sapply(data, function(list)
          max(list$u + list$Mstar)) # Max of min number of individuals known to be alive at each site
    }
    else if(submodel=="mr")
      minN <- sapply(data,function(X) X$n)
    
    Ninit <- minN / (1 - (1 - .4) ^ model$T[1])
    
    lambda <- mean(Ninit)
    
    alpha <- lambda^2/(max(c(1.1*lambda,var(Ninit))) - lambda)
      
    ## p
    p <- .4
    
    ## Transform
    inits <- c(model$lambda$link$linkfun(lambda), 
               log(alpha),
               model$p$link$linkfun(p))
  }
  else
    stop("Only the Poisson and Negative Binomial mixtures model are defined so far.\n")
  
  ## Set control list for optim
  control <- list(fnscale = -1) ## Maximize!
  
  if(!is.null(trace))
    control$trace <- trace
  if(!is.null(maxit))
    control$maxit <- maxit
  
  ## Maximize log-likelihood
  opt_out <-
    optim(
      inits,
      lhd_wrap,
      control = control,
      method=method,
      hessian = TRUE,
      model = model,
      data = data,
      upper_limit = upper_limit
    )
  
  ## Compute standard errors and CIs on link scale
  ses <- sqrt(-1*diag(solve(opt_out$hessian)))
  ci <- opt_out$par + 1.96*outer(ses,c(-1,1))
  
  estimates_link <- data.frame(Estimate=opt_out$par,
                               SE=ses,
                               Lower95=ci[,1],
                               Upper95=ci[,2])
  
  ## Backtransform
  if (model$mixture == "Poisson") {
    estimates <- rbind(model$lambda$link$linkinv(estimates_link[1,c(1,3,4)]),
                       model$p$link$linkinv(unlist(estimates_link[2,c(1,3,4)])))
  }
  else if (model$mixture == "Negative Binomial") {
    estimates <- rbind(model$lambda$link$linkinv(estimates_link[1,c(1,3,4)]),
                       exp(estimates_link[2,c(1,3,4)]),
                       model$p$link$linkinv(unlist(estimates_link[3,c(1,3,4)])))
  }
  
  ## Return output
  return(list(optim = opt_out, 
              inits=inits, 
              link=estimates_link, 
              natural = estimates))
}