#' Likelihood for Unmarked  Model
#'
#' @param model List of model components
#' @param data List of data components
#' @param pars List of parameters
#' @param upper_limit Vector of upper limits on summation over abundance for each site
#' @param log If TRUE compute log-likelihood
#'
#' @return
#' @export
#'
#' @examples
lhd_um <- function(model,data,pars,upper_limit,log=TRUE){
  ## Add contributions for individual sites
  tmp <- sum(sapply(1:model$K,lhd_um_site,model=model,data=data,pars=pars,upper_limit=upper_limit,log=TRUE))    
  
  if(log)
    return(tmp)
  else
    return(exp(tmp))
}

#' Likelihood for Unmarked  Model (alternative)
#'
#' @param beta Vector of parameter values
#' @param model List of model components
#' @param data List of data components
#' @param upper_limit Vector of upper limits on summation over abundance for each site
#'
#' @details This form of the likelihood accepts parameters as a vector instead of in list form. 
#'   Parameters are assumed to be transformed to a scale on which their support is the entire real
#'   line. This function is intended to be called by other functions for computing estimates and 
#'   should not be used directly.
#' 
#' @return
#' @export
#'
#' @examples
lhd_um_wrap <- function(beta,model,data,upper_limit,log=TRUE){
  
  ## Map parameter vector to parameter list
  ## Abundance
  eta <- model$Xlambda %*% beta[1:ncol(model$Xlambda)]
  pars$lambda <- model$lambda$link$linkinv(eta)
    
  if(model$mixture=="Negative Binomial"){
    pars$alpha <- exp(beta[ncol(model$Xlambda)+1])
  }

  ## Detection
  if(model$mixture=="Poisson")
    eta <- model$Xp %*% beta[-(1:ncol(model$Xlambda))]
  else if(model$mixture=="Negative Binomial")
      eta <- model$Xp %*% beta[-(1:(ncol(model$Xlambda)+1))]
    
  l <- 0
  for(k in 1:model$K){
    pars$p[[k]] <- model$p$link$linkinv(eta[l + (1:model$T[k])])
    l <- l + model$T[k]
  }
 
  ## Compute lhd
  lhd_um(model,data,pars,upper_limit,log = log)
}

#' Single Site Component of Likelihood for Unmarked  Model
#'
#' @param model List of model components
#' @param data List of data components
#' @param pars List of parameters
#' @param upper_limit Vector of upper limits on summation over abundance for each site
#' @param log If TRUE compute log-likelihood
#'
#' @return
#' @export
#'
#' @examples
lhd_um_site <- function(k,model,data,pars,upper_limit,log=TRUE){
  ## Compute lower limit of sum
  lower_limit <- max(data[[k]]$y)
  
  ## Sum over complete data likelihood
  tmp <- sum(sapply(lower_limit:upper_limit[k],function(N){
    data[[k]]$N <- N
    cdl_um_site(k,model,data,pars,log=FALSE)
  }))
  
  if(log)
    return(log(tmp))
  else
    return(tmp)
}

#' Complete Data Likelihood for Unmarked  Model
#'
#' @param model List of model components
#' @param data List of data components
#' @param pars List of parameters
#' @param log If TRUE compute log-likelihood
#'
#' @return
#' @export
#'
#' @examples
cdl_um <- function(model,data,pars,log=TRUE){
  ## Compute complete data likelihood for each site
  cdl <- sapply(1:model$K,cdl_um_site,model=model,data=data,pars=pars,log=TRUE)
  
  if(log)
    return(sum(cdl))
  else
    return(exp(sum(cdl)))
}

#' Single Site Component of Complete Data Likelihood for Unmarked 
#'
#' @param k Site number
#' @param model List of model components.
#' @param data List of data components.
#' @param pars List of parameters.
#' @param log If true compute log-likelihood
#'
#' @return
#' @export
cdl_um_site <- function(k,model,data,pars,log=TRUE){
    
  # 1. Abundance
  cdl1 <- model$dN(data[[k]]$N,k,pars,log=TRUE)
  
  # 2. Detections
  cdl2 <- sum(lchoose(data[[k]]$N,data[[k]]$y) + 
          data[[k]]$y * log(pars$p[[k]]) +
          (data[[k]]$N - data[[k]]$y) * log(1-pars$p[[k]]))
  
  if(log)
    return(cdl1 + cdl2)
  else
    return(exp(cdl1 + cdl2))
}