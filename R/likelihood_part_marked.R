#' Likelihood for N-Mix/MR Model
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
lhd_pm <- function(beta=NULL,model,data,pars=NULL,upper_limit,log=TRUE){
  ## Map regression coefficients to natural parameters if necessary
  if(is.null(pars)){
    if(is.null(beta))
      stop("You must either supply beta or pars as arguments.")
    else
      pars <- map_parameters_pm(beta,model,data)
  }
  
  ## Add contributions for individual sites
  tmp <- sum(sapply(1:model$K,lhd_pm_site_fast,model=model,data=data,pars=pars,upper_limit=upper_limit,log=TRUE))    
  
  if(log)
    return(tmp)
  else
    return(exp(tmp))
}

#' Map Coefficients to Natural Parameters (Part marked)
#'
#' @param beta Vector of parameter values
#' @param model List of model components
#' @param data List of data components
#'
#' @details 
#' 
#' @return
#' @export
#'
#' @examples
map_parameters_pm <- function(beta,model,data){
  ## Abundance
  eta <- model$lambda$X %*% beta[1:ncol(model$lambda$X)]
  pars <- list(lambda=model$lambda$link$linkinv(eta))
    
  if(model$mixture=="Negative Binomial"){
    pars$alpha <- exp(beta[ncol(model$lambda$X)+1])
  }

  ## Detection
  if(model$mixture=="Poisson")
    eta <- model$p$X %*% beta[-(1:ncol(model$lambda$X))]
  else if(model$mixture=="Negative Binomial")
      eta <- model$p$X %*% beta[-(1:(ncol(model$lambda$X)+1))]
    
  l <- 0
  for(k in 1:model$K){
    pars$p[[k]] <- model$p$link$linkinv(eta[l + (1:model$T[k])])
    l <- l + model$T[k]
  }

  return(pars)
}

#' Single Site Component of Likelihood for N-Mix/MR Model
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
lhd_pm_site <- function(k,model,data,pars,upper_limit,log=TRUE){
  ## Compute lower limit of sum
  lower_limit <- max(data[[k]]$u + data[[k]]$Mstar)
  
  ## Sum over complete data likelihood
  tmp <- sum(sapply(lower_limit:upper_limit[k],function(N){
    data[[k]]$N <- N
    cdl_pm_site(k,model,data,pars,log=FALSE)
  }))
  
  if(log)
    return(log(tmp))
  else
    return(tmp)
}

#' Complete Data Likelihood for N-Mix/MR Model
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
cdl_pm <- function(model,data,pars,log=TRUE){
  ## Compute complete data likelihood for each site
  cdl <- sapply(1:model$K,cdl_pm_site,model=model,data=data,pars=pars,log=TRUE)
  
  if(log)
    return(sum(cdl))
  else
    return(exp(sum(cdl)))
}

#' Single Site Component of Complete Data Likelihood for N-Mix/MR
#'
#' @param k Site number
#' @param model List of model components.
#' @param data List of data components.
#' @param pars List of parameters.
#' @param log If true compute log-likelihood
#'
#' @return
#' @export
cdl_pm_site <- function(k,model,data,pars,log=TRUE){
    
  # 1. Abundance
  cdl1 <- model$dN(data[[k]]$N,k,pars,log=TRUE)
  
  # 2. Detections
  cdl2 <- sum(lchoose(data[[k]]$N-data[[k]]$Mstar,data[[k]]$u) + 
                data[[k]]$y * log(pars$p[[k]]) +
                (data[[k]]$N - data[[k]]$y) * log(1-pars$p[[k]]))
  if(log)
    return(cdl1 + cdl2)
  else
    return(exp(cdl1 + cdl2))
}

#' Single Site Component of Likelihood for Unmarked Model 
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
lhd_pm_site_fast <- function(k,model,data,pars,upper_limit,log = TRUE,debug = FALSE){
  
  if(debug)
    browser()
  
  ## Compute lower limit of sum
  lower_limit <- max(data[[k]]$u + data[[k]]$Mstar)
  
  ## Compute individual likelihood components
  # Abundance
  cdl1 <- model$dN(lower_limit:upper_limit[k],k,pars,log = TRUE)
  
  # Combinatorial term
  cdl2 <- sum(lchoose(lower_limit-data[[k]]$Mstar,data[[k]]$u)) +
    #cumsum(c(0,model$T[k]*log((lower_limit + 1):upper_limit[k]-data[[k]]$Mstar))) - # Need an outer here too
    cumsum(c(0,apply(log(outer((lower_limit + 1):upper_limit[k],data[[k]]$Mstar,"-")),1,sum))) - 
    cumsum(c(0,apply(log(outer((lower_limit + 1):upper_limit[k],data[[k]]$u+data[[k]]$Mstar,"-")),1,sum)))
  
  # Detections
  cdl3 <- sum(data[[k]]$y * log(pars$p[[k]])) +
    outer((lower_limit:upper_limit[k]),data[[k]]$y,"-") %*% log(1 - pars$p[[k]])
  
  tmp <- sum(exp(cdl1 + cdl2 + cdl3))
  
  if(log)
    return(log(tmp))
  else
    return(tmp)
}

