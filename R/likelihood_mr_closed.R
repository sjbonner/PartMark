#' Likelihood for Model M_t
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
lhd_mr <- function(beta=NULL,model,data,pars=NULL,upper_limit,log=TRUE){
  ## Map regression coefficients to natural parameters if necessary
  if(is.null(pars)){
    if(is.null(beta))
      stop("You must either supply beta or pars as arguments.")
    else
      pars <- map_parameters_mr(beta,model,data)
  }
  
  ## Add contributions for individual sites
  tmp <- sum(sapply(1:model$K,lhd_mr_site,model=model,data=data,pars=pars,upper_limit=upper_limit,log=TRUE))    
  
  if(log)
    return(tmp)
  else
    return(exp(tmp))
}

#' Map Coefficients to Natural Parameters (Fully marked)
#'
#' @param beta Vector of parameter values
#' @param model List of model components
#' @param data List of data components
#'
#' @details 
#' @return
#' @export
#'
#' @examples
map_parameters_mr <- function(beta,model,data){
 
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

  return(pars)
}

#' Single Site Component of Likelihood for Model M_t
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
lhd_mr_site <- function(k,model,data,pars,upper_limit,log=TRUE){
  ## Compute lower limit of sum
  lower_limit <- data[[k]]$n
  
  ## Sum over complete data likelihood
  tmp <- sum(sapply(lower_limit:upper_limit[k],function(N){
    data[[k]]$N <- N
    cdl_mr_site(k,model,data,pars,log=FALSE)
  }))
  
  if(log)
    return(log(tmp))
  else
    return(tmp)
}

#' Complete Data Likelihood for Model M_t
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
cdl_mr <- function(model,data,pars,log=TRUE){
  ## Compute complete data likelihood for each site
  cdl <- sapply(1:model$K,cdl_mr_site,model=model,data=data,pars=pars,log=TRUE)
  
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
cdl_mr_site <- function(k,model,data,pars,log=TRUE){
    
  # 1. Abundance
  cdl1 <- model$dN(data[[k]]$N,k,pars,log=TRUE)
  
  # 2. Detections
  cdl2 <- lchoose(data[[k]]$N,data[[k]]$n) +
    sum(data[[k]]$y * log(pars$p[[k]]) +
          (data[[k]]$N - data[[k]]$y) * log(1-pars$p[[k]]))
  if(log)
    return(cdl1 + cdl2)
  else
    return(exp(cdl1 + cdl2))
}
