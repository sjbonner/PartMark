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
lhd_mt <- function(model,data,pars,upper_limit,log=TRUE){
  ## Add contributions for individual sites
  tmp <- sum(sapply(1:model$K,lhd_mt_site,model=model,data=data,pars=pars,upper_limit=upper_limit,log=TRUE))    
  
  if(log)
    return(tmp)
  else
    return(exp(tmp))
}

#' Likelihood for Model M_t (alternative)
#'
#' @param model List of model components
#' @param data List of data components
#' @param parsvec Vector of parameter values
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
lhd_m0_wrap <- function(parsvec,model,data,upper_limit,log=TRUE){
 
  ## Map parameter vector to parameter list
  if(model$mixture=="Poisson"){
    ## First parameter is log(\lambda)
    pars$lambda <- exp(parsvec[1])
    
    ## Further parameters are logit capture probabilitlies in site x visit order
    for(k in 1:model$K){
      pars$p[[k]] <- rep(ilogit(parsvec[2]),model$T[k])
    }
  }
  else
    stop("Only the Poisson mixture model is defined so far.\n")
  
  ## Compute lhd
  lhd_mt(model,data,pars,upper_limit,log = log)
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
lhd_mt_site <- function(k,model,data,pars,upper_limit,log=TRUE){
  ## Compute lower limit of sum
  lower_limit <- data[[k]]$n
  
  ## Sum over complete data likelihood
  tmp <- sum(sapply(lower_limit:upper_limit[k],function(N){
    data[[k]]$N <- N
    cdl_mt_site(k,model,data,pars,log=FALSE)
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
cdl <- function(model,data,pars,log=TRUE){
  ## Compute complete data likelihood for each site
  cdl <- sapply(1:model$K,cdl_mt_site,model=model,data=data,pars=pars,log=TRUE)
  
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
cdl_mt_site <- function(k,model,data,pars,log=TRUE){
    
  # 1. Abundance
  cdl1 <- model$dN(data[[k]]$N,pars,log=TRUE)
  
  # 2. Detections
  cdl2 <- lfactorial(data[[k]]$N) -
    lfactorial(data[[k]]$N-data[[k]]$n) -
    lfactorial(data[[k]]$n) +
    sum(data[[k]]$y * log(pars$p[[k]]) +
          (data[[k]]$N - data[[k]]$y) * log(1-pars$p[[k]]))


  if(log)
    return(cdl1 + cdl2)
  else
    return(exp(cdl1 + cdl2))
}
