#' Define N-Mix/MR model
#'
#' @param mixture Name of distribution for abundance
#' @param K Number of sites (scalar)
#' @param T Number of visits to each site (vector of length K)
#' @param M Number of marks applied per visit at each site (vector of length T assumed constant across sites)
#'
#' @details The argument M allows specification of the intended number of animals marked on each visit
#'   to each site. This is only important for simulations in which data is generated. The number of 
#'   individuals marked per visit is recycled across sites and is truncated if less than M[t] individuals
#'   are captured on a given visit at a given site. Simulating data will create the vectors M and Mstar
#'   for each site recording the actual number of marked on each visit and the number marked prior to each
#'   visit. If you are working with real data then these values need to be defined directly in the data
#'   list and not through the model list.
#' @return
#' @export
#'
#' @examples
make_model <- function(mixture="Poisson",K,T,M=NULL){
  
  ## Initialize model list
  model <- list(mixture=mixture,K=K,T=T)
  
  ## Add M if defined
  if(!is.null(M))
    model$M <- M
  
  ## Add density and simulation functions for given mixture
  if(mixture == "Poisson"){
    model$dN <- function(N,pars,log) dpois(N,pars$lambda,log=log)
    model$rN <- function(pars) rpois(1,pars$lambda)
  }
  else
    stop("Only the Poisson mixture model is defined so far.\n")
  
  return(model)
}