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
make_model <- function(mixture="Poisson",K,T,M=NULL,lambda=NULL,p=NULL){
  
  ## Initialize model list
  model <- list(mixture=mixture,K=K,T=T)
  
  ## Add M if defined
  if(!is.null(M))
    model$M <- M
  
  ## Add density and simulation functions for given mixture
  if(mixture == "Poisson"){
    model$dN <- function(N,k,pars,log) dpois(N,pars$lambda[k],log=log)
    model$rN <- function(k,pars) rpois(1,pars$lambda[k])
  }
  else if(mixture == "Negative Binomial"){
    model$dN <- function(N,k,pars,log) dnbinom(N,mu=pars$lambda[k],size=pars$alpha[k],log=log)
    model$rN <- function(k,pars) rnbinom(N,mu=pars$lambda[k],size=pars$alpha[k])
  }
  else
    stop("Only the Poisson and Negative Binomial mixture models are defined so far.\n")
  
  ## Create link functions for abundance and detection
  ## Default
  model$lambda <- list(formula=~1,link=make.link("log"))
  model$p <- list(formula=~ 1,link=make.link("logit"))
  
  ## Override if alternatives are specified
  if(!is.null(lambda$formula))
    model$lambda$formula <- lambda$formula
  if(!is.null(lambda$link))
    model$lambda$link <- make.link(lambda$link)

  if(!is.null(p$formula))
    model$p$formula <- p$formula
  if(!is.null(p$link))
    model$p$link <- make.link(p$link)
  
  return(model)
}