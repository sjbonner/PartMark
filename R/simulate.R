#' Simulate data from the N-Mix/MR model
#'
#' @param model List of model components.
#' @param data List of data components.
#' @param pars List of parameters.
#'
#' @return
#' @export
#'
simulate_Nmix_MR <- function(model,pars){
  ## Simulate data on a site by site basis
  data <- lapply(1:model$K,function(k) 
    simulate_Nmix_MR_site(k,model,pars))
  
  ## Reconfigure data
  list(N=sapply(data,function(X) X$N),
       marked=lapply(data,function(X) X$marked),
       unmarked=lapply(data,function(X) X$unmarked),
       partmarked=lapply(data,function(X) X$partmarked))
}

simulate_Nmix_MR_site <- function(k,model,pars){
  
  ## Simulate abundance
  N <- model$rN(k,pars)

  ## Generate capture histories
  W <- t(replicate(N,rbinom(model$T[k],1,pars$p[[k]])))
  
  ## Mark individuals
  marked <- data.frame(ID=c(),Occasion=c())
  unmarked <- 1:N
  
  for(t in which(model$M[1:model$T[k]]>0)){
    # Select individuals to mark
    caps <- which(W[unmarked,t]==1)
    
    if(length(caps)<=model$M[t])
      m <- caps
    else
      m <- sort(sample(caps,size=model$M[t]))
    
    # Add to data frame
    marked <- rbind(marked,data.frame(ID=unmarked[m],Occasion=t))
    
    # Remove from list of unmarked individuals
    unmarked <- unmarked[-m]
  }

  # 1) Extract capture histories of marked individuals after time of marking
  if(nrow(marked)>0){
    W_marked <- t(sapply(1:nrow(marked),function(i){
      c(rep(0,marked[i,"Occasion"]-1),
        W[marked[i,"ID"],
          marked[i,"Occasion"]:model$T[k]])
    }))
    
    ## Count time of marking for marked individuals
    M <- sapply(1:model$T[k],function(t) sum(marked$Occasion==t))

    ## Compute counts of previously marked individuals
    v <- apply(W_marked,2,sum) - M
  }
  else{
    M <- v <- rep(0,model$T[[k]])
  }
  
  ## Compute counts of unmarked individuals
  u <- apply(W,2,sum) - v

  ## Summarize data for unmarked analysis
  unmarked.data <- list(y=apply(W,2,sum))
  
  ## Summarize data for marked analysis
  # 1) Subsample captures so that number of captured individuals is approximatley \sum(M)
  caps1 <- which(W==1)
  ncap <- sum(W)
  ncap2 <- round(N*model$T[k]*(1-(1-sum(M)/N)^(1/model$T[k])))
  remove <- sample(caps1,ncap - ncap2)
  W2 <- W
  W2[remove] <- 0
  
  marked.data <- list(n=sum(apply(W2,1,sum)>0),y=apply(W2,2,sum))
  
  ## Sumarized data for part marked analysis
  partmarked.data <- list(u=u,v=v,y=u+v,M=M,Mstar=c(0,cumsum(M[-model$T[k]])))
  
  ## Return data
  return(list(N=N,
              marked=marked.data,
              unmarked=unmarked.data,
              partmarked=partmarked.data))
}
