fitted_pm <- function(output,model,data){
  
  ## Compute fitted values on link scale
  ncoeff_lambda <- ncol(model$lambda$X)
  index <- 1:ncoeff_lambda
  fit_lambda_link <- as.vector(model$lambda$X %*% opt_out$par[index])
  se_lambda_link <- diag(model$lambda$X %*% opt_out$var[index,index] %*% t(model$lambda$X))
  ci_lambda_link <- fit_lambda_link + 1.96 * outer(se_lambda_link,c(-1,1))
  
  if(model$mixture=="Poisson")
    ncoeff_alpha <- 0
  else if(model$mixture=="Negative Binomial")
    ncoeff_alpha <- 1
  
  ncoeff_p <- ncol(model$p$X)
  index <- (ncoeff_lambda+ncoeff_alpha-1) + 1:ncoeff_p
  fit_p_link <- as.vector(model$p$X %*% opt_out$par[index])
  se_p_link <- diag(model$p$X %*% opt_out$var[index,index] %*% t(model$p$X))
  ci_p_link <- fit_p_link + 1.96 * outer(se_p_link,c(-1,1))
  
  ## Backtransform fitted values
  fitted_lambda <- model$lambda$link$linkinv(cbind(fit_lambda_link,ci_lambda_link))
  fitted_p <- model$lambda$link$linkinv(cbind(fit_p_link,ci_p_link))
  
  colnames(fitted_lambda) <- colnames(fitted_p) <- c("Estimate","Lower95","Upper95")
  
  return(list(lambda=fitted_lambda,p=fitted_p))
}

