#' Get credible intervals for the fixed parameters of the Bayesian model
#'
#' @export
#' @param fit The stanfit object return by BayesianFit.
#' @param alpha Error level of the credible intervals
#' @return A data frame with the credible intervals
credint <- function(fit, alpha = .05) {
  # Extract the fixed parameters
  pars <- rstan::extract(fit, pars = c("p", "q", "theta"))
  # Compute the credible intervals
  out <- lapply(pars, function(x) quantile(x, c(alpha/2, 1-alpha/2)))
  out <- t(as.data.frame(out))
  out <- as.data.frame(out)
  colnames(out) <- c("lower", "upper")
  rownames(out) <- c("p", "q", "theta")
  return(out)
}

#' Get the posterior probability that T=1 from the Bayesian model
#'
#' @export
#' @param fit The stanfit object return by BayesianFit.
#' @return A vector of posterior probabilities that T=1
postprobT <- function(fit){
  # Extract the latent variable
  Ti <- rstan::extract(fit, pars = "Ti")$Ti
  # Compute the posterior probability
  out <- colMeans(Ti)
  return(out)
}
