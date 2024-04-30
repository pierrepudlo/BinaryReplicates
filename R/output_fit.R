#' Get credible intervals for the fixed parameters of the Bayesian model,
#' the Bayesian scorings, and the Bayesian prevalence estimate.
#'
#' @name bayesian_computations
#' @title Bayesian computations
#' @param fit The stanfit object return by BayesianFit
#' @param alpha Error level of the credible intervals
#' @return
#' The `credint` function returns the credible interval bounds for the fixed
#' parameters of the Bayesian model. The default error level is 0.05.
#'
#' The `bayesian_scoring` function returns the vector of Bayesian scorings.
#' These scorings are the posterior probabilities
#' that the true latent T[i] is equal to 1.
#'
#' The `bayesian_prevalence_estimate` function returns the posterior mean of the
#' posterior distribution on the prevalence of T = 1.
#'
#' @seealso \code{\link{classification_from_scoring}}
#'
#' @examples
#' data("periodontal")
#' theta <- mean(periodontal$ti)
#' fit <- BayesianFit(periodontal$ni, periodontal$si, chains = 2, iter = 5000)
#' credint(fit)
#' Y_B <- bayesian_scorings(fit)
#' T_B <- classification_from_scoring(Y_B, .4, .6)
#' theta_B <- bayesian_prevalence_estimate(fit)
#' cat("The Bayesian prevalence estimate is ", theta_B, "\n")
#' cat("The prevalence in the data is ", theta, "\n")
#'
#' @rdname bayesian_computations
#' @export
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


#' @rdname bayesian_computations
#' @export
bayesian_scorings <- function(fit){
  # Extract the latent variable
  Ti <- rstan::extract(fit, pars = "Ti")$Ti
  # Compute the posterior probability
  out <- colMeans(Ti)
  return(out)
}

#' @rdname bayesian_computations
#' @export
bayesian_prevalence_estimate <- function(fit){
  # Extract the latent variable
  theta <- rstan::extract(fit, pars = "theta")$theta
  # Compute the posterior mean
  out <- mean(theta)
  return(out)
}
