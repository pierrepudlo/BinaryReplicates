#' Get credible intervals, Bayesian scores, and prevalence estimate
#' from the `stanfit` object returned by [BayesianFit].
#'
#' @name bayesian_computations
#' @title Bayesian computations
#' @param fit The `stanfit` object return by [BayesianFit]
#' @param level Posterior probability of the credible intervals
#' @details See [BayesianFit] for details on the Bayesian model.
#' @return
#' The `credint` function returns the credible interval bounds for the fixed
#' parameters of the Bayesian model. The default posterior probability is 90%.
#'
#' The `bayesian_scoring` function returns the Bayesian scores.
#' These scores are the posterior probabilities
#' that the true latent \eqn{T_i}'s are equal to 1.
#'
#' The `bayesian_prevalence_estimate` function returns the posterior mean of the
#' posterior distribution on the prevalence of \eqn{T_i = 1}.
#'
#' @seealso \code{\link{classify_with_scores}}, \code{\link{BayesianFit}}
#'
#' @examples
#' data("periodontal")
#' theta <- mean(periodontal$ti)
#' fit <- BayesianFit(periodontal$ni, periodontal$si, chains = 2, iter = 5000)
#' credint(fit)
#' Y_B <- bayesian_scoring(fit)
#' T_B <- classify_with_scores(Y_B, .4, .6)
#' theta_B <- bayesian_prevalence_estimate(fit)
#' cat("The Bayesian prevalence estimate is ", theta_B, "\n")
#' cat("The prevalence in the data is ", theta, "\n")
#'
#' @rdname bayesian_computations
#' @importFrom stats quantile
#' @export
credint <- function(fit, level = .9) {
  # Extract the fixed parameters
  pars <- rstan::extract(fit, pars = c("p", "q", "theta"))
  # Compute the credible intervals
  alpha <- 1 - level
  out <- lapply(pars, function(x) quantile(x, c(alpha/2, 1-alpha/2)))
  out <- t(as.data.frame(out))
  out <- as.data.frame(out)
  colnames(out) <- c("lower", "upper")
  rownames(out) <- c("p", "q", "theta")
  return(out)
}


#' @rdname bayesian_computations
#' @export
bayesian_scoring <- function(fit){
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
