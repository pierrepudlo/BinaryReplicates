#' Get credible intervals, Bayesian scores, and prevalence estimate
#' from the `stanfit` object returned by [BayesianFit].
#'
#' @name bayesian_computations
#' @title Bayesian computations
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param fit The `stanfit` object returned by [BayesianFit]
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
#' Y_B <- bayesian_scoring(periodontal$ni, periodontal$si, fit)
#' T_B <- classify_with_scores(Y_B, .4, .6)
#' theta_B <- bayesian_prevalence_estimate(fit)
#' cat("The Bayesian prevalence estimate is ", theta_B, "\n")
#' cat("The prevalence in the data is ", theta, "\n")
#'
#' @rdname bayesian_computations
#' @importFrom stats quantile
#' @importFrom rstan extract
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
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#' @importFrom magrittr %>%
#' @export
bayesian_scoring <- function(ni, si, fit){
  # Extract the latent variable
  pi <- rstan::extract(fit, pars = "posterior_prob")$posterior_prob
  # Compute the posterior probability
  yi <- colMeans(pi)
  # Pairs (ni, si) that are equal should have the same posterior probability
  # We need to average the Monte Carlo estimates with a group_by
  results <- data.frame(ni = ni, si = si, yi = yi) %>%
    dplyr::group_by(ni, si) %>%
    dplyr::summarise(yi = mean(yi), .groups="keep")
  yi2 <- data.frame(ni = ni, si = si) %>%
    dplyr::left_join(results, by = c("ni", "si")) %>%
    dplyr::pull(yi)
  return(yi2)
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
