# Save this file as `R/BayesianFit.R`

#' Fit the Bayesian model for Binary Replicates
#'
#' @export
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param prior A list of prior parameters for the model, see Details.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @details The model is a Bayesian model for binary replicates.
#' The prior distribution is as follows:
#'
#' \itemize{
#' \item The false positivity rate: \eqn{p \sim \text{Beta}(a_{FP}, b_{FP})}
#' \item The false negativity rate: \eqn{q \sim \text{Beta}(a_{FN}, b_{FN})}
#' \item The prevalence: \eqn{\theta \sim \text{Beta}(a_T, b_T)}
#' }
#'
#' The statistical model considers that the true status are the latent
#' \deqn{T_i \sim \text{Bernoulli}(\theta).}
#'
#' And, given the true status \eqn{T_i}, the number of positive replicates is
#' \deqn{S_i \sim \text{Binomial}(n_i, T_i(1-q)+(1-T_i)p).}
#' @return An object of class `stanfit` returned by `rstan::sampling`.
#'
#' The Stan model samples the posterior distribution of the fixed parameters \eqn{p},
#' \eqn{q} and \eqn{\theta}. It also generates the latent variables \eqn{T_i} according
#' to their predictive distribution.
#'
#' @seealso [credint], [bayesian_scoring], [classify_with_scores],
#' [bayesian_prevalence_estimate]
BayesianFit <- function(ni, si,
                    prior = list(a_FP=2, b_FP=2,
                                 a_FN=2, b_FN=2,
                                 a_T=.5, b_T=.5),...) {
  # Check the data
  if (!is.numeric(ni) | !is.numeric(si)) {
    stop("ni and si must be numeric vectors")
  }
  if (length(ni) != length(si)) {
    stop("ni and si must have the same length")
  }
  if (!is.integer(ni)) ni <- as.integer(ni)
  if (!is.integer(si)) si <- as.integer(si)

  n <- length(ni)
  standata <- list(n = n, ni = ni, si = si,
                   a_FP = prior$a_FP, b_FP = prior$b_FP,
                   a_FN = prior$a_FN, b_FN = prior$b_FN,
                   a_T = prior$a_T, b_T = prior$b_T)
  out <- rstan::sampling(stanmodels$BayesianModel, data = standata, ...)
  return(out)
}
