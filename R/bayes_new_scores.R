#' Compute predictive Bayesian scores
#'
#' @param fit The `stanfit` object return by [BayesianFit]
#' @param newdata_ni Numeric vector of the total numbers of replicates per individuals
#' @param newdata_si Numeric vector of the numbers of positive replicates per individuals
#'
#' @importFrom rstan extract
#' @importFrom dplyr summarise group_by ungroup
#' @return
#' The `predict_scores` function returns the predictive Bayesian scores in a numeric vector.
#' The predictive Bayesian scores are the posterior probabilities that the true
#' latent \eqn{T_i}'s are equal to 1 on new data, averaged over the posterior distribution.
#' @export
#' @examples
#' data("periodontal")
#' theta <- mean(periodontal$ti)
#' fit <- BayesianFit(periodontal$ni, periodontal$si, chains = 2, iter = 5000)
#' credint(fit)
#' newdata <- data.frame(ni = rep(4, 5), si = 0:4)
#' newdata$Y_BP <- predict_scores(fit, newdata$ni, newdata$si)
#' newdata
predict_scores <- function(fit, newdata_ni, newdata_si) {
  if (length(newdata_ni) != length(newdata_si)) {
    stop("newdata_ni and newdata_ni must have the same length")
  }
  n <- length(newdata_ni)
  # Extract the latent variable
  post <- rstan::extract(fit, pars = c("theta", "p", "q"))
  n_post <- length(post$theta)
  # Compute the posterior probability
  newdata <- data.frame(
    id = rep(1:n, each = n_post),
    ni = rep(newdata_ni, each = n_post),
    si = rep(newdata_si, each = n_post),
    theta = rep(post$theta, times = n),
    p = rep(post$p, times = n),
    q = rep(post$q, times = n)
  )
  newdata$Y_B <- likelihood_scoring(newdata$ni, newdata$si, newdata$theta,
                                    newdata$p, newdata$q)
  out <- dplyr::summarise(dplyr::group_by(newdata, newdata$id,
                                          newdata$ni, newdata$si),
                          Y_B = mean(newdata$Y_B))
  out <- dplyr::ungroup(out)

  return(out$Y_B)
}
