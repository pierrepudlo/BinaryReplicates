#' Compute predictive Bayesian scores
#'
#' @param newdata_ni Numeric vector of the total numbers of replicates per individuals
#' @param newdata_si Numeric vector of the numbers of positive replicates per individuals
#' @param fit The `stanfit` object return by [BayesianFit]
#'
#' @details
#' The `predict_scores` function computes the predictive Bayesian scores.
#' It makes the empirical estimator, for a new individual \eqn{n+1}$, of the following integral:
#' \deqn{Y_{B,n+1} = \int Y_{L,n+1}(\theta_T, p, q) \pi(\theta_T, p, q|S_1,...,S_{n})\text{d}\theta_T\text{d}p\text{d}q}
#' where \eqn{\pi(\theta, p, q|S_1,...,S_{n})} is the posterior distribution
#' of the parameters \eqn{\theta}, \eqn{p} and \eqn{q} given the data
#' \eqn{S_1,...,S_{n}} and \eqn{Y_{L,n+1}(\theta_T, p, q)} is given
#' by the function [likelihood_scoring], such as
#' \deqn{Y_{L,n+1}(\theta_T, p, q) = \boldsymbol{P}(T_{n+1}=1|S_{n+1}=s_{n+1}) = \frac{\theta_T q^{n_{n+1}-S_{n+1}} {(1-q)}^{S_{n+1}}}{\theta_T q^{n_{n+1}-S_{n+1}} (1-q)^{S_{n+1}} + (1-\theta_T)p^{S_{n+1}} {(1-p)}^{n_{n+1}-S_{n+1}}}.}
#' Thus the estimator is given by
#' \deqn{\hat{Y}_{B,n+1} = \frac{1}{K} \sum_{k=1}^K Y_{L,n+1}(\theta_{T,k}, p_k, q_k),}
#' where each parameter \eqn{(\theta_{T,k}, p_k, q_k)_k} is sampled from the
#' posterior distribution, output of the function [BayesianFit]. \eqn{K} is the total number of sampled parameters.
#' @importFrom rstan extract
#' @importFrom stats aggregate
#' @return
#' The `predict_scores` function returns the predictive Bayesian scores in a numeric vector.
#' The predictive Bayesian scores are the posterior probabilities that the true
#' latent \eqn{T_i}'s are equal to 1 on new data, averaged over the posterior distribution.
#' @export
#' @examples
#' data("periodontal")
#' theta <- mean(periodontal$ti)
#' fitBay <- BayesianFit(periodontal$ni, periodontal$si, chains = 2, iter = 500)
#' fitMAP <- EMFit(periodontal$si,periodontal$ni)
#'
#' ## Comparaison Bayesian <--> MAP
#' ni <- 200
#' Ni <- rep(ni,ni+1)
#' Si <- 0:ni
#' scores <- cbind(predict_scores(Ni,Si,fitBay),
#'                 likelihood_scoring(Ni,Si,fitMAP$parameters_hat))
#' matplot(Si,scores,type = "l",lty = 1,col = 1:2,
#'         ylab = "Scores",xlab = "Number of Successes",main = "")
predict_scores <- function(newdata_ni, newdata_si, fit) {
  if (length(newdata_ni) != length(newdata_si)) {
    stop("newdata_ni and newdata_si must have the same length")
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

  # Calculate Y_B
  newdata$Y_B <- likelihood_scoring(newdata$ni, newdata$si,
                                    list(theta = newdata$theta,
                                         p = newdata$p,
                                         q = newdata$q))

  # Calculate mean Y_B by id using base R
  out <- aggregate(Y_B ~ id, data = newdata, FUN = mean)
  return(out$Y_B)
}
