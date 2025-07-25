% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes_new_scores.R
\name{predict_scores}
\alias{predict_scores}
\title{Compute predictive Bayesian scores}
\usage{
predict_scores(newdata_ni, newdata_si, fit)
}
\arguments{
\item{newdata_ni}{Numeric vector of the total numbers of replicates per individuals}

\item{newdata_si}{Numeric vector of the numbers of positive replicates per individuals}

\item{fit}{The \code{stanfit} object return by \link{BayesianFit}}
}
\value{
The \code{predict_scores} function returns the predictive Bayesian scores in a numeric vector.
The predictive Bayesian scores are the posterior probabilities that the true
latent \eqn{T_i}'s are equal to 1 on new data, averaged over the posterior distribution.
}
\description{
Compute predictive Bayesian scores
}
\details{
The \code{predict_scores} function computes the predictive Bayesian scores.
It makes the empirical estimator, for a new individual \eqn{n+1}$, of the following integral:
\deqn{Y_{B,n+1} = \int Y_{L,n+1}(\theta_T, p, q) \pi(\theta_T, p, q|S_1,...,S_{n})\text{d}\theta_T\text{d}p\text{d}q}
where \eqn{\pi(\theta, p, q|S_1,...,S_{n})} is the posterior distribution
of the parameters \eqn{\theta}, \eqn{p} and \eqn{q} given the data
\eqn{S_1,...,S_{n}} and \eqn{Y_{L,n+1}(\theta_T, p, q)} is given
by the function \link{likelihood_scoring}, such as
\deqn{Y_{L,n+1}(\theta_T, p, q) = \boldsymbol{P}(T_{n+1}=1|S_{n+1}=s_{n+1}) = \frac{\theta_T q^{n_{n+1}-S_{n+1}} {(1-q)}^{S_{n+1}}}{\theta_T q^{n_{n+1}-S_{n+1}} (1-q)^{S_{n+1}} + (1-\theta_T)p^{S_{n+1}} {(1-p)}^{n_{n+1}-S_{n+1}}}.}
Thus the estimator is given by
\deqn{\hat{Y}_{B,n+1} = \frac{1}{K} \sum_{k=1}^K Y_{L,n+1}(\theta_{T,k}, p_k, q_k),}
where each parameter \eqn{(\theta_{T,k}, p_k, q_k)_k} is sampled from the
posterior distribution, output of the function \link{BayesianFit}. \eqn{K} is the total number of sampled parameters.
}
\examples{
data("periodontal")
theta <- mean(periodontal$ti)
fitBay <- BayesianFit(periodontal$ni, periodontal$si, chains = 2, iter = 500)
fitMAP <- EMFit(periodontal$si,periodontal$ni)

## Comparaison Bayesian <--> MAP
ni <- 200
Ni <- rep(ni,ni+1)
Si <- 0:ni
scores <- cbind(predict_scores(Ni,Si,fitBay),
                likelihood_scoring(Ni,Si,fitMAP$parameters_hat))
matplot(Si,scores,type = "l",lty = 1,col = 1:2,
        ylab = "Scores",xlab = "Number of Successes",main = "")
}
