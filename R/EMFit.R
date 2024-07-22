getLogLikelihood <- function(ni, si, p, q, theta, N = 20) {
  sum(log(theta * dbinom(si, ni, 1 - q) + (1 - theta) * dbinom(si, ni, p)))
}
# Fit the model
emBin <- function(ni, si, ti = NULL, maxIter = 1e3, errorMin = 1e-7) {
  n <- length(si)
  if (!is.null(ti) & !any(is.na(ti))) {
    theta_hat <- mean(ti)
    p_hat <- mean(si[ti < .5]) / mean(ni[ti < .5])
    q_hat <- mean((ni - si)[ti >= .5]) / mean(ni[ti >= .5])
    out <- list(
      ti_est = ti,
      ti_clas = ti,
      parameters_hat = data.frame(list(
        theta = theta_hat, p = p_hat, q = q_hat
      ))
    )
  } else {
    if (is.null(ti)) {
      ti_est <- rbeta(n, shape1 = si + .5, shape2 = ni - si + .5)
      id_na <- 1:n
    } else {
      id_na <- which(is.na(ti))
      n_na <- length(id_na)
      ti_est <- ti
      ti_est[id_na] <- rbeta(n, shape1 = si + .5, shape2 = ni - si + .5)[id_na]
    }
    iter <- 1
    error <- 2 * errorMin
    theta_hat <- mean(ti_est)
    p_hat <- mean(si[ti_est < .5]) / mean(ni[ti_est < .5])
    q_hat <- mean((ni - si)[ti_est >= .5]) / mean(ni[ti_est >= .5])
    while (error > errorMin & iter < maxIter) {
      ## Estimation
      ti_est_0 <- ti_est
      numer <- theta_hat * q_hat^(ni - si) * (1 - q_hat)^si
      denom_right <- (1 - theta_hat) * p_hat ^ si * (1 - p_hat)^(ni - si)
      denom <- numer + denom_right
      ti_est[id_na] <- numer[id_na] / denom[id_na]
      if (any(is.na(ti_est)))
        browser()
      ## Maximization
      theta_hat_0 <- theta_hat
      p_hat_0 <- p_hat
      q_hat_0 <- q_hat
      theta_hat <- mean(ti_est)
      id_left <- ti_est < .5
      p_hat <- mean(si[id_left]) / mean(ni[id_left])
      if (is.na(p_hat))
      {
        cat("p_hat NA")
        p_hat <- runif(1, 0, .5)
      }
      id_right <- ti_est >= .5
      q_hat <- mean((ni - si)[id_right]) / mean(ni[id_right])
      if (is.na(q_hat))
        q_hat <- runif(1, 0, .5)
      if (p_hat > .5 & q_hat > .5) {
        p_hat <- 1 - q_hat
        q_hat <- 1 - p_hat
        theta_hat <- 1 - theta_hat
        ti_est[id_na] <- 1 - ti_est[id_na]
      }
      error <- mean(abs(ti_est_0[id_na] - ti_est[id_na]))
      iter <- iter + 1
    }
    ti_clas <- (ti_est > .5) * 1
    out <- list(
      ti_est = ti_est,
      ti_clas = ti_clas,
      parameters_hat = data.frame(list(
        theta = theta_hat, p = p_hat, q = q_hat
      ))
    )
  }
  out
}


#' Compute the maximum likelihood estimate with the EM algorithm
#'
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param ti Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual
#'  If \code{NULL}, the EM algorithm is used to estimate the parameters. Default to \code{NULL}. See details
#' @param NN The number of runs of the EM algorithm if \code{ti} is not provided. Default to 20
#' @param maxIter The maximum number of iterations if EM algorithm is used. Default to 1e3
#' @param errorMin The minimum error computed if EM algorithm is used. Default to 1e-7
#' @return A list with the following components:
#' \itemize{
#' \item{ti_est}{The estimated values of \eqn{t_i}}
#' \item{ti_clas}{The estimated values of \eqn{t_i} classified as 0 or 1}
#' \item{parameters_hat}{The estimated values of the parameters \eqn{\theta}, \eqn{p} and \eqn{q}}
#' \item{MC}{If \code{ti} is not provided, a matrix with the estimated values of the parameters and the Monte Carlo mean and error}
#' }
#'
#' @details
#' This function chooses its algorithm according to what is provided in the \code{ti} argument:
#' \itemize{
#'  \item{\code{ti} is fully provided}{the function computes the maximum likelihood
#'  estimate, with an explicit formula.}
#'  \item{\code{ti} is not provided}{the function uses the EM algorithm
#'  to estimate the parameters.}
#'  \item{\code{ti} is partially provided}{the function uses the EM algorithm to
#'  estimate the parameters.}
#' }
#'
#'
#' @examples
#' data("periodontal")
#' periodontal_ml <- EMFit(periodontal$si,periodontal$ni,periodontal$ti)
#' periodontal_EM <- EMFit(periodontal$si,periodontal$ni,ti = NULL)
#'
#' @export
EMFit <- function(ni, si, ti = NULL, NN = 20, maxIter = 1e3, errorMin = 1e-7) {
  if (!is.null(ti) & !any(is.na(ti))) {
    out <- emBin(ni, si, ti)
  } else {
    result <- data.frame(
      theta = rep(0, NN),
      p = rep(0, NN),
      q = rep(0, NN),
      logLik = rep(0, NN)
    )
    allRes <- list()
    for (i in 1:NN) {
      allRes[[i]] <- emBin(ni, si)
      tii <- allRes[[i]]$ti_clas
      parai <- allRes[[i]]$parameters_hat
      log_likeli <- getLogLikelihood(ni,
                                     si,
                                     p = parai$p,
                                     q = parai$q,
                                     theta = parai$theta)
      result[i, ] <- c(as.numeric(parai), log_likeli)
    }
    i_opt <- which.max(result$logLik)
    out <- c(theta = result$theta[i_opt],
             p = result$p[i_opt],
             q = result$q[i_opt])
  }
  out
}
