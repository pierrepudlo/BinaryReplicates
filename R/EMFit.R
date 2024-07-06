#' Compute the maximum likelihood estimate with the EM algorithm
#'
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param ti Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual
#'  If \code{NULL}, the EM algorithm is used to estimate the parameters. Default to \code{NULL}. See details
#' @param NN The number of Monte Carlo samples if \code{ti} is not provided. Default to 20
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
EMFit <- function(si,ni,ti=NULL,NN=20,maxIter=1e3,errorMin=1e-7){
  ### Inner functions
  getLikelihood <- function(ni,si,ti,p,q,theta,N=20){
    n <- length(ti)
    #L <- BinaryReplicates::likelihood_scoring(ni, si, theta, p, q)
    L <- rep(theta, n)
    likelihood_rand <- rep(NA,N)
    for(t in 1:N){
      rand_ti <- rbinom(n,1,L)
      denom <- sum(dbinom(rand_ti,1,L,log=TRUE))
      numer <- sum(dbinom(si,ni,ti*(1-q)+(1-ti)*p,log = TRUE))
      numer <- numer +
        sum(dbinom(rand_ti,1,theta,log=TRUE))
      likelihood_rand[t] <- numer - denom
    }
    maxSample <- max(likelihood_rand)
    likelihood_rand <- exp(likelihood_rand-maxSample)
    log_likelihood_meanMC <- log(mean(likelihood_rand)) + maxSample
    log_likelihood_sdMC <- log(sd(likelihood_rand)) + maxSample
    c(MCmean=log_likelihood_meanMC,MCerror=log_likelihood_sdMC-log(N)/2)
  }

  # Fit the model
  emBin <- function(si,ni,ti=NULL,maxIter=1e3,errorMin=1e-7){
    n <- length(si)
    if(!is.null(ti) & !any(is.na(ti))){
      theta_hat <- mean(ti)
      p_hat <- mean(si[ti<1/2])/mean(ni[ti<1/2])
      q_hat <- mean((ni-si)[ti>=1/2])/mean(ni[ti>=1/2])
      out <- list(ti_est=ti,ti_clas=ti,
                  parameters_hat=data.frame(list(theta=theta_hat,p=p_hat,q=q_hat)))
    } else {
      if(is.null(ti)){
        ti_est <- rbeta(n,shape1 = si,shape2 = ni-si)#sample(c(0,1),n,replace=T)
        id_na <- 1:n
      } else {
        id_na <- which(is.na(ti))
        n_na <- length(id_na)
        ti_est <- ti
        ti_est[id_na] <- rbeta(n,shape1 = si,shape2 = ni-si)[id_na]#(si/ni)[id_na]#sample(c(0,1),n_na,replace=T)
      }
      iter <- 1
      error= 2*errorMin
      theta_hat <- mean(ti_est)
      p_hat <- mean(si[ti_est<1/2])/mean(ni[ti_est<1/2])
      q_hat <- mean((ni-si)[ti_est>=1/2])/mean(ni[ti_est>=1/2])
      while(error>errorMin & iter<maxIter){
        ## Estimation
        ti_est_0 <- ti_est
        numer <- theta_hat*q_hat^(ni-si)*(1-q_hat)^si
        denom_right <- (1-theta_hat)*p_hat^si*(1-p_hat)^(ni-si)
        denom <- numer+denom_right
        ti_est[id_na] <- numer[id_na]/denom[id_na]
        if(any(is.na(ti_est))) browser()
        ## Maximization
        theta_hat_0 <- theta_hat
        p_hat_0 <- p_hat
        q_hat_0 <- q_hat
        theta_hat <- mean(ti_est)
        id_left <- ti_est<1/2
        # if(length(id_left) != 0){
        p_hat <- mean(si[id_left])/mean(ni[id_left])
        # } else {
        if(is.na(p_hat)) p_hat <- runif(1,0,1/2)
        # }
        id_right <- ti_est>=1/2
        # if(length(id_right) != 0){
        q_hat <- mean((ni-si)[id_right])/mean(ni[id_right])
        # } else {
        if(is.na(q_hat))q_hat <- runif(1,0,1/2)
        # }
        if(is.na(p_hat)) browser()
        if(p_hat>1/2 & q_hat>1/2){
          p_hat <- 1-p_hat
          q_hat <- 1-q_hat
          theta_hat <- 1-theta_hat
          ti_est[id_na] <- 1-ti_est[id_na]
        }
        error <- mean(abs(ti_est_0[id_na]-ti_est[id_na]))
        iter <- iter + 1
      }
      ti_clas <- (ti_est>1/2)*1
      out <- list(ti_est=ti_est,ti_clas=ti_clas,
                  parameters_hat=data.frame(list(theta=theta_hat,p=p_hat,q=q_hat)))
    }
    out
  }
  ### ---------------

  if(!is.null(ti) & !any(is.na(ti))){
    out <- emBin(si,ni,ti)
  } else {
    result <- matrix(NA,NN,5)
    colnames(result) <- c("theta","p","q","MCmean","MCerror")
    allRes <- list()
    for(i in 1:NN){
      allRes[[i]] <- emBin(si,ni)
      tii <- allRes[[i]]$ti_clas
      parai <- allRes[[i]]$parameters_hat
      log_likeli <- getLikelihood(ni,si,tii,
                                  p=parai$p,q=parai$q,theta=parai$theta,
                                  N=20)
      result[i,] <- c(as.numeric(parai),log_likeli[1:2])
    }
    i_opt <- which.max(result[,4])
    out <- allRes[[i_opt]]
    out$MC <- result
  }
  out
}
