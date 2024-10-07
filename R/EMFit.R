#' Compute the maximum likelihood estimate with the EM algorithm
#'
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param ti Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual
#'  If \code{NULL}, the EM algorithm is used to estimate the parameters. Default to \code{NULL}. See details
#' @param vL Numeric, the lower bound for the classification of the scores. Default to 0.5
#' @param vU Numeric, the upper bound for the classification of the scores. Default to 0.5
#' @param N_init The number of initializations if \code{ti} is not provided. Default to 20
#' @param maxIter The maximum number of iterations if EM algorithm is used. Default to 1e3
#' @param errorMin The minimum error computed if EM algorithm is used. Default to 1e-7
#' @param correction Whether or not we want to correct the ML to avoid going to 0. Default to TRUE
#' @return A list with the following components:
#' \itemize{
#' \item{score}{ The estimated values of the scores}
#' \item{parameters_hat}{ The estimated values of the parameters \eqn{\theta}, \eqn{p} and \eqn{q}}
#' }
#'
#' @details
#' This function chooses its algorithm according to what is provided in the \code{ti} argument:
#' \itemize{
#'  \item{\emph{\code{ti} is fully provided}: }{the function computes the maximum likelihood
#'  estimate, with an explicit formula.}
#'  \item{\emph{\code{ti} is not provided}: }{the function uses the EM algorithm
#'  to estimate the parameters.}
#'  \item{\emph{\code{ti} is partially provided}: }{the function uses the EM algorithm to
#'  estimate the parameters.}
#' }
#'
#'
#' @examples
#' data("periodontal")
#' periodontal_ml <- EMFit(periodontal$si,periodontal$ni,periodontal$ti)
#' periodontal_EM <- EMFit(periodontal$si,periodontal$ni,ti = NULL)
#'
#'@seealso [classify_with_scores]
#'
#' @export
EMFit <- function(si,ni,ti=NULL,vL=0.5,vU=0.5,
                  N_init=20,maxIter=1e3,errorMin=1e-7,
                  correction=TRUE){
  ### Inner functions
  getLikelihood_MC <- function(ni,si,ti,p,q,theta,N=20){
    n <- length(ti)
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

  getLogLikelihood <- function(ni,si,p,q,theta){
    n <- length(ni)
    Likelihoods <- (theta*(1-q)^si*q^(ni-si) + (1-theta)*p^si*(1-p)^(ni-si))
    sum(log(Likelihoods))
  }

  # Fit the model
  emBin <- function(si,ni,ti=NULL,maxIter=1e3,errorMin=1e-7, correction=correction){
    n <- length(si)
    if(!is.null(ti) & !any(is.na(ti))){
      theta_hat <- mean(ti)
      p_hat <- mean(si[ti<1/2])/mean(ni[ti<1/2])
      q_hat <- mean((ni-si)[ti>=1/2])/mean(ni[ti>=1/2])
      out <- list(score=ti,ti_clas=ti,
                  parameters_hat=data.frame(list(theta=theta_hat,p=p_hat,q=q_hat)))
    } else {
      if(is.null(ti)){
        score <- rbeta(n,shape1 = si,shape2 = ni-si)#sample(c(0,1),n,replace=T)
        id_na <- 1:n
      } else {
        id_na <- which(is.na(ti))
        n_na <- length(id_na)
        score <- ti
        score[id_na] <- rbeta(n,shape1 = si,shape2 = ni-si)[id_na]#(si/ni)[id_na]#sample(c(0,1),n_na,replace=T)
      }
      iter <- 1
      error= 2*errorMin
      theta_hat <- mean(score)
      if (correction){
        # p_hat <- sum(si[score<1/2])/sum(ni[score<1/2])
        p_hat <- (sum(si*(1-score)) + 1)/(sum(ni*(1-score)) + 2)
        # q_hat <- sum((ni-si)[score>=1/2])/sum(ni[score>=1/2])
        q_hat <- (sum((ni-si)*score) + 1)/(sum(ni*score) + 2)

      } else {
        # p_hat <- sum(si[score<1/2])/sum(ni[score<1/2])
        p_hat <- sum(si*(1-score))/sum(ni*(1-score))
        # q_hat <- sum((ni-si)[score>=1/2])/sum(ni[score>=1/2])
        q_hat <- sum((ni-si)*score)/sum(ni*score)
      }
      ps <- p_hat
      while(error>errorMin & iter<maxIter){
        ## Estimation
        score_0 <- score
        numer <- theta_hat*q_hat^(ni-si)*(1-q_hat)^si
        denom_right <- (1-theta_hat)*p_hat^si*(1-p_hat)^(ni-si)
        denom <- numer+denom_right
        score[id_na] <- numer[id_na]/denom[id_na]
        if(any(is.na(score))) browser()
        ## Maximization
        theta_hat_0 <- theta_hat
        p_hat_0 <- p_hat
        q_hat_0 <- q_hat
        theta_hat <- mean(score)
        # id_left <- score<1/2
        # p_hat <- sum(si[id_left])/sum(ni[id_left])
        if (correction) {
          p_hat <- (sum(si*(1-score)) + 1)/(sum(ni*(1-score)) + 2)

        } else {
          p_hat <- sum(si*(1-score))/sum(ni*(1-score))
        }
        ps <- c(ps,p_hat)
        if(is.na(p_hat))
        {
          cat("p_hat NA")
          p_hat <- runif(1,0,1/2)
        }
        # id_right <- score>=1/2
        # q_hat <- sum((ni-si)[id_right])/sum(ni[id_right])
        if (correction) {
          q_hat <- (sum((ni-si)*score) + 1)/(sum(ni*score) + 2)
        } else {
          q_hat <- sum((ni-si)*score)/sum(ni*score)
        }
        if(is.na(q_hat))q_hat <- runif(1,0,1/2)
        if(p_hat>1/2 & q_hat>1/2){
          p_hat <- 1-p_hat
          q_hat <- 1-q_hat
          theta_hat <- 1-theta_hat
          score[id_na] <- 1-score[id_na]
        }
        error <- mean(abs(score_0[id_na]-score[id_na]))
        iter <- iter + 1
      }
      # plot(ps)
      # print(p_hat)
      # ti_clas <- (score<vL)*0 + (score>=vL & score<=vU)*(1/2) + (score>vU)*1
      out <- list(score=score,#ti_clas=ti_clas,
                  parameters_hat=data.frame(list(theta=theta_hat,p=p_hat,q=q_hat)))
    }
    out
  }
  ### ---------------

  if(!is.null(ti) & !any(is.na(ti))){
    out <- emBin(si,ni,ti)
  } else {
    # result <- matrix(NA,N_init,5)
    result <- matrix(NA,N_init,4)
    # colnames(result) <- c("theta","p","q","MCmean","MCerror")
    colnames(result) <- c("theta","p","q","Likelihood")
    allRes <- list()
    for(i in 1:N_init){
      allRes[[i]] <- emBin(si,ni)
      tii <- allRes[[i]]$ti_clas
      parai <- allRes[[i]]$parameters_hat
      # log_likeli <- getLikelihood_MC(ni,si,tii,
      #                             p=parai$p,q=parai$q,theta=parai$theta,
      #                             N=20)
      # result[i,] <- c(as.numeric(parai),log_likeli[1:2])
      loglikeli <- getLogLikelihood(ni,si,p=parai$p,q=parai$q,theta=parai$theta)
      result[i,] <- c(as.numeric(parai),loglikeli)
    }
    # i_opt <- which.max(result[,4])
    i_opt <- which.max(result[,3])
    out <- allRes[[i_opt]]
    # out$MC <- result
  }
  out
}
