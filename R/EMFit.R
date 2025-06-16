#' Compute the *Maximum-A-Posteriori* estimate with the EM algorithm
#'
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param ti Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual
#'  If \code{NULL}, the EM algorithm is used to estimate the parameters. Default to \code{NULL}. See details
#' @param N_init The number of initializations if \code{ti} is not provided. Default to 20
#' @param maxIter The maximum number of iterations if EM algorithm is used. Default to 1e3
#' @param errorMin The minimum error computed if EM algorithm is used. Default to 1e-7
#' @param prior A list of prior parameters for the model. The prior distribution is as follows:
#'
#' \itemize{
#' \item The false positivity rate: \eqn{p \sim \text{Beta}(a_{FP}, b_{FP})}
#' \item The false negativity rate: \eqn{q \sim \text{Beta}(a_{FN}, b_{FN})}
#' }
#' @return A list with the following components:
#' \itemize{
#' \item{score}{ The estimated values of the scores}
#' \item{parameters_hat}{ The estimated values of the parameters \eqn{\theta}, \eqn{p} and \eqn{q}}
#' }
#'
#' @details
#' This function chooses its algorithm according to what is provided in the \code{ti} argument:
#' \itemize{
#'  \item{\emph{\code{ti} is fully provided}: }{the function computes the *Maximum-A-Posteriori*
#'  estimate, with an explicit formula.}
#'  \item{\emph{\code{ti} is not provided}: }{the function uses the EM algorithm
#'  to estimate the parameters.}
#'  \item{\emph{\code{ti} is partially provided}: }{the function uses the EM algorithm to
#'  estimate the parameters.}
#' }
#'
#' @importFrom stats rbinom dbinom
#'
#' @examples
#' data("periodontal")
#' periodontal_ml <- EMFit(periodontal$si,periodontal$ni,periodontal$ti)
#' periodontal_EM <- EMFit(periodontal$si,periodontal$ni,ti = NULL)
#'
#'@seealso [classify_with_scores]
#'
#' @export
EMFit <- function(si,ni,ti=NULL,prior=list(a_FP=1, b_FP=1,
                                           a_FN=1, b_FN=1),
                  N_init=20,maxIter=1e3,errorMin=1e-7){
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
  emBin <- function(si,ni,ti=NULL,maxIter=1e3,errorMin=1e-7, prior=prior){
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
      # if (correction){
      #   # p_hat <- sum(si[score<1/2])/sum(ni[score<1/2])
      #   # p_hat <- (sum(si*(1-score)) + 1)/(sum(ni*(1-score)) + 2)
      #   p_hat <- (sum(si*(1-score)) + 1)/(sum(ni*(1-score)) + 6)
      #   # q_hat <- sum((ni-si)[score>=1/2])/sum(ni[score>=1/2])
      #   # q_hat <- (sum((ni-si)*score) + 1)/(sum(ni*score) + 2)
      #   q_hat <- (sum((ni-si)*score) + 1)/(sum(ni*score) + 6)
      # } else {
      #   # p_hat <- sum(si[score<1/2])/sum(ni[score<1/2])
      #   p_hat <- sum(si*(1-score))/sum(ni*(1-score))
      #   # q_hat <- sum((ni-si)[score>=1/2])/sum(ni[score>=1/2])
      #   q_hat <- sum((ni-si)*score)/sum(ni*score)
      # }
      p_hat <- (sum(si*(1-score)) + prior$a_FP)/
        (sum(ni*(1-score)) + prior$a_FP + prior$b_FP)
      q_hat <- (sum((ni-si)*score) + prior$a_FN)/
        (sum(ni*score) + prior$a_FN + prior$b_FN)
      # p_hat=q_hat <- 1/2
      p_hat <- rbeta(1,2,2)
      q_hat <- rbeta(1,2,2)
      ps <- p_hat
      qs <- q_hat
      thetas <- theta_hat
      while(error>errorMin & iter<maxIter){
        ## Estimation
        score_0 <- score
        numer <- theta_hat*q_hat^(ni-si)*(1-q_hat)^si
        denom_right <- (1-theta_hat)*p_hat^si*(1-p_hat)^(ni-si)
        denom <- numer+denom_right
        score[id_na] <- numer[id_na]/denom[id_na]
        ## Maximization
        theta_hat_0 <- theta_hat
        p_hat_0 <- p_hat
        q_hat_0 <- q_hat
        theta_hat <- mean(score)
        # if (correction) {
        #   # p_hat <- (sum(si*(1-score)) + 1)/(sum(ni*(1-score)) + 2)
        #   p_hat <- (sum(si*(1-score)) + 1)/(sum(ni*(1-score)) + 6)
        # } else {
        #   p_hat <- sum(si*(1-score))/sum(ni*(1-score))
        # }
        p_hat <- (sum(si*(1-score)) + prior$a_FP)/
          (sum(ni*(1-score)) + prior$a_FP + prior$b_FP)
        q_hat <- (sum((ni-si)*score) + prior$a_FN)/
          (sum(ni*score) + prior$a_FN + prior$b_FN)
        ps <- c(ps,p_hat)
        qs <- c(qs,q_hat)
        thetas <- c(thetas,theta_hat)
        # if (correction) {
        #   # q_hat <- (sum((ni-si)*score) + 1)/(sum(ni*score) + 2)
        #   q_hat <- (sum((ni-si)*score) + 1)/(sum(ni*score) + 6)
        # } else {
        #   q_hat <- sum((ni-si)*score)/sum(ni*score)
        # }
        p_hat <- (sum(si*(1-score)) + prior$a_FP)/
          (sum(ni*(1-score)) + prior$a_FP + prior$b_FP)
        q_hat <- (sum((ni-si)*score) + prior$a_FN)/
          (sum(ni*score) + prior$a_FN + prior$b_FN)
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
      out <- list(score=score,
                  parameters_hat=data.frame(list(theta=theta_hat,p=p_hat,q=q_hat)),
                  all_paras=list(theta=thetas,p=ps,q=qs))
    }
    out
  }
  ### ---------------

  if(!is.null(ti) & !any(is.na(ti))){
    out <- emBin(si,ni,ti,prior=prior,
                 errorMin = errorMin,maxIter = maxIter)
  } else {
    # result <- matrix(NA,N_init,5)
    result <- matrix(NA,N_init,4)
    # colnames(result) <- c("theta","p","q","MCmean","MCerror")
    colnames(result) <- c("theta","p","q","LogPosterior")
    allRes <- list()
    all_initializations <- list()
    for(i in 1:N_init){
      allRes[[i]] <- emBin(si,ni,prior=prior,
                           errorMin = errorMin,maxIter = maxIter)
      tii <- allRes[[i]]$ti_clas
      parai <- allRes[[i]]$parameters_hat
      loglikeli <- getLogLikelihood(ni,si,p=parai$p,q=parai$q,theta=parai$theta)
      logPoster <- loglikeli + log(parai$p*(1-parai$p)*parai$q*(1-parai$q))
      result[i,] <- c(as.numeric(parai),logPoster)#loglikeli)
      allRes[[i]]$all_paras$logPosterior <- logPoster
      all_initializations[[i]] <- allRes[[i]]$all_paras
    }
    # i_opt <- which.max(result[,4])
    i_opt <- which.max(result[,3])
    out <- allRes[[i_opt]]
    out$all_paras <- NULL
    out$all_initializations <- all_initializations
  }
  out
}

#' Function to perform cross-validation in the EM algorithm
#'
#' @param ni Numeric vector of \eqn{n_i}'s, the total numbers of replicates for each individual
#' @param si Numeric vector of \eqn{s_i}'s, the numbers of replicates equal to 1 for each individual
#' @param ti Numeric vector of \eqn{t_i}'s, the true values of the binary variable for each individual
#'  If \code{NULL}, the EM algorithm is used to estimate the parameters. Default to \code{NULL}. See details
#' @param N_cv The number of folds. Default to 20
#' @param N_init The number of initializations if \code{ti} is not provided. Def
#' ault to 20 if  observations otherwise default to the number of observations,
#' corresponding to leave-one-out cross-validation.
#' @param maxIter The maximum number of iterations if EM algorithm is used. Default to 1e3
#' @param errorMin The minimum error computed if EM algorithm is used. Default to 1e-7
#' @param prior A list of prior parameters for the model. The prior distribution is as follows:
#'
#' \itemize{
#' \item The false positivity rate: \eqn{p \sim \text{Beta}(a_{FP}, b_{FP})}
#' \item The false negativity rate: \eqn{q \sim \text{Beta}(a_{FN}, b_{FN})}
#' }
#'
#' @return A list with the following components:
#' \itemize{
#' \item{score}{ The estimated values of the scores}
#' \item{parameters_hat}{ The estimated values of the parameters \eqn{\theta}, \eqn{p} and \eqn{q}}
#' }
#'
#' @details
#' This function chooses its algorithm according to what is provided in the \code{ti} argument:
#' \itemize{
#'  \item{\emph{\code{ti} is fully provided}: }{the function computes the *Maximum-A-Posteriori*
#'  estimate, with an explicit formula.}
#'  \item{\emph{\code{ti} is not provided}: }{the function uses the EM algorithm
#'  to estimate the parameters.}
#'  \item{\emph{\code{ti} is partially provided}: }{the function uses the EM algorithm to
#'  estimate the parameters.}
#' }
#'
#'
#' @return A list with the following components:
#' \itemize{
#' \item{models}{ A list of the models for each fold}
#' \item{predictions}{ A list of the predictions for each fold}
#' }
#' @export
#'
#' @examples
#' data("periodontal")
#' modelCV <- cvEM(periodontal$ni,periodontal$si)
#'
#'@seealso [classify_with_scores,EMFit]
cvEM <- function(ni,si,ti=NULL,N_cv=NULL,
                 N_init=20,maxIter=1e3,errorMin=1e-7,
                 prior=list(a_FP=1, b_FP=1,
                            a_FN=1, b_FN=1)){
  n <- length(ni)
  if(is.null(N_cv)) N_cv <- min(20,n)
  if(N_cv < 2)stop("Choose at least 2 folds")
  if(N_cv%%1 != 0)stop("Select an entire number of folds")
  if(N_cv > n)stop("Select at most n folds")
  id_cv <- unlist(lapply(1:N_cv,rep,floor(n/N_cv)))
  if(n%%N_cv>0){
    id_cv <- c(id_cv,1:(n%%N_cv))
  }
  id_cv <- sample(id_cv,size = n,replace = FALSE)
  MODELS <- list()
  PREDICTIONS <- list()
  for(i in 1:N_cv){
    id_test <- id_cv==i
    id_train <- !id_test
    ni_train <- ni[id_train]
    si_train <- si[id_train]
    ni_test <- ni[id_test]
    si_test <- si[id_test]
    if(!is.null(ti)){
      ti_train <- ti[id_train]
      ti_test <- ti[id_test]
    } else {
      ti_train <- NULL
      ti_test <- NULL
    }
    MODELS[[i]] <- EMFit(si = si_train,ni = ni_train,ti = ti_train,
                         N_init=N_init,maxIter=maxIter,errorMin=errorMin,
                         prior=prior)
    PREDICTIONS[[i]] <- list()
    PREDICTIONS[[i]]$scores_predicted <- likelihood_scoring(ni_test,si_test,MODELS[[i]]$parameters_hat)
    names(PREDICTIONS[[i]]$scores_predicted) <- which(id_test)
    if(!is.null(ti)){
      PREDICTIONS[[i]]$class_observed <- ti_test
      names(PREDICTIONS[[i]]$class_observed) <- which(id_test)
    } else {
      PREDICTIONS[[i]]$class_observed <- NULL
    }
  }
  out <- list(models=MODELS,predictions=PREDICTIONS)
  class(out) <- "cvEM"
  out
}

#' Perform classification on the scores for each fold of a cvEM object.
#' Note that if \eqn{t_i} is provided, the empirical risk is estimated with
#' \eqn{a=v_L}.
#'
#' @param object An object of class cvEM
#' @param vL The lower threshold for classification. Default to 0.5
#' @param vU The upper threshold for classification. Default to 0.5
#'
#' @return A cvEM object with the following components:
#' \itemize{
#' \item{predictions}{ A list of the predictions for each fold}
#' \item{risk}{ The empirical risk if \eqn{t_i} is provided}
#' }
#' @export
#'
#' @examples
#' data(periodontal)
#' modelCV <- cvEM(periodontal$ni,periodontal$si)
#' modelCV2 <- classify_with_scores_cvEM(modelCV,vL=0.4)
classify_with_scores_cvEM <- function(object,ti=NULL,vL=0.5,vU=0.5){
  get_risk <- function(ti_hat,ti,a){
    n <- length(ti)
    risk <- 0
    for(i in 1:n){
      if(abs(ti[i]-ti_hat[i])==1 ) risk <- risk + 1
      if(abs(ti[i]-ti_hat[i])==1/2 ) risk <- risk + a
    }
    risk
  }
  if(!inherits(object,"cvEM"))stop("The object must be of class cvEM")
  N_cv <- length(object$predictions)
  if(!is.null(ti)) RISK <- 0
  for(i in 1:N_cv){
    object$predictions[[i]]$classes_predicted <- classify_with_scores(object$predictions[[i]]$scores_predicted,vL=vL,vU=vU)
    id_i <- as.numeric(names(object$predictions[[i]]$scores_predicted))
    names(object$predictions[[i]]$classes_predicted) <-  id_i
    if(!is.null(ti)){
      RISK <- RISK + get_risk(object$predictions[[i]]$classes_predicted,ti[id_i],a=vL)
    }
  }
  if(!is.null(ti)) object$risk <- RISK
  object
}
