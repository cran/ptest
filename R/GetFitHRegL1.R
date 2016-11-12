#' Compute loglikelihood ratio test statistic under the laplace noises
#' 
#' @description  
#' The -2 loglikelihood ratio test statistic, -2LLR, is computed 
#' for testing for periodicity under the laplace noises.
#' 
#' @param y Vector containing the series.
#'
#' @return  Vector of length 2 containing 
#' the -2 loglikelihood ratio test statistic under the laplace noises 
#' and the estimated frequency.
#' 
#' @details 
#' To compute the likelihood ratio statistic, 
#' a harmonic regression model is fitted by 
#' the maximum likelihood estimations (MLE)
#' according to the selected model. In particular, the frequency f is found by 
#' the grid search among 
#' \eqn{En = {j/101 | j=1,\dots,50 and j/101 \ge 1/n}}.
#' The MLE is equivalent to 
#' the least absolute estimation. 
#' The computation is completed by the Exterior Point Methods 
#' with theimported function \code{\link{rq.fit.br}} 
#' from the package \code{quantreg}.
#' 
#' @author Yuanhao Lai
#' 
#' @references  
#' Li, T. H. (2010). A nonlinear method for robust spectral analysis. 
#' Signal Processing, IEEE Transactions on, 58(5), 2466-2474.
#'                                                         
#' @keywords internal

GetFitHRegL1 <- function(y) {
  ##Set potential frequencies
  n <- length(y)
  t <- 1:n
  K <- 50
  aF <- 2.0*K+1
  lambda <- (ceiling(aF/n):K)/aF
  nlambda <- length(lambda)
  
  ##Set the design matrix
  X <- array(0,dim = c(n,3,nlambda))
  for(i in 1:nlambda){
    X1 <- cos(2*pi*lambda[i]*t)
    X2 <- sin(2*pi*lambda[i]*t) 
    X[,,i] <- cbind(1,X1,X2)  
  }
  
  ####################Computet the statistic######################
  SSE1 <- sum(abs(y-median(y)))   #SSE under NUll
  SSE2 <- numeric(nlambda)              #SSE under Alternative
  
  #Compute the statistic -2LLR under the laplace distribution 
  for(i in 1:nlambda){
    L1m2 <- rq.fit.br(x = X[,,i], y = y, tau = 0.5)
    SSE2[i] <- sum(abs(L1m2$residuals))
  }
  
  Index <- which.min(SSE2)
  LR <- n*log(SSE1/SSE2[Index])
  ans <- c(LR=LR,freq=lambda[Index]) #statistic and frequency
  ################################################################
  names(ans) <- c("L1","freq")
  return(ans)
}



