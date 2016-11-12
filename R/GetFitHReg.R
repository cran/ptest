#' Compute loglikelihood ratio test statistic
#' 
#' @description  
#' The -2 loglikelihood ratio test statistic, -2LLR, is computed 
#' for testing for periodicity.
#' 
#' @param z Vector containing the series.
#' @param t Vector of corresponding time points.
#'
#' @return Vector of length 2 containing -2 loglikelihood ratio 
#' test statistic and the estimated frequency.
#' 
#' @details 
#' To compute the likelihood ratio statistic, 
#' a harmonic regression model is fitted by 
#' the maximum likelihood estimations (MLE)
#' according to the selected model. In particular, the frequency f is found by 
#' the grid search among 
#' \eqn{En = {j/101 | j=1,\dots,50 and j/101 \ge 1/n}}.
#' The MLE is equivalent to 
#' the least square estimation. 
#' 
#' @author A.I. McLeod
#' 
#' @references  Islam, M.S. (2008).  Peridocity, 
#' Change Detection and Prediction in Microarrays. 
#' Ph.D. Thesis, The University of Western Ontario. 
#'                                                         
#' @keywords internal
GetFitHReg <-
function (y, t = 1:length(y)) 
{
    theta <- numeric(2)
    theta[1] <- length(y)
    ans <- .C("GetHReg", y = as.double(y), t = as.double(t), 
        theta = as.double(theta), PACKAGE = "ptest")
    ans <- ans$theta
    names(ans) <- c("LS","freq")
    ans
}
