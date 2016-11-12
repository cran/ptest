#' The cumulative distribution function for Fisher'g
#' 
#' @description  
#' Compute the the cumulative distribution function for Fisher'g
#' from the selected method. See details.
#' 
#' @param g Fisher'g statistic.
#' @param n length of the series.
#' @param method method to compute the distribution function.
#'
#' @return the distribution function.
#'
#' @details method provides the following choices:
#' 
#' \code{exact}: Use the exact distribution.
#' 
#' \code{RSR}: Distribution fitted by the response surface regression method.
#' 
#' @author Yuanhao
#' 
#' @seealso \code{\link{ptestg}} 
#' 
#' @examples 
#' #Given the Fisher's g statistic, find the cumulative probability
#' pFisherg(g=0.3,n=10, method = "exact")
#' pFisherg(g=0.3,n=10, method = "RSR")
#' 
#' @references 
#' Fisher, R.A. (1929). Tests of significance in harmonic analysis. 
#' Proc. Roy. Soc. A, 125, 54-59. 
#' 
#' MacKinnon, James (2001) : Computing numerical distribution functions in 
#' econometrics, Queen's Economics Department Working Paper, No. 1037.
#'                                                         
#' @keywords internal
#' 
pFisherg <- function(g,n,
                     method = c("exact","RSR")){
  if(n<=0 | !n%%1==0 ) stop("n should be a positive integer")
  method <- match.arg(method)
  tablegEven <- NULL 
  tablegOdd <- NULL
  
  if(method=="exact"){
    #Define the exact cumulative density function (CDF) for the Fisher's g
    Fg <- function(g,n){
      m <- ifelse(n%%2==0,(n-2)/2,(n-1)/2)
      p <- floor(1/g)
      i <- 1:p
      1-sum(choose(m,i)*(-1)^(i-1) *(1-i*g)^(m-1))
    }
    CDF <- sapply(g,Fg,n=n)
    
  }else if(method=="RSR"){
    if(identical(n%%2,0)){
      data("tablegEven", package = "ptest", envir = environment())
      CDF <- 1-pvalueRSR(n, g, tablegEven)
    }else{
      data("tablegOdd", package = "ptest", envir = environment())
      CDF <- 1-pvalueRSR(n, g, tablegOdd)
    }
    
  }
  CDF
}