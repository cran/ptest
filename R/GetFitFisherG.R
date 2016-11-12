#' Compute the Fisher's g test statistic
#' 
#' @description  
#' The Fisher's g test statistic is computed.
#' 
#' @param z Vector containing the series.
#'
#' @return Vector of length 2 containing the modified Fisher's G test statistic 
#' and the frequency for the maximum periodogram.
#' 
#' @author Yuanhao Lai
#' 
#' @references  Fisher, R.A. (1929). Tests of significance in harmonic analysis. 
#' Proc. Roy. Soc. A, 125, 54-59.
#'                                                         
#' @keywords internal
#' 
GetFitFisherG <- function(z) {
  n <- length(z)
  m <- ifelse(n%%2==0,(n-2)/2,(n-1)/2)
  Ip <- pgram(z)[,2]
  if( n%%2 ==0 )
    Ip<-Ip[-(m+1)]
  maxL <- which.max(Ip)
  g <- Ip[maxL]/sum(Ip)
  ans <- c(g=g,freq=maxL/n)
  ans
}








