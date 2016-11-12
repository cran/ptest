#' Fisher's g test
#' 
#' @description  
#' Calculates p-value using Fisher's g-test for a single frequency.
#' 
#' @param z Input series, no missing values.
#' 
#' @return Vector of length 3 containing the g test statistic, the p-value, 
#' and the frequency for the maximum periodogram.
#' 
#' @details Program is interfaced to C for efficient computation.
#' 
#' @author A.I. McLeod
#' 
#' @references Fisher, R.A. (1929). Tests of significance in harmonic analysis. 
#' Proc. Roy. Soc. A, 125, 54-59. 
#'                                                         
#' @keywords internal

FisherGTest <-
function(z){
    n <- length(z)
    m <- ifelse(n%%2==0,(n-2)/2,(n-1)/2)
    Ip <- pgram(z)[,2]
    if( n%%2 ==0 )
        Ip<-Ip[-(m+1)]
    maxL <- which.max(Ip)
    g <- Ip[maxL]/sum(Ip)
    p <- floor(1/g)
    i <- 1:p
    pvalue <- sum(choose(m,i)*(-1)^(i-1) *(1-i*g)^(m-1))
    ans <- c(gstat=g,pvalue=pvalue,freq=maxL/n)
    ans
 }

