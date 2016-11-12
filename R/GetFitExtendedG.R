#' Compute the extended Fisher's g test statistic
#' 
#' @description  
#' The extended Fisher's g te  statistic is computed for testing for periodicity.
#' 
#' @param y Vector containing the series.
#'
#' @return Vector of length 2 containing the modified Fisher's G test statistic 
#' and the frequency for the maximum periodogram.
#' 
#' @details Extend the Fisher's g test by enlarging the searching region 
#' of the frequency from the fourier frequencies to be 
#' \eqn{En = {j/101 | j=1,\dots,50 and j/101 \ge 1/n}}.
#' 
#' @author Yuanhao Lai
#' 
#' @references  Fisher, R.A. (1929). Tests of significance in harmonic analysis. 
#' Proc. Roy. Soc. A, 125, 54-59.
#'                                                         
#' @keywords internal

GetFitExtendedG <- function(y) {
  #Compute the extended Fisher's G test statistic 
  n <- length(y)

  # Mean removal
  y <- y-mean(y)
  
    
  t0 <- ceiling(101/n)
  Dpgram <- pgram(y, fr=101)
  Dpgram <- Dpgram[t0:50,]

  maxL <- which.max(Dpgram[,2])
  D <- sum(Dpgram[,2])
  gm <- Dpgram[maxL,2]/D
  ans <- c(gm,Dpgram[maxL,1]) #statistic and frequency
  names(ans) <- c("gm","freq")
  ans
}








