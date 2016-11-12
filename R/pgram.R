#' Periodogram computation
#' 
#' @description  
#' The periodogram is computed.
#' 
#' @param z time series vector of length n, say.
#' @param fr use "default" for usual Fourier frequencies, {1/n, ..., floor(n/2)/n}. 
#' Set fr = N, to evaluate the periodogram at the Fourier frequencies 
#' corresponding to a time series of length N. Finally set fr to any
#' desired set of frequencies. Note frequencies are in cycles per unit
#' time sometimes called temoral frequency to distinguish from angular
#' frequency. Both are widely used in time series.
#' @param method either periodogram or regression
#' 
#' @return Periodogram evaluated at the Fourier frequencies or R-square.
#' 
#' @details Uses FFT.
#' So if the length of z is a highly composite number, the computation
#' is very efficient. Otherwise the usual DFT is used.
#' 
#' @author A.I. McLeod and Yuanhao Lai
#' 
#' @examples 
#' z<-sunspot.year
#' n<-length(z) 
#' I<-pgram(z)
#' f<-I[,1]
#' I <- I[,2]
#' plot(f, I, xlab="f", ylab="f", type="l") 
#' title(main="Periodogram for Annual Sunpots, 1700-1988") 
#' #
#' z<-c(0.42, 0.89, 1.44, 1.98, 2.21, 2.04, 0.82, 0.62, 0.56, 0.8, 1.33)
#' fr <- (1:50)/101
#' pgram(z)
#' pgram(z, fr=101)
#' pgram(z, fr=fr)
#' pgram(z, method="regression")
#' pgram(z, method="regression", fr=101)
#' pgram(z, method="regression", fr=fr)
#'   
#'                                                         
#' @keywords ts

pgram <-
  function(z, fr="default", method=c("periodogram", "regression"))
  {
    methd = match.arg(method)
    n <- length(z)
    if (methd=="periodogram") {
      if (identical(fr, "default")) {
        ans <- cbind((1:floor(n/2))/n, (Mod(fft(z))^2/n)[2:(n %/% 2 + 1)])
      } else {
        stopifnot(is.numeric(fr))
        if (length(fr) == 1) {
          madj <- fr-n
          x <- c(z, rep(0, madj))
          ans <- cbind((1:floor(fr/2))/fr, (Mod(fft(x))^2 /n)[2:(fr%/%2+1)])
        } else {
          ans <- cbind(fr, Re(
            sapply(fr, function(f0) {
              y <- sum(z*exp(1i * 2*pi*f0 * seq(0,n-1)))
              y*Conj(y)}
            )/length(z))
          )
        }
      }
      colnames(ans) <- c("frequency", "periodogram")
      return(ans)
    } else {#regression - F ratio
      SSTot <- sum((z-mean(z))^2)
      sw <- c("a","b","c")[c(identical(fr, "default"), 
                             is.numeric(fr)&&length(fr)==1, 
                             is.numeric(fr)&&length(fr)>1)]
      freq <- switch(sw,
                     a = (1:floor(n/2))/n,
                     b = (1:floor(fr/2))/fr,
                     c = fr
      )
      SSErr <- NULL
      t <- 1:n
      for (k in 1:length(freq)) {
        A <- cbind(rep(1,n), cos(2*pi*freq[k]*t), sin(2*pi*freq[k]*t))
        res <- qr.resid(qr(A), z)
        SSErr[k] <- sum(res^2)
      }
      ans <- cbind(freq, 1-SSErr/SSTot)
      colnames(ans) <- c("frequency", "RSq")
      return(ans)
    }
  }

