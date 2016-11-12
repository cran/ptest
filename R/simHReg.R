#' Simulate harmonic regression models
#' 
#' @description  
#' Simulates a harmonic regression.  Possible types of models are 
#' normal, t(5), Laplace, cubic and AR1.
#' 
#' @param n Length of series.
#' @param f Frequency.
#' @param A Cosine amplitude.
#' @param B Sine amplitude.
#' @param model The model used for generating the error term. See details.
#' @param phi Only used if AR1 error distribution is selected.
#' @param sig The standard error of the series.
#'
#' @return Vector of length n, simulated harmonic series.
#' 
#' @details Generate a harmonic series y with length n, 
#' where \eqn{y_t = A*cos(2*pi*f*t)+B*sin(2*pi*f*t)+sig*e_t,\ t=1,...,n,}
#' and e comes from one of the following specified distributions 
#' with mean 0 and standard error 1:
#' 
#' \code{Gaussian}: A standard normal distribution (i.i.d.).
#' 
#' \code{t5}: A t distribution with 5 degrees of freedom 
#' (i.i.d., standardized to mean 0 and variance 1).
#' 
#' \code{Laplace}: A Laplace (double exponential) distribution
#' (i.i.d., standardized to mean 0 and variance 1).
#' 
#' \code{cubic}: A standard normal distribution for e, 
#' but \eqn{y=y^3} this time.
#' 
#' \code{AR1}: An AR(1) series with autocorrelation paramater phi
#' (standardized to mean 0 and variance 1).
#' 
#' @author A.I. McLeod and Yuanhao Lai
#' 
#' @references McLeod, A.I., Yu, Hao and Krougly, Z. (2007),  
#' Algorithms for Linear Time 
#' Series Analysis: With R Package, Journal of Statistical Software  23, 5 1-26.
#' 
#' @seealso \code{\link{fitHReg}}, \code{\link{ptestReg}}
#' 
#' @examples 
#' #Simulate the harmonic regression model with standard Gaussian error terms
#' z <- simHReg(10, f=2/10, 1, 2, model="Gaussian",sig=1) #Fourier Frequency
#' plot(1:10,z,type="b")
#' 
#' #Simulate the AR(1) errors
#' z <- simHReg(10, f=0/10, 0,0, model="AR1",phi=0.2,sig=1) 
#' acf(z)
#' 
#' @keywords ts
simHReg <- function(n, f, A, B, 
                    model = c("Gaussian", "t5", "Laplace", "cubic","AR1"), 
                    phi = 0, sig = 1){
  model <- match.arg(model)
  if(phi>=1 | phi<=-1) stop("phi should be in the range (-1,1) for stationarity!")
  if(sig<=0) stop("sig should be positive!")
  
  t<-1:n
  x<-A*cos(2*pi*f*t)+B*sin(2*pi*f*t)
  
  if(model=="cubic"){
    e <- rnorm(n)
    e <- sig*e
    z<-(x+e)^3
  }else{
    e<-switch(model,
              Gaussian=rnorm(n), 
              t5=rt(n,df = 5)/sqrt(5/3), 
              Laplace=simLaplace(n),
              AR1=simAR1(n,phi) )
    e <- sig*e
    z<-x+e
  }
  
  return(z)
}
