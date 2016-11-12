#' Fits Three Parameter Harmonic Regression
#' 
#' @description  
#' Estimates A, B and f in the harmonic regression,
#' y(t)=mu+A*cos(2*pi*f*t)+B*sin(2*pi*f*t)+e(t). 
#' The default algorithm is enumerative but an exact non-linear LS
#' option is also provided.
#'
#' @param y series.
#' @param t Time points.
#' @param algorithm method for the optimization
#' 
#' @return Object of class "HReg" produced.
#' 
#' @details Program is interfaced to C for efficient computation.
#' 
#' @author A.I. McLeod and Yuanhao Lai
#' 
#' @examples 
#' set.seed(193)
#' z <- simHReg(10,f=2.5/10,1,1)
#' ans <- fitHReg(z)
#' ans$freq #optimal frequency = 0.2376238 
#' #
#' #ORF06806 in Cc dataset. 
#' z<-c(0.42, 0.89, 1.44, 1.98, 2.21, 2.04, 0.82, 0.62, 0.56, 0.8, 1.33)
#' ans2 <- fitHReg(z, algorithm="exact")
#' sum(resid(ans2)^2) #0.2037463
#' ans1 <- fitHReg(z)
#' sum(resid(ans1)^2) #0.242072
#' #compare with nls()
#' t <- 1:length(z)
#' ans <- nls(z ~ mu+alpha*cos(2*pi*lambda*t+phi), 
#'               start=list(mu=1, alpha=1, lambda=0.1, phi=0.0))
#' coefficients(ans)
#' sum(resid(ans)^2) #0.2037
#' 
#'                                                         
#' @keywords ts

fitHReg <-
  function (y, t = 1:length(y), algorithm=c("enumerative", "exact")) 
  {
    alg <- match.arg(algorithm)
    n <- length(y)
    
    theta <- numeric(2)
    theta[1] <- length(y)
    ansH <- .C("GetHReg", y = as.double(y), t = as.double(t), 
               theta = as.double(theta), PACKAGE = "ptest")
 #   LR <- (ansH$theta)[1]
    fopt <- (ansH$theta)[2]
    
    if(alg=="exact"){
      #for ORF06806 example this converges to wrong answer!
      # fopt <- optimize(f=function(fr) sum(qr.resid(
      #        qr(cbind(1, cos(2*pi*fr*t), sin(2*pi*fr*t))),y)^2),
      #        interval=c(1/n, 0.5))$minimum
      # alternate formulation using nls()
      ansNls <- nls(y ~ mu+alpha*cos(2*pi*lambda*t+phi),
                    control=list(maxiter=100,minFactor=1/5000),
                    start=list(mu=1, alpha=1, lambda=fopt, phi=0.0))
      fopt <- coef(ansNls)["lambda"]
      # alternative using optim() with initial value
      # fopt <- optim(par = fopt, fn=function(fr){sum(qr.resid(
      #               qr(cbind(1, cos(2*pi*fr*t), sin(2*pi*fr*t))),y)^2)},
      #               method = "L-BFGS-B",lower = 1/n, upper = 0.5
      #              )$par
      # theta <- numeric(2)
      # theta[1] <- length(y)
      # ansH <- .C("GetHRegExact", y = as.double(y), t = as.double(t), 
      #            theta = as.double(theta), PACKAGE = "ptest")
      # #   LR <- (ansH$theta)[1]
      # fopt <- (ansH$theta)[2]
    }

    x1<-cos(2*pi*fopt*t)
    x2<-sin(2*pi*fopt*t)
    ansLM<-lm(y~x1+x2)
    co<-coef(ansLM)
    names(co)<-c("mu","A","B")
    re <- residuals(ansLM)
    a<-summary(ansLM)
    v<-a$cov.unscaled
    Rsq<-a$r.squared
    fstat<- a$fstatistic
    sig<-a$sigma
    SSE1 <- var(y)*(n-1)
    SSE2 <- sum(re^2)
    LR <- n*log(SSE1/SSE2)
    ans<-list(coefficients=co, residuals=re, Rsq=Rsq, 
              fstatistic=fstat, sigma=sig, freq=fopt, LRStat=LR)
    class(ans)<-"HReg"
    ans
  }

