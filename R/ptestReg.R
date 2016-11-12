#' Test short time series for periodicity with maximum likelihood ratio tests
#' 
#' @description  
#' This function is used to test the existence of the periodicity for 
#' a short time series (length<=100). Likelihood ratio tests under 
#' the Gaussian or the Laplace assumptions are provided with 
#' the response surface method implemented for efficiently obtaining 
#' accurate p-values.
#' 
#' @param z A series or a matrix containg series as columns
#' @param method The statistical test to be used. See details for more information.
#' @param multiple Indicating whether z contains multiple series.
#' 
#' @return Object of class "Htest" produced.
#' 
#' An object of class "Htest" is a list containing the following components:
#' 
#' \item{obsStat}{Vector containing the observed test statistics.}
#' \item{pvalue}{Vector containing the p-values of the selected tests.}
#' \item{freq}{Vector containing the estimated frequencies.}
#' 
#' @details The null hypothesis is set as no peridicities, H0: f=0. 
#' Discriptions of different test statistics (methods) are as follow:
#' 
#' \code{LS}: The -2 loglikelihood ratio test statistic based on 
#' the likelihood ratio test with normal noises, 
#' where the p-values are efficiently computed by 
#' the response surface method.
#' 
#' \code{L1}: The -2 loglikelihood ratio test statistic based on 
#' the likelihood ratio test with Laplace noises, 
#' where the p-values are efficiently computed by 
#' the response surface method.
#' 
#' @author Yuanhao Lai and A.I. McLeod
#' 
#' @references
#' 
#' Islam, M.S. (2008). Peridocity, Change Detection and Prediction in Microarrays. 
#' Ph.D. Thesis, The University of Western Ontario. 
#' 
#' Li, T. H. (2010). A nonlinear method for robust spectral analysis. 
#' Signal Processing, IEEE Transactions on, 58(5), 2466-2474.
#' 
#' MacKinnon, James (2001) : 
#' Computing numerical distribution functions 
#' in econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @seealso \code{\link{fitHReg}, \link{ptestg}}
#' 
#' @examples 
#' # Simulate the harmonic regression model with standard Gaussian error terms
#' set.seed(193)
#' # Non-Fourier frequency
#' z <- simHReg(n = 14, f=2/10, A = 2, B = 1, model="Gaussian",sig=1) 
#' ptestReg(z,method = "LS") #Normal likelihood ratio test
#' ptestReg(z,method = "L1") #Laplace likelihood ratio test  
#' fitHReg(z, algorithm="exact") #the nls fitted result 
#'     
#'                                            
#' # Performe tests on the alpha factor experiment
#' data(alpha)
#' ## Eliminate genes with missing observations
#' alpha.nonNA <- alpha[complete.cases(alpha),]
#' ## Using the multiple option to do the test for all the genes
#' ## Transpose the data set so that each column stands for a gene
#' alpha.nonNA <- t(alpha.nonNA)
#' result <- ptestReg(alpha.nonNA, method = "LS",multiple=TRUE) 
#' str(result)       
#'
#'
#' # The movtivating example: gene ORF06806 in Cc
#' data(Cc)
#' x <- Cc[which(rownames(Cc)=="ORF06806"),]
#' plot(1:length(x),x,type="b", main="ORF06806",
#'      xlab="time",ylab="Gene expression")
#' ptestg(x,method="Fisher") #Fail to detect the periodicity
#' ptestReg(x,method="LS") #The periodicity is significantly not zero
#' ptestReg(x,method="L1") #The periodicity is significantly not zero
#'
#' @keywords ts
ptestReg <- function(z, method=c("LS", "L1"),
                  multiple=FALSE){
  method <- match.arg(method)
  
  if(multiple == FALSE){
    n <- length(z)
    m <- 1
    z <- as.matrix(z,ncol=1,drop=FALSE)
  }else{
    nm <- dim(z)
    if(is.null(nm)) stop("There is no multiple series")
    n <- nm[1]
    m <- nm[2]
    if(n<=1 | m<=1) stop("There is no multiple series")
  }

  if(n>100){ stop("length is larger than 100! Fisher's g test is suggested.") }
  obsStat <- numeric(m)
  freq <- numeric(m)
  
  if(method=="LS"){ #the likelihood ratio test with normal noises
    tableRegLs <- NULL
    for (i in 1:m) {
      obsStatFreq <- GetFitHReg(z[,i])
      obsStat[i] <- obsStatFreq[1]
      freq[i] <- obsStatFreq[2]
    }
    data("tableRegLs", package = "ptest", envir = environment())
    pvalue <- pvalueRSR(n, obsStat, tableRegLs)
    
  }else if(method=="L1"){ #the likelihood ratio test with Laplace noises
    tableRegL1Even <- NULL
    tableRegL1Odd <- NULL
    for (i in 1:m) {
      obsStatFreq <- GetFitHRegL1(z[,i])
      obsStat[i] <- obsStatFreq[1]
      freq[i] <- obsStatFreq[2]
    }
    
    if(identical(n%%2,0)){
      data("tableRegL1Even", package = "ptest", envir = environment())
      pvalue <- pvalueRSR(n,  obsStat, tableRegL1Even)
    }else{
      data("tableRegL1Odd", package = "ptest", envir = environment())
      pvalue <- pvalueRSR(n, obsStat, tableRegL1Odd)
    }
  }
  
  ans <- list(obsStat=obsStat, pvalue=pvalue, freq=freq)
  class(ans)<-"Htest"
  return(ans)
}
