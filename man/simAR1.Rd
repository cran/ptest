% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simAR1.R
\name{simAR1}
\alias{simAR1}
\title{Simulate AR(1) series}
\usage{
simAR1(n, phi = 0.3)
}
\arguments{
\item{n}{Length of series.}

\item{phi}{Autocorrelation parameter.}
}
\value{
Series of length n.
}
\description{
An AR(1) series with mean zero and variance 1 and
with autocorrelation paramater phi is simulated.
}
\details{
The model equation is:
z[t] = phi*z[t-1]+a[t],
where z[1] is N(0,1) and a[t] are NID(0, siga),
\eqn{siga=\sqrt(1/(1-phi^2))}.
}
\author{
A.I. McLeod
}
\references{
McLeod, A.I., Yu, Hao and Krougly, Z. (2007),  
Algorithms for Linear Time 
Series Analysis: With R Package, Journal of Statistical Software  23, 5 1-26.
}
\keyword{internal}

