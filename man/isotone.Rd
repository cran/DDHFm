\name{isotone}
\alias{isotone}
\title{Performs Isotone regression using "pool-adjacent-violators" algorithm}
\description{
  It is applied in order to estimate the h() function in the Data-Driven Haar-Fisz algorithm
}
\usage{
isotone(x, wt = rep(1, length(x)), increasing = FALSE, Ccode=TRUE)
}
\arguments{
  \item{x}{the vector that will be fitted with the regression}
  \item{wt}{a vector of weights}
  \item{increasing}{if TRUE the curve is set to be increasing, else FALSE}
  \item{Ccode}{if TRUE then faster C code is used rather than R}
}
\value{
  \item{isotone}{the regression results}
  }
\references{
Johnstone, I.M. and Silverman, B.W. (2005), EbayesThresh: R and S-Plus programs for empirical Bayes 
thresholding, Journal of Statistical Software, to appear
}
\author{Bernard Silverman, with C modifications by GPN}
\keyword{regression}
