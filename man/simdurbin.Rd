\name{simdurbin}
\alias{simdurbin}
\title{Gene intensities simulator}
\description{
  Simulates gene intensities from the two components model as in
	the papers of Durbin and Rocke.
}
\usage{
simdurbin(n, alpha, mu, seta, seps)
}
\arguments{
  \item{n}{number of intensities to be simulated}
  \item{alpha}{background mean of gene intensities}
  \item{mu}{the means where the gene intensities are simulated from}
  \item{seta}{standard deviation of high-level gene intensities}
  \item{seps}{standard deviation of low-level gene intensities}
}
\details{This function generates \code{n} intensities from the
  Durbin and Rocke (2001) gene intensity model with parameters
  alpha, seta and seps. 
}
\value{A vector of the appropriate intensities.
}
\seealso{\code{\link{simdurbin2}}}
\references{
Rocke, D.M. and Durbin, B.P. (2001), A model for measurement error for gene expression arrays, Journal of 
Computational Biology, 8, 557-569
}
\author{Guy Nason <g.p.nason@bris.ac.uk>}
\keyword{datagen}
