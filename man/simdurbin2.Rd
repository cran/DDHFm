\name{simdurbin2}
\alias{simdurbin2}
\title{Gene intensities simulator}
\description{
  Simulates gene intensities from the two components model of
  Durbin and Rocke.
}
\usage{
simdurbin2(mu, alpha, seta, seps)
}
\arguments{
  \item{mu}{the means where the gene intensities are simulated from}
  \item{alpha}{background mean of gene intensities}
  \item{seta}{standard deviation of high-level gene intensities}
  \item{seps}{standard deviation of low-level gene intensities}
}
\value{Like \code{\link{simdurbin}} this function simulates gene
 intensities from the Durbin and Rocke two component model (see
 Durbin and Rocke, 2001). Unlike \code{\link{simdurbin}} this function
 generates each intensity with a mean specified by an entry in
 \code{mu}. So, \code{length(mu)} intensities are generated with
 intensity \code{i} having mean \code{mu[i]}. 
}
\seealso{\code{\link{simdurbin}}}
\references{
Rocke, D.M. and Durbin, B.P. (2001), A model for measurement error for gene expression arrays, Journal of
Computational Biology, 8, 557-569
}
\examples{
#
# Generate data with mu equal to 80 and mu equal to 80000 with 10 replicates each
#
#
v <- simdurbin2( c(rep(80, 10), rep(80000,10)), alpha=24800, seta=0.227, seps=4800)
#
# Let's look at the mean of the first 10 and the last 10
#
mean(v[1:10])
#[1] 27361.95
#
mean(v[11:20])
#[1] 92455.66
}
\author{Guy Nason <g.p.nason@bris.ac.uk>}
\keyword{datagen}
