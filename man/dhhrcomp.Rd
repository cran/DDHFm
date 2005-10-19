\name{dhhrcomp}
\alias{dhhrcomp}
\title{Simulated genes, apply DDHFm then compute and return variance,
	skewness and kurtosis values}
\description{An example gene intensity simulator with given mean values mu0
 which then returns the variance, skew and kurtosis of the
 DDHF transformed intensities. The gene intensity
 simulation is performed using \code{\link{simdurbin2}}
}
\usage{
dhhrcomp(nsims = 1024, nmu = 4, mu0 = c(0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 65, 70, 80, 100, 200, 300, 500, 1000) * 1000)
}
\arguments{
  \item{nsims}{Number of replicates for each mu}
  \item{nmu}{Selects first \code{nmu} genes from \code{mu} for use}
  \item{mu0}{The possible gene intensity means}
  }
\details{This function simulates some gene intensities. A list of
	possible intensity means are supplied in \code{mu0}.
	The first \code{nmu} of these are selected. Then for each
	of the selected means \code{nsims} gene intensities are
	generated.

	The intensities are them subjected to DDHF transformation.

}
\value{
  \item{mu}{A vector of length \code{nmu} of the mean intensities considered}
  \item{v.hft}{A vector of \code{nmu} variances, one for each mean in
	\code{mu0}.
	The variance is the variance of the DDHFm transformed simulated
	intensity data.}
  \item{s.hft}{A vector of \code{nmu} skewnesses, one for each mean in
	\code{mu0}.
	The skewness is the skewness of the DDHFm transformed simulated
	intensity data.}
  \item{k.hft}{A vector of \code{nmu} kurtoses, one for each mean in
	\code{mu0}.
	The kurtosis is the kurtosis of the DDHFm transformed simulated
	intensity data.}
	
  }
\seealso{\code{\link{simdurbin2}}}
\references{
Durbin, B.P., Hardin, J.S., Hawkins, D.M. and Rocke, D.M. (2002), A variance-stabilizing transformation 
for gene expression microarray data, Bioinformatics, 18, S105-S110}
\author{Guy Nason <g.p.nason@bris.ac.uk.}
\keyword{datagen}
