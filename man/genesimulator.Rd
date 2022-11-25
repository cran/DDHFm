\name{genesimulator}
\alias{genesimulator}
\title{Gene means simulator}
\description{Simulates means of the gene intensities for "nps" genes, each
	replicated "nreps" times. The mean is drawn from a gamma distribution
	with shape parameter "shape" and scale parameter "scale" }
\usage{
genesimulator(nreps = 3, nps = 100, shape = 4, scale = 100)
}
\arguments{
  \item{nreps}{Number of replicates}
  \item{nps}{Number of genes}
  \item{shape}{Shape parameter}
  \item{scale}{Scale parameter}
}
\details{For many problems a set of reasonable gene mean intensities is
	useful for testing algorithms that later draw actual gene intensities
	with distributions that possess a gene mean intensity.
}
\value{
  A matrix containing nreps*nps rows and 3
	columns. The first col contains all the gene mean
	intensities. The second and third col contain the gene
	replicate number and gene number respectively.
}
\author{Guy Nason <g.p.nason@bris.ac.uk>}
\keyword{datagen}
