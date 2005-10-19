\name{cDNAdata}
\alias{cDNAdata}
\title{Samples from the cDNA data vector}
\description{Arranges a vector of intensities into a form amenable to analysis
	by the DDHF and also can restrict the number of genes analysed and
	also obtain a random sample}
}
\usage{
 cDNAdata(data.vect,cdnalength,datasize,ng,nrep)
}
\arguments{
  \item{data.vect}{The data vector of intensities}
  \item{cdnalength}{Only considers the first \code{cdnalength} genes for analysis}
  \item{datasize}{Needs to be a power of two. The number of genes that get
	randomly sampled from the first \code{cdnalengthq} and then
	subsequently transformed.}
  \item{ng}{The total number of genes described by \code{data.vect}} 
  \item{nrep}{The number of replicates (should be a power of two) contained
	in \code{data.vect}}
}
\details{
   The $J = ng \times nrep$ vector data.vect should contain first the
	intensities of the first replicate of the \code{ng} genes',
	then the second replicate of all of the \code{ng} genes in the
	same order, and so on.

   This function first puts the one dimensional \code{data.vect} into
   a matrix form with \code{ng} rows and \code{nrep} columns (so the
   row number indices the gene and the column number the replicate).

   Then the first \code{cdnalength} rows are extracted and kept, the remaining
   rows are discarded.

   Then \code{datasize} rows at random are extracted and kept and
	the remaining rows are discarded.
} 
\value{
  \item{cDNAdata}{The data vector in the proper format to perform Data-Driven Haar-Fisz algorithm}
  }
\author{Efthimios Motakis <e.motakis@bris.ac.uk>}
\keyword{datagen}
