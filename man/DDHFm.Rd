\name{DDHFm}
\alias{DDHFm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data-driven Haar-Fisz for microarrays}
\description{Takes a matrix containing the resuls of a (e.g.) microarray
	experiment (rows are different genes, columns are replicates)
	and performs a multiscale variance stabilization and Gaussianization
	transform. Note, you have to have replicates to make this work.
}
\usage{
DDHFm(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix containing n rows (genes) and k columns (replicates)}
}
\details{The data-drive Haar-Fisz transform is a multiscale variance
	stabilization and Gaussianization transform.
	The algorithm will run on a one-column matrix but the results
	will NOT be good. For one column of data use some other
	variance stabilizer (in which case why are you using one column,
	eh?)

	However, if you have replicates then you can use this
	powerful multiscale algorithm which will stabilize variance
	and also draw the distribution of the intensities towards
	Gaussian extremely effectively (so skew and kurtosis will look
	good too).

	This function forms a single vector of intensities by the
	following procedure. First, it computes the mean over replicates
	for each gene (this is why you need replicates, really).
	Then it orders the rows of the matrix into increasing mean
	order. Then it converts that matrix into a vector (which is therefore
	a sequence of intensities which are grouped according to gene
	(according to increasing mean intensity within each group)).
	Then this vector is padded up to the next power-of-two in length
	with zeroes. Then the padded vector is subjected to the
	\code{\link{ddhft.np.2}} function which actually does the
	data-driven Haar-Fisz transform. The results of this are
	unpadded and then rearranged into a matrix with an identical
	structure to the input matrix. This final matrix contains
	stabilized and Gaussianized (is this a word?) intensities.

	It is important that you have replicates (ie, the number of
	columns of the matrix has to be greater than 1). Also, it
	is scientifically important that the replicates be just that
	(ie, although the rows might refer to the same genes the columns
	must be results from repeats experiments held under identical
	conditions). For example, you shouldn't have results for the
	gene on one subject and then the same gene on another subject
	unless I suppose you are considering them to be replicates in some way.
}
\value{A single matrix with the same dimensions as the input matrix.
	Each entry in the resultant matrix corresponds to a variance
	stabilized and Gaussianized version of the entry in the input
	matrix.
}
\references{Motakis, E.S., Nason, G.P., Fryzlewicz, P. and Rutter, G.A.
	(2005)
	Variance stabilization and normalization for one-color
	microarray data using a data-driven multiscale approach.
	\emph{Technical Report}, 05:16, Statistics Group,
	Department of Mathematics, University of Bristol, UK}
\author{Guy Nason}

\seealso{\code{\link{ddhft.np.2}}}
\examples{
#
# First, let's make up some "true" intensities (these are gamma)
#
TrueInt <- genesimulator(nreps=5, nps=100)
#
# Now let us simulate data from the Durbin and Rocke (2001) model conditioning
# on the TrueInt
#
MAInts <- simdurbin2(TrueInt,alpha=24800, seta=0.227, seps=4800)
#
# Now put these intensities into the correct format for DDHFm. That is
# replicates across columns
#
MAm <- matrix(MAInts, ncol=5, byrow=TRUE)
#
# Apply DDHFm transform
#
MAmDDHFm <- DDHFm(MAm)
#
# Now look at the variances of each
#
diag(var(MAm))
#
#[1] 22443437 25065729 23430160 22470058 21023183
#
diag(var(MAm))/ mean(diag(var(MAm)))
#
#[1] 0.9806403 1.0952183 1.0237540 0.9818035 0.9185839
#
# Now for the DDHFm version
#
diag(var(MAmDDHFm))
#
# [1] 10938.71 11369.41 11169.50 11144.74 10959.94
#
diag(var(MAmDDHFm))/ mean(diag(var(MAmDDHFm)))
#
# [1] 0.9840099 1.0227549 1.0047711 1.0025441 0.9859200
#
# So stabilization is better with DDHFm 
#
}
\keyword{manip}% at least one, from doc/KEYWORDS
