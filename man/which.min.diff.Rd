\name{which.min.diff}
\alias{which.min.diff}
\title{Find index where two vectors are closest}
\description{Finds the index point where two index vectors are closest
 in value}
\usage{
which.min.diff(a, vect)
}
\arguments{
  \item{a}{input vector}
  \item{vect}{input vector}
  }
\value{The index where \code{abs(a-vect)} is smallest, ie where \code{a}
	and \code{vect} are smallest.
  }
\examples{
#
# Make up two vectors
#
a <- c(1,2,3)
vect <- c(3,2,1)
#
# Now see on which index are the two closest in value 
#
which.min.diff(a,vect)
#[1] 2
#
# ie its the second index where both vectors are actually 2
}
\author{Piotr Fryzlewicz <p.fryzlewicz@imperial.ac.uk>}
\keyword{methods}
