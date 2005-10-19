\name{dhhrss}
\alias{dhhrss}
\title{Tabulates variance, skewness and kurtosis coefficients from
	the output of dhhrcomp}
\description{
  Tabulates the variance, skewness and kurtosis coefficients from the
	output of \code{\link{dhhrcomp}}
}
\usage{
dhhrss(dhhrobj)
}
\arguments{
  \item{dhhrobj}{The results to be tabulated. The output of a call
	to \code{\link{dhhrcomp}} }
}
\details{This routine merely runs summary statistics calculations
	on the output results from \code{\link{dhhrcomp}}
}
\value{None, prints a table.
  }
\author{Guy Nason <g.p.nason@bris.ac.uk.}
\keyword{manip}
