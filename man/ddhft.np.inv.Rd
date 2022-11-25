\name{ddhft.np.inv}
\alias{ddhft.np.inv}
\title{Inverse Data-Driven Haar-Fisz transformation}
\description{
 Inverts the Data-Driven Haar-Fisz transform to obtain the raw
	(untransformed) data
}
\usage{
ddhft.np.inv(hft.obj)
}
\arguments{
	\item{hft.obj}{An object with the same structure as that returned
	by the \code{\link{ddhft.np.2}} function.}
}
\details{Merely performs the inverse of the DDHF transform see
	\code{\link{ddhft.np.2}}
}
\value{The inverted transform}

\seealso{\code{\link{ddhft.np.2}}}
 
\references{ 
Delouille, V., Fryzlewicz, P. and Nason, G.P. (2005), A data-driven Haar-Fisz transformation 
for multiscale variance stabilization. Technical Report 05:06,
Statistics Group, Department of Mathematics, University of Bristol.
}
\author{Piotr Fryzlewicz <p.fryzlewicz@imperial.ac.uk>
}
\keyword{manip}
