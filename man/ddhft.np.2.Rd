\name{ddhft.np.2}
\alias{ddhft.np.2}
\title{Data-Driven Haar-Fisz transformation}
\description{Forward Data-Driven Haar-Fisz transform}

\usage{
ddhft.np.2(data, Ccode=TRUE)
}
\arguments{
  \item{data}{A vector of size $2^J$ containing the data to variance
	stabilize and Gaussianize. }
  \item{Ccode}{If TRUE then fast C code is used, otherwise R code is used}
}
\details{Performs the data-driven Haar-Fisz transform on sequence data.
	This consists of (i) the Haar wavelet transform of sequence;
	(ii) estimation of mean-variance relationship between finest
	level smoothing and detail wavelet coefficients using
	isotonic regression (see \code{\link{isotone}}); (iii)
	divide wavelet detail coefficients by smooth ones subjected
	to the estimated mean-variance relationship; (iv) perform
	the inverse Haar wavelet transform of the modified coefficients.

	The aim is to variance stabilize and Gaussianize the sequence
	data which is only assumed to be positive and possess an
	underlying increasing mean-variance relationship.
}
\value{A list containing the following components
	\item{hft}{The data-driven Haar-Fisz transform of the input sequence}
	\item{mu}{The mu's obtained from the input sequence used for
		estimating the mean-variance relationship} 
	\item{sigma}{The estimated standard deviation as a function of
		the mean, the result of the isotonic regression fit of
		sigma2 on sigma}
	\item{sigma2}{The local multiscale standard deviations associated
		with each mean}
	\item{factors}{The numbers that divide the detail coefficients
		to standardize variance (obtained from the mean-variance
		estimation)}
}
\seealso{\code{\link{ddhft.np.inv}}}

\references{
Fisz, M. (1955), The limiting distribution of a function of two independent random variables 
and its statistical application, Colloquium Mathematicum, 3, 138-146.

Delouille, V., Fryzlewicz, P. and Nason, G.P. (2005), A data-driven Haar-Fisz transformation 
for multiscale variance stabilization. Technical Report, 05:06,
Statistics Group, Department of Mathematics, University of Bristol }
\examples{
#
# Generate example Poisson data set.
#
# Intensity function is steps from 1 to 32 in steps of 4 with each intensity
# lasting for 128 observations. Then sample Poisson with these intensities
#
v <- rpois(1024, lambda=rep(seq(from=1, to=32, by=4), rep(1024/8,8)))
#
# Let's take a look at this
#
\dontrun{ts.plot(v)}
#
# Ok. So mean of intensity clear increasing, but variance increasing too
#
# Now do data-driven Haar-Fisz
#
vhft <- ddhft.np.2(v)
#
# Now plot the variance stabilized series
#
\dontrun{ts.plot(vhft$hft)}
#
# The variance of the observations is much closer to 1. For example, let's
# look a the variance of the original series and the transformed one
#
# For the first intensity of 1
#
var(v[1:128])
#[1] 0.6628322
#
var(vhft$hft[1:128])
#[1] 1.025151
#
#
# And for second intensity of 5
#
#
var(v[129:256])
#[1] 4.389518
var(vhft$hft[129:256])
#[1] 1.312953
#
# So both transformed variances near to 1
#
# Now plot the estimated variance-mean relationship
#
\dontrun{plot(vhft$mu, vhft$sigma)}
\dontrun{lines(vhft$mu, sqrt(vhft$mu))}
#
# This is an approximately square root function (because you expect the
# sd of Poisson to be the square root of the mean).
#
}
\author{Piotr Fryzlewicz <p.fryzlewicz@imperial.ac.uk>
}
\keyword{manip}
