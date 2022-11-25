\name{hftrialdatagen}
\alias{hftrialdatagen}
\title{Gene intensities simulator and DDHFm tester}
\description{Simulates gene intensities and also applies DDHFm to them}

\usage{
hftrialdatagen(nreps = 4, nps = 128, plot.it = FALSE, uvp = 0.8)
}
\arguments{
  \item{nreps}{Number of replicates}
  \item{nps}{Number of genes}
  \item{plot.it}{Takes TRUE to activate the command of the respective plot and FALSE to deactivate it}
  \item{uvp}{a parameter for the denoising}
 }
\details{The code is well commented for further information.

  First, \code{\link{genesimulator}} is called to obtain a vector of
  mean gene intensities (for a number of genes and a number of replicates
  for each gene.

  Then \code{link{simdurbin2}} simulates a series of gene intensities
  using the (log-normal type) model as described in Durbin and Rocke  
  (2001,2002).

  Then for each gene the mean of replicates for that gene is computed.

  Optionally, if \code{plot.it} is \code{TRUE} then the mean 
  is plotted against its standard deviation (over
  replicates).

  Then the intensities are sorted according to increasing replicate mean.

  Optionally, if \code{plot.it} is \code{TRUE} then a plot of the
  intensities values as a vector (sorted according to increasing
  replicate mean) is plotted in black, and then the true mean plotted
  in colour 2 (on my screen this is red) and the computed replicate
  mean plotted in green.

  The DDHF transform of the sorted intensities is computed.

  Optionally, if \code{plot.it} is \code{TRUE} then a plot of the
  transformed means versus the transformed standard deviations is plotted.
  Followed by a time series plot of the transformed sorted intensities.
  These can be studied to see how well DDHF has done the transformation.

  Then two smoothing methods are applied the the DDHF transformed data.
  One method is translation invariant, Haar wavelet universal thresholding.
  The other method is the classical smoothing spline. If \code{plot.it} is
  \code{TRUE} then these smoothed estimates are plotted in different
  colours.

  Then the mean estimated intensity for each gene is  computed and this
  is returned as the first column of a two-column matrix (\code{ansm}).
  The second column is the true underlying mean. The object
  \code{hftssq} contains a measure of error between the estimated
  and true gene means. 
}

\value{
  \item{ansm}{Two column matrix containing the estimated gene intensities
	and the true ones}
  \item{hftssq}{Sum of squares between estimated means and true means}
  \item{yhf}{Simulated gene intensities}
  }
\examples{
#
# First run hftrialdatagen
#
\dontrun{v <- hftrialdatagen()}
#
# Now plot the Haar-Fisz transformed intensities.
#
\dontrun{ts.plot(v$yhf)}
#
# Now plot the denoised intensities 
#
# Note that above we have 128 genes and 4 replicates and so there are
# 4*128 = 512 intensities to plot.
#
# However, there are only 128 gene intensities, and estimates. So, for this
# plot we choose to plot the noisy intensities and then for each replicate
# group (which are colocated on the plot) plot the (necessarily constant)
# true and estimated intensities (ie we plot each true/estimated intensity
# 4 times, once for each replicate).
#
# First estimates...
#
\dontrun{lines(1:512, rep(v$ansm[,1], rep(4,128)), col=2)}
#
# Now plot the truth
#
\dontrun{lines(1:512, rep(v$ansm[,2], rep(4,128)), col=3)}
}
\author{Guy Nason <g.p.nason@bris.ac.uk>}
\keyword{datagen}
\keyword{manip}
