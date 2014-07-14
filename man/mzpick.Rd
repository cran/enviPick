\name{mzpick}
\alias{mzpick}

\title{EIC peak picking}

\description{
Peak-picking within individual EIC cluster formed by \code{\link[enviPick]{mzclust}} without assuming a certain peak shape. 
Includes a baseline subtraction step.
}

\usage{
mzpick(MSlist, minpeak = 4, drtsmall = 20, drtfill = 10,  drttotal = 200, recurs = 4, 
weight = 2, SB = 3, SN=2, minint = 1E4, maxint = 1e+07, ended = 2, progbar = FALSE, 
from = FALSE, to = FALSE)}

\arguments{
  \item{MSlist}{An MSlist returned by \code{\link[enviPick]{mzclust}}}
  \item{minpeak}{Peak definition - minimum number of measurements required within the RT window of \code{drtsmall}}
  \item{drtsmall}{Peak definition - RT window of a peak; cp. \code{minpeak} argument}
  \item{drtfill}{Maximum RT gap length to be filled, cp. details}
  \item{drttotal}{Peak definition - Maximum RT length of a single peak}
  \item{recurs}{Maximum number of peaks within one EIC, cp. details}
  \item{weight}{Weight for assigning measurements to a peak, cp. details}
  \item{SB}{Peak definition - signal-to-base ratio}
  \item{SN}{Peak definition - signal-to-noise ratio}
  \item{minint}{Peak definition - minimum intensity of a peak}
  \item{maxint}{Peaks above this intensity are always retained, regardless of other checks}
  \item{ended}{Within the peak detection recursion set by argument \code{recurs}, how often can a peak detection fail to end the recursion? }
  \item{progbar}{For debugging, ignore}
  \item{from}{For debugging, ignore}
  \item{to}{For debugging, ignore}
}


\details{
In a first step, RT gaps between measurements in an EIC not larger than \code{drtfill} are filled by linear interpolation.
Subsequently, peaks are assigned over a number of \code{recurs} recursions not interrupted by more than \code{ended} times of failed peak detections. 

At each recursion, the most intense EIC measurement not yet assigned to a peak is selected as peak apex and neighbouring unassigned measurements at lower and higher RT are evaluated 
for forming the peak. To this end, increases (lower RT) and decreases (higher RT) in intensity of consecutive measurements over a maximum RT width of \code{drtdens} are summed and 
penalized by a factor of \code{weight} for intensity reversions. The measurements with optimum values are then selected to define the start and end
measurement of the peak.

Thereupon, the candidate peak is checked to  
(a) have at least \code{minpeaks} within a RT window of \code{drtsmall}, 
(b) be larger than the minimum peak intensity \code{minint} and
(c) have a minimum \code{SB} ratio (the ratio between the most intense measurement and the mimimum intensity of the first or last peak measurement). 
Candidate peaks failing in any of the aspects (a) to (c) are discarded (adding to \code{ended}), unless they are higher in intensity than \code{maxint}.

Next, all measurements assigned to peaks are removed from the EIC and the resulting gaps linearly interpolated and smoothed by a moving window average to form a baseline.
The latter is then subtracted from the assigned peaks. 

In a last step, peaks are checked for their signal-to-noise \code{SN} ratio in relation to the baseline measurements (if present).
Herein, \code{SN} is defined as the ratio between the most intense (baseline-corrected) peak measurement and the median of the difference 
between the non-peak measurements (if any) and the baseline. 
}



\value{
	Returns the argument MSlist, with entries made:
	
	\item{Parameters}{MSlist[[2]]: saves the parameter settings.}
	\item{Scans}{MSlist[[4]]: matrix with raw measurements and tags for picked peaks within EICs within partitions.}
	\item{Peak_index}{MSlist[[7]]: Index assigning picked peaks to sections in the raw measurment of MSlist[[4]]. Required for fast access.}
	\item{Peaklist}{MSlist[[8]]: matrix with picked peak characteristics, with columns: 
		m/z (mean m/z of peak measurements), 
		var_m/z (m/z variation of peak measurements),
		max_int (base-line corrected maximum intensity),
		sum_int (sum of all base-line corrected peak measurement intensities),
		RT (retention time at maximum intensity),
		minRT (start peak RT),
		maxRT (end peak RT),
		peak# (peak ID number),
		EIC# (EIC ID number) and
		Score (not yet implemented).
		}
}

\section{imbecile}{
\code{ended} must be smaller than \code{recurs}.
\code{minpeak} and \code{drtsmall} should be congruent in both and \code{\link[enviPick]{mzclust}} and \code{\link[enviPick]{mzpick}}.
}

\author{Martin Loos}

\seealso{
	\code{\link[enviPick]{plotMSlist}}
	\code{\link[enviPick]{writePeaklist}}
}

