% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net.scores.R
\name{net.scores}
\alias{net.scores}
\title{Network Scores}
\usage{
net.scores(
  data,
  A,
  wc,
  loading.method = c("original", "revised"),
  rotation = NULL,
  scores = c("Anderson", "Bartlett", "components", "Harman", "network", "tenBerge",
    "Thurstone"),
  loading.structure = c("simple", "full"),
  impute = c("mean", "median", "none"),
  ...
)
}
\arguments{
\item{data}{Matrix or data frame.
Should consist only of variables to be used in the analysis}

\item{A}{Network matrix, data frame, or \code{\link[EGAnet]{EGA}} object}

\item{wc}{Numeric or character vector (length = \code{ncol(A)}).
A vector of community assignments.
If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
then \code{wc} is automatically detected}

\item{loading.method}{Character (length = 1).
Sets network loading calculation based on implementation
described in \code{"original"} (Christensen & Golino, 2021) or
the \code{"revised"} (Christensen et al., 2024) implementation.
Defaults to \code{"revised"}}

\item{rotation}{Character.
A rotation to use to obtain a simpler structure.
For a list of rotations, see \code{\link[GPArotation]{rotations}} for options.
Defaults to \code{NULL} or no rotation.
By setting a rotation, \code{scores} estimation will be
based on the rotated loadings rather than unrotated loadings}

\item{scores}{Character (length = 1).
How should scores be estimated?
Defaults to \code{"network"} for network scores.
Set to other scoring methods which will be computed using
\code{\link[psych]{factor.scores}} (see link for arguments
and explanations for other methods)}

\item{loading.structure}{Character (length = 1).
Whether simple structure or the saturated loading matrix
should be used when computing scores.
Defaults to \code{"simple"}

\code{"simple"} structure more closely mirrors sum scores and CFA;
\code{"full"} structure more closely mirrors EFA

Simple structure is the more "conservative" (established) approach
and is therefore the default. Treat \code{"full"} as experimental
as proper vetting and validation has not been established}

\item{impute}{Character (length = 1).
If there are any missing data, then imputation can be implemented.
Available options:

\itemize{

\item \code{"none"} --- Default. No imputation is performed

\item \code{"mean"} --- The mean value of each variable is used to replace missing data
for that variable

\item \code{"median"} --- The median value of each variable is used to replace missing data
for that variable

}}

\item{...}{Additional arguments to be passed on to
\code{\link[EGAnet]{net.loads}} and
\code{\link[psych]{factor.scores}}}
}
\value{
Returns a list containing:

\item{scores}{A list containing the standardized (\code{std.scores})
rotated (\code{rot.scores}) scores. If \code{rotation = NULL}, then
\code{rot.scores} will be \code{NULL}}

\item{loadings}{Output from \code{\link[EGAnet]{net.loads}}}
}
\description{
This function computes network scores computed based on
each node's \code{strength} within each community in the network
(see \code{\link[EGAnet]{net.loads}}). These values are used as "network loadings"
for the weights of each variable.

Network scores are computed as a formative composite rather than a reflective factor.
This composite representation is consistent with no latent factors that psychometric
network theory proposes.

Scores can be computed as a "simple" structure, which is equivalent to a weighted
sum scores or as a "full" structure, which is equivalent to an EFA approach.
Conservatively, the "simple" structure approach is recommended until further
validation
}
\examples{
# Load data
wmt <- wmt2[,7:24]

# Estimate EGA
ega.wmt <- EGA(
  data = wmt,
  plot.EGA = FALSE # No plot for CRAN checks
)

# Network scores
net.scores(data = wmt, A = ega.wmt)

}
\references{
\strong{Original implementation and simulation for loadings} \cr
Christensen, A. P., & Golino, H. (2021).
On the equivalency of factor and network loadings.
\emph{Behavior Research Methods}, \emph{53}, 1563-1580.

\strong{Preliminary simulation for scores} \cr
Golino, H., Christensen, A. P., Moulder, R., Kim, S., & Boker, S. M. (2021).
Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
\emph{Psychometrika}.

\strong{Revised network loadings} \cr
Christensen, A. P., Golino, H., Abad, F. J., & Garrido, L. E. (2024).
Revised network loadings.
\emph{PsyArXiv}.
}
\author{
Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
}
