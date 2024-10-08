% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynEGA.R
\name{dynEGA}
\alias{dynEGA}
\title{Dynamic Exploratory Graph Analysis}
\usage{
dynEGA(
  data,
  id = NULL,
  group = NULL,
  n.embed = 5,
  tau = 1,
  delta = 1,
  use.derivatives = 1,
  level = c("individual", "group", "population"),
  corr = c("auto", "cor_auto", "pearson", "spearman"),
  na.data = c("pairwise", "listwise"),
  model = c("BGGM", "glasso", "TMFG"),
  algorithm = c("leiden", "louvain", "walktrap"),
  uni.method = c("expand", "LE", "louvain"),
  ncores,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{data}{Matrix or data frame.
Participants and variable should be in long format such that
row \emph{t} represents observations for all variables at time point
\emph{t} for a participant. The next row, \emph{t + 1}, represents
the next measurement occasion for that same participant. The next
participant's data should immediately follow, in the same pattern,
after the previous participant

\code{data} should have an ID variable labeled \code{"ID"}; otherwise, it is
assumed that the data represent the population

For groups, \code{data} should have a Group variable labeled \code{"Group"};
otherwise, it is assumed that there are no groups in \code{data}

Arguments \code{id} and \code{group} can be specified to tell the function
which column in \code{data} it should use as the ID and Group variable, respectively

A measurement occasion variable is not necessary and should be \emph{removed}
from the data before proceeding with the analysis}

\item{id}{Numeric or character (length = 1).
Number or name of the column identifying each individual.
Defaults to \code{NULL}}

\item{group}{Numeric or character (length = 1).
Number of the column identifying group membership.
Defaults to \code{NULL}}

\item{n.embed}{Numeric (length = 1).
Defaults to \code{5}.
Number of embedded dimensions (the number of observations to
be used in the \code{\link[EGAnet]{Embed}} function). For example,
an \code{"n.embed = 5"} will use five consecutive observations
to estimate a single derivative}

\item{tau}{Numeric (length = 1).
Defaults to \code{1}.
Number of observations to offset successive embeddings in
the \code{\link[EGAnet]{Embed}} function.
Generally recommended to leave "as is"}

\item{delta}{Numeric (length = 1).
Defaults to \code{1}.
The time between successive observations in the time series (i.e, lag).
Generally recommended to leave "as is"}

\item{use.derivatives}{Numeric (length = 1).
Defaults to \code{1}.
The order of the derivative to be used in the analysis.
Available options:

\itemize{

\item \code{0} --- No derivatives; consistent with moving average

\item \code{1} --- First-order derivatives; interpreted as "velocity" or
rate of change over time

\item \code{2} --- Second-order derivatives; interpreted as "acceleration" or
rate of the rate of change over time

}

Generally recommended to leave "as is"}

\item{level}{Character vector (up to length of 3).
A character vector indicating which level(s) to estimate:

\itemize{

\item \code{"individual"} --- Estimates \code{\link[EGAnet]{EGA}} for each individual in \code{data}
(intraindividual structure; requires an \code{"ID"} column, see \code{data})

\item \code{"group"} --- Estimates \code{\link[EGAnet]{EGA}} for each group in \code{data}
(group structure; requires a \code{"Group"} column, see \code{data})

\item \code{"population"} --- Estimates \code{\link[EGAnet]{EGA}} across all \code{data}
(interindividual structure)

}}

\item{corr}{Character (length = 1).
Method to compute correlations.
Defaults to \code{"auto"}.
Available options:

\itemize{

\item \code{"auto"} --- Automatically computes appropriate correlations for
the data using Pearson's for continuous, polychoric for ordinal,
tetrachoric for binary, and polyserial/biserial for ordinal/binary with
continuous. To change the number of categories that are considered
ordinal, use \code{ordinal.categories}
(see \code{\link[EGAnet]{polychoric.matrix}} for more details)

\item \code{"cor_auto"} --- Uses \code{\link[qgraph]{cor_auto}} to compute correlations.
Arguments can be passed along to the function

\item \code{"pearson"} --- Pearson's correlation is computed for all
variables regardless of categories

\item \code{"spearman"} --- Spearman's rank-order correlation is computed
for all variables regardless of categories

}

For other similarity measures, compute them first and input them
into \code{data} with the sample size (\code{n})}

\item{na.data}{Character (length = 1).
How should missing data be handled?
Defaults to \code{"pairwise"}.
Available options:

\itemize{

\item \code{"pairwise"} --- Computes correlation for all available cases between
two variables

\item \code{"listwise"} --- Computes correlation for all complete cases in the dataset

}}

\item{model}{Character (length = 1).
Defaults to \code{"glasso"}.
Available options:

\itemize{

\item \code{"BGGM"} --- Computes the Bayesian Gaussian Graphical Model.
Set argument \code{ordinal.categories} to determine
levels allowed for a variable to be considered ordinal.
See \code{?BGGM::estimate} for more details

\item \code{"glasso"} --- Computes the GLASSO with EBIC model selection.
See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details

\item \code{"TMFG"} --- Computes the TMFG method.
See \code{\link[EGAnet]{TMFG}} for more details

}}

\item{algorithm}{Character or
\code{igraph} \code{cluster_*} function (length = 1).
Defaults to \code{"walktrap"}.
Three options are listed below but all are available
(see \code{\link[EGAnet]{community.detection}} for other options):

\itemize{

\item \code{"leiden"} --- See \code{\link[igraph]{cluster_leiden}} for more details

\item \code{"louvain"} --- By default, \code{"louvain"} will implement the Louvain algorithm using
the consensus clustering method (see \code{\link[EGAnet]{community.consensus}}
for more information). This function will implement
\code{consensus.method = "most_common"} and \code{consensus.iter = 1000}
unless specified otherwise

\item \code{"walktrap"} --- See \code{\link[igraph]{cluster_walktrap}} for more details

}}

\item{uni.method}{Character (length = 1).
What unidimensionality method should be used?
Defaults to \code{"louvain"}.
Available options:

\itemize{

\item \code{"expand"} --- Expands the correlation matrix with four variables correlated 0.50.
If number of dimension returns 2 or less in check, then the data
are unidimensional; otherwise, regular EGA with no matrix
expansion is used. This method was used in the Golino et al.'s (2020)
\emph{Psychological Methods} simulation

\item \code{"LE"} --- Applies the Leading Eigenvector algorithm
(\code{\link[igraph]{cluster_leading_eigen}})
on the empirical correlation matrix. If the number of dimensions is 1,
then the Leading Eigenvector solution is used; otherwise, regular EGA
is used. This method was used in the Christensen et al.'s (2023)
\emph{Behavior Research Methods} simulation

\item \code{"louvain"} --- Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
on the empirical correlation matrix. If the number of dimensions is 1,
then the Louvain solution is used; otherwise, regular EGA is used.
This method was validated Christensen's (2022) \emph{PsyArXiv} simulation.
Consensus clustering can be used by specifying either
\code{"consensus.method"} or \code{"consensus.iter"}

}}

\item{ncores}{Numeric (length = 1).
Number of cores to use in computing results.
Defaults to \code{ceiling(parallel::detectCores() / 2)} or half of your
computer's processing power.
Set to \code{1} to not use parallel computing

If you're unsure how many cores your computer has,
then type: \code{parallel::detectCores()}}

\item{verbose}{Boolean (length = 1).
Should progress be displayed?
Defaults to \code{TRUE}.
Set to \code{FALSE} to not display progress}

\item{...}{Additional arguments to be passed on to
\code{\link[EGAnet]{auto.correlate}},
\code{\link[EGAnet]{network.estimation}},
\code{\link[EGAnet]{community.detection}},
\code{\link[EGAnet]{community.consensus}}, and
\code{\link[EGAnet]{EGA}}}
}
\value{
A list containing:

\item{Derivatives}{A list containing:

\itemize{

\item \code{Estimates} --- A list the length of the unique IDs containing
data frames of zero- to second-order derivatives for each ID in \code{data}

\item \code{EstimatesDF} --- A data frame of derivatives across all IDs containing
columns of the zero- to second-order derivatives as well as \code{id} and
\code{group} variables (\code{group} is automatically set to \code{1}
for all if no \code{group} is provided)

}

}

\item{dynEGA}{A list containing:

\itemize{

\item \code{population} --- If \code{level} includes \code{"populaton"}, then
the \code{\link[EGAnet]{EGA}} results for the entire sample

\item \code{group} --- If \code{level} includes \code{"group"}, then
a list containing the \code{\link[EGAnet]{EGA}} results for each \code{group}

\item \code{individual} --- If \code{level} includes \code{"individual"}, then
a list containing the \code{\link[EGAnet]{EGA}} results for each \code{id}

}

}
}
\description{
Estimates dynamic communities in multivariate time series
(e.g., panel data, longitudinal data, intensive longitudinal data) at multiple
time scales and at different levels of analysis:
individuals (intraindividual structure), groups, and population (interindividual structure)
}
\details{
Derivatives for each variable's time series for each participant are
estimated using generalized local linear approximation (see \code{\link[EGAnet]{glla}}).
\code{\link[EGAnet]{EGA}} is then applied to these derivatives to model how variables
are changing together over time. Variables that change together over time are detected
as communities
}
\examples{
# Population structure
simulated_population <- dynEGA(
  data = sim.dynEGA, level = "population"
  # uses simulated data in package
  # useful to understand how data should be structured
)

# Group structure
simulated_group <- dynEGA(
  data = sim.dynEGA, level = "group"
  # uses simulated data in package
  # useful to understand how data should be structured
)

\dontrun{
# Individual structure
simulated_individual <- dynEGA(
  data = sim.dynEGA, level = "individual",
  ncores = 2, # use more for quicker results
  verbose = TRUE # progress bar
)

# Population, group, and individual structure
simulated_all <- dynEGA(
  data = sim.dynEGA,
  level = c("individual", "group", "population"),
  ncores = 2, # use more for quicker results
  verbose = TRUE # progress bar
)

# Plot population
plot(simulated_all$dynEGA$population)

# Plot groups
plot(simulated_all$dynEGA$group)

# Plot individual
plot(simulated_all$dynEGA$individual, id = 1)

# Step through all plots
# Unless `id` is specified, 4 random IDs
# will be drawn from individuals
plot(simulated_all)}

}
\references{
\strong{Generalized local linear approximation} \cr
Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010)
Generalized local linear approximation of derivatives from time series. In S.-M. Chow, E. Ferrer, & F. Hsieh (Eds.),
\emph{The Notre Dame series on quantitative methodology. Statistical methods for modeling human dynamics: An interdisciplinary dialogue},
(p. 161-178). \emph{Routledge/Taylor & Francis Group}.

Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
Using derivative estimates to describe intraindividual variability at multiple time scales.
\emph{Psychological Methods}, \emph{14(4)}, 367-386.

\strong{Original dynamic EGA implementation} \cr
Golino, H., Christensen, A. P., Moulder, R. G., Kim, S., & Boker, S. M. (2021).
Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
\emph{Psychometrika}.

\strong{Time delay embedding procedure} \cr
Savitzky, A., & Golay, M. J. (1964).
Smoothing and differentiation of data by simplified least squares procedures.
\emph{Analytical Chemistry}, \emph{36(8)}, 1627-1639.
}
\seealso{
\code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{EGAnet}
}
\author{
Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
}
