#' EGAnet--package
#' 
#' @description Implements the Exploratory Graph Analysis (\code{\link[EGAnet]{EGA}}; Golino & Epskamp, 2017; Golino, Shi, et al., 2020) framework for dimensionality and psychometric assessment.
#' EGA is part of a new area called \emph{network psychometrics} that uses undirected network models for the assessment
#' of psychometric properties. EGA estimates the number of dimensions (or factors) using graphical lasso \code{\link[qgraph]{EBICglasso}} or
#' Triangulated Maximally Filtered Graph (\code{\link[NetworkToolbox]{TMFG}}) and a weighted network community detection algorithm (Christensen & Golino, under review). A bootstrap method for
#' verifying the stability of the dimensions and items in those dimensions is available (\code{\link[EGAnet]{bootEGA}}; Christensen & Golino, 2019). The fit of the structure suggested by EGA
#' can be verified using Entropy Fit Indices (\code{\link[EGAnet]{entropyFit}}, \code{\link[EGAnet]{tefi}}; Golino, Moulder, et al., 2020). A novel approach called Unique Variable Analysis (\code{\link[EGAnet]{UVA}}) can be used to
#' identify and reduce redundant variables in multivariate data (Christensen, Garrido, & Golino, under review). Network loadings (\code{\link[EGAnet]{net.loads}}),
#' which are roughly equivalent to factor loadings when the data generating model is a factor model, are available (Christensen & Golino, 2021).
#' Network scores (\code{\link[EGAnet]{net.scores}}) can also be computed using the network loadings. Finally, dynamic EGA (\code{\link[EGAnet]{dynEGA}}) will estimate dimensions from time series data for individual, group, and sample levels (Golino, Christensen, et al., under review).
#'
#' @references
#' Christensen, A. P., Garrido, L. E., & Golino, H. (under review).
#' Unique Variable Analysis: A novel approach to detect redundant variables in multivariate data.
#' \emph{PsyArXiv}.
#' \doi{10.31234/osf.io/4kra2} \cr
#' # Related functions: \code{\link[EGAnet]{UVA}}
#' 
#' Christensen, A. P., Garrido, L. E., & Golino, H. (2021).
#' Comparing community detection algorithms in psychological data: A Monte Carlo simulation.
#' \emph{PsyArXiv}.
#' \doi{10.31234/osf.io/hz89e} \cr
#' # Related functions: \code{\link[EGAnet]{EGA}}
#' 
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{Psych}.
#' \doi{10.31234/osf.io/9deay} \cr
#' # Related functions: \code{\link[EGAnet]{bootEGA}}, \code{\link[EGAnet]{dimensionStability}},
#' # and \code{\link[EGAnet]{itemStability}}
#' 
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}, \emph{53}, 1563-1580.
#' \doi{10.3758/s13428-020-01500-6} \cr
#' # Related functions: \code{\link[EGAnet]{LCT}} and \code{\link[EGAnet]{net.loads}}
#' 
#' Christensen, A. P., & Golino, H. (2021).
#' Factor or network model? Predictions from neural networks.
#' \emph{Journal of Behavioral Data Science}, \emph{1}(1), 85-126.
#' \doi{10.35566/jbds/v1n1/p5} \cr
#' # Related functions: \code{\link[EGAnet]{LCT}}
#' 
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}, 1095-1108.
#' \doi{10.1002/per.2265} \cr
#' # Related functions: \code{\link[EGAnet]{bootEGA}}, \code{\link[EGAnet]{dimensionStability}},
#' # \code{\link[EGAnet]{EGA}}, \code{\link[EGAnet]{itemStability}}, and \code{\link[EGAnet]{UVA}}
#' 
#' Golino, H., Christensen, A. P., Moulder, R. G., Kim, S., & Boker, S. M. (under review).
#' Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
#' \emph{PsyArXiv}.
#' \doi{10.31234/osf.io/tfs7c} \cr
#' # Related functions: \code{\link[EGAnet]{dynEGA}} and \code{\link[EGAnet]{simDFM}}
#' 
#' Golino, H., & Demetriou, A. (2017).
#' Estimating the dimensionality of intelligence like data using Exploratory Graph Analysis.
#' \emph{Intelligence}, \emph{62}, 54-70.
#' \doi{10.1016/j.intell.2017.02.007} \cr
#' # Related functions: \code{\link[EGAnet]{EGA}}
#' 
#' Golino, H., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PLoS ONE}, \emph{12}, e0174035.
#' \doi{10.1371/journal.pone.0174035} \cr
#' # Related functions: \code{\link[EGAnet]{EGA}}
#' 
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#' \doi{10.1080/00273171.2020.1779642} \cr
#' # Related functions: \code{\link[EGAnet]{entropyFit}}, \code{\link[EGAnet]{tefi}}, and \code{\link[EGAnet]{vn.entropy}} 
#'
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., Thiyagarajan, J. A., & Martinez-Molina, A. (2020).
#' Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors:
#' A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}, 292-320. 
#' \doi{10.1037/met0000255} \cr
#' # Related functions: \code{\link[EGAnet]{EGA}}
#' 
#' Golino, H., Thiyagarajan, J. A., Sadana, M., Teles, M., Christensen, A. P., & Boker, S. M. (under review).
#' Investigating the broad domains of intrinsic capacity, functional ability, and environment:
#' An exploratory graph analysis approach for improving analytical methodologies for measuring healthy aging.
#' \emph{PsyArXiv}.
#' \doi{10.31234/osf.io/hj5mc} \cr
#' # Related functions: \code{\link[EGAnet]{EGA.fit}} and \code{\link[EGAnet]{tefi}}
#'
#' @author Hudson Golino <hfg9s@virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#' 
#' @importFrom utils packageDescription
#'
"_PACKAGE"
#> [1] "_PACKAGE"
#EGAnet----
