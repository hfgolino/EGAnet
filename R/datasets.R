#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Datasets for EGAnet // Updated 24.02.2022
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Bootstrap EGA WMT-2 Data----
#' Bootstrap EGA WMT-2 Data
#'
#' \code{\link[EGAnet]{bootEGA}} Results of \code{\link[EGAnet]{wmt2}}Data
#'
#' \code{\link[EGAnet]{bootEGA}} results using the \code{"glasso"} model and 500 iterations
#' of the Wiener Matrizen-Test 2 (WMT-2)
#'
#' @name boot.wmt
#'
#' @docType data
#'
#' @usage data(boot.wmt)
#'
#' @format A list with 8 objects (see \code{\link[EGAnet]{bootEGA}})
#'
#' @keywords datasets
#'
#' @examples
#' data("boot.wmt")
#'
NULL

#Depression Data----
#' Depression Data
#'
#' A response matrix (n = 574) of the Beck Depression Inventory, Beck Anxiety Inventory and
#' the Athens Insomnia Scale.
#'
#' @name depression
#'
#' @docType data
#'
#' @usage data(depression)
#'
#' @format A 574x78 response matrix
#'
#' @keywords datasets
#'
#' @examples
#' data("depression")
#'
NULL

#' Loadings Comparison Test Deep Learning Neural Network Weights
#'
#' A list of weights from four different neural network models:
#' random vs. non-random model (\code{r_nr_weights}),
#' low correlation factor vs. network model (\code{lf_n_weights}),
#' high correlation with variables less than or equal to factors vs. network model (\code{hlf_n_weights}), and
#' high correlation with variables greater than factors vs. network model (\code{hgf_n_weights})
#'
#' @name dnn.weights
#'
#' @docType data
#'
#' @usage data(dnn.weights)
#'
#' @format A list of with a length of 4
#'
#' @keywords datasets
#'
#' @examples
#' data("dnn.weights")
#'
NULL

#EGA WMT-2 Data----
#' EGA WMT-2 Data
#'
#' \code{\link[EGAnet]{EGA}} Network of \code{\link[EGAnet]{wmt2}}Data
#'
#' An \code{\link[EGAnet]{EGA}} using the \code{"glasso"} model of the
#' Wiener Matrizen-Test 2 (WMT-2)
#'
#' @name ega.wmt
#'
#' @docType data
#'
#' @usage data(ega.wmt)
#'
#' @format A 17 x 17 adjacency matrix
#'
#' @keywords datasets
#'
#' @examples
#' data("ega.wmt")
#'
NULL

#Intelligence Data----
#' Intelligence Data
#'
#' A response matrix (n = 1152) of the International Cognitive Ability Resource (ICAR)
#' intelligence battery developed by Condon and Revelle (2016).
#'
#' @name intelligenceBattery
#'
#' @docType data
#'
#' @usage data(intelligenceBattery)
#'
#' @format A 1185x125 response matrix
#'
#' @keywords datasets
#'
#' @examples
#' data("intelligenceBattery")
#'
NULL

#Optimism Data----
#' Optimism Data
#'
#' A response matrix (n = 282) containing responses to 10 items of the Revised Life
#' Orientation Test (LOT-R), developed by Scheier, Carver, & Bridges (1994).
#'
#' @name optimism
#'
#' @docType data
#'
#' @usage data(optimism)
#'
#' @format A 282x10 response matrix
#'
#' @keywords datasets
#'
#' @references
#' Scheier, M. F., Carver, C. S., & Bridges, M. W. (1994).
#' Distinguishing optimism from neuroticism (and trait anxiety, self-mastery, and self-esteem): a reevaluation of the Life Orientation Test.
#' \emph{Journal of Personality and Social Psychology}, \emph{67}, 1063-1078.
#'
#' @examples
#' data("optimism")
#'
NULL

#Prime Numbers ----
#' Prime Numbers through 100,000
#'
#' Numeric vector of primes generated from the primes package. Used in
#' the function \code{[EGAnet]{ergoInfo}}. Not for general use
#'
#' @name prime.num
#'
#' @docType data
#'
#' @usage data(prime.num)
#'
#' @format A 1185x24 response matrix
#'
#' @keywords datasets
#'
#' @examples
#' data("prime.num")
#'
NULL

#sim.dynEGA Data ----
#' sim.dynEGA Data
#'
#' A simulated (multivariate time series) data with 20 variables, 200 individual observations,
#' 50 time points per individual and 2 groups of individuals.
#'
#' @name sim.dynEGA
#'
#' @docType data
#'
#' @usage data(sim.dynEGA)
#'
#' @format A 10000x22 multivariate time series
#'
#' @keywords datasets
#'
#' @examples
#' data("sim.dynEGA")
#'
NULL

#Toy Example Data----
#' Toy Example Data
#'
#' A simulated dataset with 2 factors, three items per factor and  n = 500.
#'
#' @name toy.example
#'
#' @docType data
#'
#' @usage data(toy.example)
#'
#' @format A 500x6 response matrix
#'
#' @keywords datasets
#'
#' @examples
#' data("toy.example")
#'
NULL

#WMT-2 Data----
#' WMT-2 Data
#'
#' A response matrix (n = 1185) of the Wiener Matrizen-Test 2 (WMT-2).
#'
#' @name wmt2
#'
#' @docType data
#'
#' @usage data(wmt2)
#'
#' @format A 1185x24 response matrix
#'
#' @keywords datasets
#'
#' @examples
#' data("wmt2")
#'
NULL
