#----------------------------------------------
## Datasets for EGAnet // Updated 25.03.2020
#----------------------------------------------

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
#' \doi{10.1037//0022-3514.67.6.1063}
#'
#' @examples
#' data("optimism")
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
