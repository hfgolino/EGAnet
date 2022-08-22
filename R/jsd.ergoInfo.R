#' Jensen-Shannon Distance Test for Ergodicity
#'
#' @description Tests the Jensen-Shannon Distance (\code{\link[EGAnet]{jsd}}) of each
#' individual's network structure to the population's network structure. Using a threshold,
#' the proportion of \code{\link[EGAnet]{jsd}} values \strong{greater} than the threshold
#' are computed. If this proportion is \strong{greater} or equal to the maximum allowable
#' proportion of individuals, then the system is determined to be nonergodic. If the
#' proportion is \strong{fewer} than the maximum allowable proportion, then the system
#' is determined to be ergodic
#'
#' @param dynEGA.object  A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.ind.pop}} object that is used to match the arguments of the EII object.
#'
#' @param method Character.
#' Method to compute Jensen-Shannon Distance.
#' Defaults to \code{"spectral"}.
#' Options:
#' 
#' \itemize{
#' 
#' \item{\code{"kld"}}
#' {Uses Kullback-Leibler Divergence}
#' 
#' \item{\code{"spectral"}}
#' {Uses eigenvalues of combinatiorial Laplacian matrix to compute
#' Von Neumann entropy}
#' 
#' }
#'
#' @param threshold Numeric.
#' Sets the threshold of the \code{\link[EGAnet]{jsd}} determined
#' to be \strong{too distant} for ergodicity to hold.
#' Defaults to \code{0.20}.
#' Values can range between 0 and 1
#' 
#' @param max.proportion Numeric.
#' Sets the proportion of \code{\link[EGAnet]{jsd}} values
#' that are \strong{greater} than the specified \code{threshold}.
#' Defaults to \code{0.10}.
#' Values can range between 0 and 1
#'
#' @examples
#' \donttest{# Dynamic EGA individual and population structures
#' dyn1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA[,-c(22)], n.embed = 5, tau = 1,
#'   delta = 1, id = 21, use.derivatives = 1,
#'   model = "glasso", ncores = 2, corr = "pearson"
#' )
#'
#' # JSD Ergodicity Test
#' testing.ergoinfo <- jsd.ergoInfo(
#'   dynEGA.object = dyn1
#' )}
#'
#' @return Returns a list containing:
#'
#' \item{boot.ergoInfo}{The values of the Ergodicity Information Index obtained in the bootstrap}
#'
#' \item{p.value}{The two-sided *p*-value of the bootstrap test for the Ergodicity Information Index.
#' The null hypothesis is that the empirical Ergodicity Information index is equal to the expected value of the EII
#' with small variation in the population structure}
#'
#' \item{effect}{Indicates wheter the empirical EII is greater or less then the bootstrap distribution of EII.}
#'
#' \item{interpretation}{How you can interpret the result of the test in plain English}
#'
#' \item{plot.dist}{Histogram of the bootstrapped ergodicity information index}
#' 
#' \item{methods}{Methods to report for print/summary S3methods and automated Methods section}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @references 
#' Golino, H., Nesselroade, J., & Christensen, A. P. (2022).
#' Toward a psychology of individuals: The ergodicity information index and a bottom-up approach for finding generalizations.
#' \emph{PsyArXiv}.
#'
#' @export
# JSD Ergodicity Test
# Updated 22.08.2022
jsd.ergoInfo <- function(
    dynEGA.object,
    method = c("kld", "spectral"),
    threshold = 0.20,
    max.proportion = 0.10
){
  
  # Check for method
  if(missing(method)){
    method <- "spectral"
  }else{method <- tolower(match.arg(method))}
  
  # Check for appropriate threshold
  if(threshold < 0 | threshold > 1){
    stop(
      paste(
        "The 'threshold' argument must be between 0 and 1.\nThe threshold provided is out of bounds:",
        formatC(threshold, format = "f", flag = 0, digits = 2)
      )
    )
  }
  
  # Check for appropriate proportion
  if(max.proportion < 0 | max.proportion > 1){
    stop(
      paste(
        "The 'max.proportion' argument must be between 0 and 1.\nThe maximum proportion provided is out of bounds:",
        formatC(max.proportion, format = "f", flag = 0, digits = 2)
      )
    )
  }
  
  # Check for class
  if(!is(dynEGA.object, "dynEGA") & !is(dynEGA.object, "dynEGA.ind.pop")){
    stop(
      paste(
        "Input into the `dynEGA.object` argument's class is not `dynEGA` or `dynEGA.ind.pop`.\n\n",
        "Class of dynEGA.object = ", paste(
          class(dynEGA.object), sep = "", collapse = ", "
        ),
        sep = ""
      )
    )
  }else if(is(dynEGA.object, "dynEGA.ind.pop")){
    dynEGA.pop <- dynEGA.object$dynEGA.pop
  }else if(is(dynEGA.object, "dynEGA")){
    dynEGA.pop <- dynEGA.object
  }
  
  # Obtain individual dynEGA objects only
  dynEGA.ind <- dynEGA.object$dynEGA.ind$dynEGA
  
  # Remove methods from dynEGA.ind
  if("methods" %in% tolower(names(dynEGA.ind))){
    dynEGA.ind <- dynEGA.ind[-which(tolower(names(dynEGA.ind)) == "methods")]
  }
  
  # Obtain individual networks
  individuals <- lapply(dynEGA.ind, function(x){
    return(x$network)
  })
  
  # Obtain population network
  population <- dynEGA.pop$dynEGA$network
  
  # Let user know JSDs are being obtained
  message("Computing Jensen-Shannon Distance...", appendLF = FALSE)
  
  # Compute JSD
  jsds <- unlist(
    lapply(individuals, function(x){
      jsd(x, population, method = method)
    })
  )
  
  # Let user know data generation is done
  message("done.")
  
  # Obtain values greater than or equal to threshold
  proportion <- mean(jsds >= threshold)
  
  # Plot
  jsd_df <- as.data.frame(jsds)
  jsd_plot <- suppressWarnings(
    suppressMessages(
      ggpubr::gghistogram(
        jsd_df, x = "jsds",
        add = "mean",
        fill = "#00AFBB",
        color = "black",
        rug = TRUE,
        ylab = "Frequency",
        xlab = "Jensen-Shannon Distance"
      ) +
        ggplot2::geom_vline(
          xintercept = threshold, color = "#00AFBB", linetype = "dotted"
        )
    )
  )
  
  # Set effect and interpretation
  ## Effect
  effect <- ifelse(
    proportion >= max.proportion,
    "nonergodic",
    "ergodic"
  )
  
  ## Interpretation
  interpretation <- switch(
    effect,
    "nonergodic" = paste(
      "The proportion, ",
      formatC(proportion, format = "f", flag = 0, digits = 2),
      ", of individuals with JSDs greater than ",
      formatC(threshold, format = "f", flag = 0, digits = 2),
      " was larger than 0.10 of the sample. ",
      "The population structure cannot adequately generalize to all individuals. ",
      "The system is nonergodic.",
      sep = ""
    ),
    "ergodic" = paste(
      "The proportion, ",
      formatC(proportion, format = "f", flag = 0, digits = 2),
      ", of individuals with JSDs greater than ",
      formatC(threshold, format = "f", flag = 0, digits = 2),
      " was less than 0.10 of the sample. ",
      "The population structure can adequately generalize to all individuals. ",
      "The system is ergodic.",
      sep = ""
    )
  )
  
  # Return results
  results <- list(
    jsd.ergoInfo = jsds,
    threshold = threshold,
    proportion = proportion,
    effect = effect,
    interpretation = interpretation,
    plot.dist = jsd_plot
  )
  
  # Set class
  class(results) <- "jsd.ergoInfo"
  
  return(results)
}
#----