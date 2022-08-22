#' Bootstrap Test for the Ergodicity Information Index
#'
#' @description Tests the Ergodicity Information Index obtained in the empirical sample with a distribution of EII 
#' obtained by bootstrap sampling. In traditional bootstrap sampling, individual participants are resampled with
#' replacement from the empirical sample. This process is time consuming when carried out across \emph{v} number
#' of variables, \emph{n} number of participants, \emph{t} number of time points, and \emph{i} number of iterations.
#' 
#' A more efficient process, the approach applied here, is to obtain a sampling distribution of EII values as if
#' all participants in the data have the population network structure. Sampling is not perfect and therefore
#' random noise is added to the edges of the population structure to simulate sampling variability. This noise
#' follows a random uniform distribution ranging from -0.10 to 0.10. In addition, a proportion of edges are
#' rewired to allow for slight variations on the population structure. The proportion of nodes that are rewired
#' is sampled from a random uniform distribution between 0.10 to 0.40. This process is carried out for each
#' participant resulting in \emph{n} variations of the population structure. Afterward, EII is computed. This
#' process is carried out for \emph{i} iterations (e.g., 100).
#' 
#' The result is a sampling distribution of EII values that would be expected if the process was ergodic. If
#' the empirical EII value is significantly less than the distribution or not significantly different, then 
#' the empirical data can be expected to be generated from an ergodic process and the population structure is 
#' sufficient to describe all individuals. If the empirical EII value is significantly greater than the distribution,
#' then the empirical data cannot be described by the population structure -- significant information is lost when
#' collapsing across to the population structure.
#'
#' @param dynEGA.object  A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.ind.pop}} object that is used to match the arguments of the EII object.
#'
#' @param EII A \code{\link[EGAnet]{ergoInfo}} object, used to estimate the Empirical Ergodicity Information Index, or the estimated value of EII estimated
#' using the \code{\link[EGAnet]{ergoInfo}} function. Inherits \code{use} from \code{\link[EGAnet]{ergoInfo}}
#'
#' @param iter Numeric integer.
#' Number of replica samples to generate from the bootstrap analysis.
#' At least \code{100} is recommended
#'
#' @param ncores Numeric.
#' Number of cores to use in computing results.
#' Defaults to \code{parallel::detectCores() / 2} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing.
#' Recommended to use maximum number of cores minus one
#'
#' If you're unsure how many cores your computer has,
#' then use the following code: \code{parallel::detectCores()}
#'
#' @examples
#' \donttest{# Dynamic EGA individual and population structures
#' dyn1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA[,-c(22)], n.embed = 5, tau = 1,
#'   delta = 1, id = 21, use.derivatives = 1,
#'   model = "glasso", ncores = 2, corr = "pearson"
#' )
#'
#' # Empirical Ergodicity Information Index
#' eii1 <- ergoInfo(dynEGA.object = dyn1, use = "edge.list")
#'
#' # Bootstrap Test for Ergodicity Information Index
#' testing.ergoinfo <- boot.ergoInfo(
#'   dynEGA.object = dyn1, EII = eii1,
#'   ncores = 2
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
# Bootstrap Test for the Ergodicity Information Index
# Updated 29.07.2022
boot.ergoInfo <- function(
    dynEGA.object,
    EII, iter = 100,
    ncores
){
  
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
  
  # Check for EII
  if(missing(EII)){
    
    # Check for dynEGA.ind.pop object
    if(is(dynEGA.object, "dynEGA.ind.pop")){
      
      # Let user know that EII will be computed
      message("`EII` argument is missing. Computing EII...", appendLF = FALSE)
      
      # Compute EII
      EII <- ergoInfo(
        dynEGA.object = dynEGA.object,
        use = "edge.list" # "weighted"
      )
      
      # Let user know EII is complete
      message("done", appendLF = TRUE)
      
    }else{
      
      # Let user know that EII cannot be computed
      stop("`EII` argument is missing. EII cannot be computed without a `dynEGA.ind.pop` object class as input to `dynEGA.object`", appendLF = FALSE)
      
    }
    
  }
  
  # Check for class of EII
  if(is(EII, "EII")){
    use <- EII$use
    EII <- EII$EII
  }
  
  if(missing(ncores)){
    ncores <- ceiling(parallel::detectCores() / 2)
  }
  
  # Obtain individual dynEGA objects only
  dynEGA.ind <- dynEGA.object$dynEGA.ind$dynEGA
  
  # Remove methods from dynEGA.ind
  if("methods" %in% tolower(names(dynEGA.ind))){
    dynEGA.ind <- dynEGA.ind[-which(tolower(names(dynEGA.ind)) == "methods")]
  }
  
  # Replace individual networks with population networks
  dynEGA.ind <- lapply(dynEGA.ind, function(x){
    x$network <- dynEGA.pop$dynEGA$network
    return(x)
  })
  
  # Let user know data are being generated
  message("Generating rewired population networks...", appendLF = FALSE)
  
  # Set up data to be consistent with `dyn.ind.pop` output
  boot.data <- lapply(1:iter, function(i){
    
    # Initialize list
    dynEGA.object <- list()
    dynEGA.object$dynEGA.pop <- list()
    dynEGA.object$dynEGA.ind <- list()
    dynEGA.object$dynEGA.pop <- dynEGA.pop
    dynEGA.object$dynEGA.ind$dynEGA <- dynEGA.ind
    dynEGA.object$dynEGA.ind$dynEGA <- lapply(
      1:length(dynEGA.ind), # same number of individuals
      function(i){
        dynEGA.ind[[i]]$network <- rewire(
          dynEGA.pop$dynEGA$network,
          min = 0.20, max = 0.40,
          noise = 0.10
        )
        return(dynEGA.ind[[i]])
      }
    )
    class(dynEGA.object) <- "dynEGA.ind.pop"
    return(dynEGA.object)
    
  })
  
  # Let user know data generation is done
  message("done.")
  
  # Let user know results are being computed
  message("Computing results...")
  
  # Set up parallelization
  cl <- parallel::makeCluster(ncores)
  
  # Export prime numbers
  parallel::clusterExport(
    cl = cl,
    varlist = "prime.num"
  )
  
  # Compute EII
  complexity.estimates <- pbapply::pblapply(
    X = boot.data, cl = cl,
    FUN = ergoInfo, use = use
  )
  
  # Stop cluster
  parallel::stopCluster(cl)
  
  # Obtain EII values
  complexity.estimates2 <- unlist(
    lapply(
      complexity.estimates, function(x){
        x$EII
      }
    )
  )
  
  ## Compute the p-value of the bootstrap test:
  N <- iter # length(dynEGA.ind)
  p.greater <- (sum(EII>=complexity.estimates2)+1)/(N +1)
  p.lower <- (sum(EII<=complexity.estimates2)+1)/(N +1)
  p.values <- c(p.greater, p.lower)
  two.sided <- 2*min(p.values)
  #one.sided <- (sum(abs(complexity.estimates2) >= abs(EII), na.rm = TRUE)+1)/(N + 1)
  
  # Plot:
  complexity.df <- data.frame(EII = complexity.estimates2, ID = 1:length(complexity.estimates2))
  plot.bootErgoInfo <- suppressWarnings(suppressMessages(
    ggpubr::gghistogram(complexity.df, x = "EII",
                        add = "mean",
                        fill = "#00AFBB",
                        color = "black",
                        rug = TRUE,
                        ylab = "Frequency",
                        xlab = "Ergodicity Information Index") +
      ggplot2::geom_vline(xintercept = EII, color = "#00AFBB", linetype = "dotted"))
  )
  
  ## Return Results:
  results <- vector("list")
  results$boot.ergoInfo <- complexity.estimates2
  results$p.value <- two.sided
  effect <- ifelse(two.sided <= .05, mean(
    EII >= complexity.estimates2, na.rm = TRUE
  ), "n.s.")
  results$effect <- ifelse(
    effect == "n.s.", "n.s.",
    ifelse(
      effect > .50, "greater", "less"
    )
  )
  interpretation <- switch(
    results$effect,
    "n.s." = "The empirical EII was not different from what would be expected from random variation in the population structure, meaning non-significant information is lost when aggregating the results into a single, population network.",
    "less" = "The empirical EII was less than what would be expected from random variation in the population structure, meaning non-significant information is lost when aggregating the results into a single, population network.",
    "greater" = "The empirical EII was greater than what would be expected from random variation in the population structure, meaning significant information is lost when aggregating the results into a single, population network."
  )
  
  results$interpretation <- interpretation
  
  results$plot.dist <- plot.bootErgoInfo
  
  # Prepare methods
  methods <- list()
  methods$use <- use
  methods$iter <- iter
  
  results$Methods <- methods
  
  # Set class
  class(results) <- "boot.ergoInfo"
  
  return(results)
}
#----