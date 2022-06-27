#' The Knife-Jack Test for the Ergodicity Information Index
#'
#' @description Tests the Ergodicity Information Index obtained in the empirical sample with the EII obtained by inverse-Jackkife sampling. In traditional Jackknife
#' sampling, one single observation is omitted in each data sample. Here, we're using the inverse of that, meaning that we analyze one individual at a time. Therefore,
#' since for each individual \code{i} the individual and population structure is the same (because is a single individual structure), the ergodicity information index obtained
#' is what you could expect if every individual in a dataset had the same structure as \code{i}. Thereforem the distribution of EII is the null distribution obtained under the hypothesis
#' that the structure of individuals are uniquely similar to a population structure. Since no individual is exactly equal to one another, there's enough variability in the Knife-Jack distribution.
#' If all individuals have the same structure, then there would be no significant difference between the empirical EII value and the null distribution of EII values obtained by the Knife-Jack approach.
#' However, if individuals have different structures, then the empirical EII value would be significantly different to the null distribution of EII values obtained via Knife-Jaccking.
#' The p-values in the Knife-Jack test can be calculated as \code{(sum(EII>=knifejack.EII)+1)/(N+1)} and as
#' \code{(sum(EII<=knifejack.EII)+1)/(N+1)}, where EII is the empirical Ergodicity Information Index, knifejack.EII is the values of the Ergodicity Information Index obtained
#' in the Knife-Jack samples, and \code{N} is the number of unique individuals in the dataset. The two-sided p-value is computed as two times the lowest p-value. In the Knife-Jack Test for the Ergodicity Information Index,
#' the null hypothesis is that the empirical value of EII is equal to the values of EII obtained in multiple individuals with the same structure as the population structure estimated
#' via \code{\link[EGAnet]{dynEGA}}.
#' Small values of p indicate that is very unlikely to obtain an EII as large as the one obtained in the empirical sample if the null hypothesis is true (i.e. all individuals have the same structure as the population structure), thus there is convincing evidence that the empirical Ergodicity Information Index is
#' different than it could be expected if all individuals had a similar latent structure.
#'
#' @param dynEGA.object  A \code{\link[EGAnet]{dynEGA}} or a
#' \code{\link[EGAnet]{dynEGA.pop.ind}} object that is used to match the arguments of the EII object.
#'
#' @param EII A \code{\link[EGAnet]{ergoInfo}} object, used to estimate the Empirical Ergodicity Information Index, or the estimated value of EII estimated
#' using the \code{\link[EGAnet]{ergoInfo}} function.
#'
#' @param use Character.
#' A string indicating what network element will be used to compute the algorithm complexity in the \code{\link[EGAnet]{ergoInfo}} function,
#' the list of edges or the weights of the network.
#' Defaults to \code{use = "edge.list"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{edge.list}}}
#' {Calculates the algorithm complexity using the list of edges.}
#'
#' \item{\strong{\code{weights}}}
#' {Calculates the algorithm complexity using the weights of the network.}
#' }
#'
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
#'
#' \dontrun{
#' \donttest{
#' dyn1 <- dynEGA.ind.pop(data = sim.dynEGA[,-c(22)], n.embed = 5, tau = 1,
#'                       delta = 1, id = 21, use.derivatives = 1,
#'                     model = "glasso", ncores = 2, corr = "pearson")
#'
#' eii1 <- ergoInfo(dynEGA.object = dyn1)$EII
#'
#' testing.ergoinfo <- knifejack.ergoInfo(
#'   dynEGA.object = dyn1,
#'   EII = eii1,
#'   use = "edge.list",
#'   ncores = 2)
#' }}
#'
#' @return Returns a list containing:
#'
#' \item{knifejack.ergoInfo}{The values of the Ergodicity Information Index obtained in the Knife-Jack sampling}
#'
#' \item{p.value.twosided}{The p-value of the Monte-Carlo test for the Ergodicity Information Index.
#' The null hypothesis is that the empirical Ergodicity Information index is equal to the expected value of the EII if the all individuals
#' had similar latent structures.}
#'
#' \item{effect}{Indicates wheter the empirical EII is greater or less then the Monte-Carlo obtained EII.}
#'
#' \item{plot.dist}{Histogram of the bootstrapped ergodicity information index}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
# Knife-Jack Test for the Ergodicity Information Index
# Updated 27.06.2022
knifejack.ergoInfo <- function(
    dynEGA.object = NULL,
    EII, use = c("edge.list", "weights"),
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

  # Obtain methods
  methods <- dynEGA.pop$dynEGA$Methods

  # Match arguments
  ## GLLA
  n.embed <- methods$glla$n.embed
  tau <- methods$glla$tau
  delta <- methods$glla$delta
  use.derivatives <- methods$glla$derivatives
  ## Model
  model <- methods$EGA$model
  model.args <- methods$EGA$model.args
  algorithm <- methods$EGA$algorithm
  algorithm.args <- methods$EGA$algorithm.args
  corr <- methods$EGA$corr

  # Check for EII
  if(missing(EII)){

    # Check for dynEGA.ind.pop object
    if(is(dynEGA.object, "dynEGA.ind.pop")){

      # Let user know that EII will be computed
      message("`EII` argument is missing. Computing EII...", appendLF = FALSE)

      # Compute EII
      EII <- ergoInfo(
        dynEGA.object = dynEGA.object,
        use = "edge.list"
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
  
  # Obtain IDs
  IDs <- names(dynEGA.ind)
  
  # Set up knifejack.data to be consistent with `dyn.ind.pop` output
  knifejack.data <- lapply(seq_along(dynEGA.ind), function(i){
    
    # Initialize list
    dynEGA.object <- list()
    dynEGA.object$dynEGA.pop <- list()
    dynEGA.object$dynEGA.ind <- list()
    dynEGA.object$dynEGA.pop$dynEGA <- dynEGA.ind[[i]]
    dynEGA.object$dynEGA.ind$dynEGA[[IDs[i]]] <- dynEGA.ind[[i]]
    class(dynEGA.object) <- "dynEGA.ind.pop"
    return(dynEGA.object)
    
  })
  
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
    X = knifejack.data, cl = cl,
    FUN = ergoInfo, use = use
  )

  # Stop cluster
  parallel::stopCluster(cl)

  complexity.estimates2 <- unlist(
    lapply(
      complexity.estimates, function(x){
        x$EII
      }
    )
  )

  ## Compute the P-value of the knifejackstrap test:
  # N <- length(dynEGA.ind)
  # p.greater <- (sum(EII>=complexity.estimates2)+1)/(N +1)
  # p.lower <- (sum(EII<=complexity.estimates2)+1)/(N +1)
  # p.values <- c(p.greater, p.lower)
  # two.sided <- 2*min(p.values)
  two.sided <- mean(abs(complexity.estimates2) >= abs(EII), na.rm = TRUE)

  # Plot:
  complexity.df <- data.frame(EII = complexity.estimates2, ID = 1:length(complexity.estimates2))
  plot.knifejackErgoInfo <- suppressWarnings(suppressMessages(
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
  results$knifejack.ergoInfo <- complexity.estimates2
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
  
  results$plot.dist <- plot.knifejackErgoInfo
  return(results)
}
#----
