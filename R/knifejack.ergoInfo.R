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
#' @param data A data frame with the variables to be used in the analysis.
#' The data frame should be in a long format (i.e. observations for the
#' same individual (for example, individual 1) are placed in order,
#' from time 1 to time t, followed by the observations from individual 2, also
#' ordered from time 1 to time t.)
#'
#' @param id Numeric.
#' Number of the column identifying each individual.
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
#' @param ncores2 Numeric.
#' Number of cores to use in computing results of the \code{dyn.pop.ind} function.
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
#' eii1 <- ergoInfo(data = dyn1)$EII
#'
#' testing.ergoinfo <- knifejack.ergoInfo(data = sim.dynEGA[,-c(22)],
#' id = 21,
#' dynEGA.object = dyn1,
#'  EII = eii1, ncores = 2, ncores2 = 2)
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
# Updated 25.06.2022
knifejack.ergoInfo <- function(
    data,
    id = NULL,
    dynEGA.object = NULL,
    EII, use = c("edge.list", "weights"),
    ncores, ncores2
){

  # MISSING ARGUMENTS HANDLING

  if(missing(dynEGA.object))
  {stop("The 'dynEGA.object' argument is missing! \n The 'dynEGA.object' must be provided to it to the arguments of the EII object.")
  }else{dynEGA.object <- dynEGA.object}

  if(missing(id))
  {stop("The 'id' argument is missing! \n The number of the column identifying each individual must be provided!")
  }else{id <- id}

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
        use = use
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
    else if(!is(EII, "EII")){
    EII <- EII$EII
  }

  if(missing(ncores)){
    ncores <- ceiling(parallel::detectCores() / 2)
  }

  if(missing(ncores2)){
    ncores2 <- ceiling(parallel::detectCores() / 2)
  }

  # Split based on ID
  data.ind <- split(data, f = data[,id])

  # Initialize data list:
  knifejack.data <- vector("list", length = length(data.ind))

  # Let user know data generation has started
  message("\nGetting the Knife-Jack samples and Estimating the Population and Individual Structures...(be patient, it takes time)")


  # Set parallelization
  cl2 <- parallel::makeCluster(ncores)

  # Export to cluster
  parallel::clusterExport(
    cl = cl2,
    varlist = c(
      "data.ind","knifejack.data"),
    envir = environment()
  )

  # Estimate population and individual structures, but for each id
  knifejack.data <- pbapply::pblapply(
    cl = cl2,
    X = data.ind,
    FUN = function(i){
     dynEGA.ind.pop(data = i,
                             n.embed = n.embed, # number of embeddings
                             id = id,
                             delta = delta,
                             use.derivatives = use.derivatives, # derivatives to use
                             model = modle, # network estimation method
                             model.args = model.args,
                             algorithm = algorithm,
                             algorithm.args = algorithm.args,
                             corr = corr,
                             ncores = ncores2 # processing cores
      )
    }

  )


  # Compute EII
  complexity.estimates <- pbapply::pblapply(X = knifejack.data, cl = cl2,
                                            FUN = ergoInfo, use = use)

  # Stop cluster
  parallel::stopCluster(cl)

  #let user know results are being computed
  message("Computing results...\n")

  complexity.estimates2 <- unlist(
    lapply(
      complexity.estimates, function(x){
        x$EII
      }
    )
  )

  ## Compute the P-value of the knifejackstrap test:
  N <- length(unique(data[,id]))
  p.greater <- (sum(EII>=complexity.estimates2)+1)/(N +1)
  p.lower <- (sum(EII<=complexity.estimates2)+1)/(N +1)
  p.values <- c(p.greater, p.lower)
  two.sided <- 2*min(p.values)

  # Plot:
  complexity.df <- data.frame(EII = complexity.estimates2, ID = 1:length(complexity.estimates2))
  plot.knifejackErgoInfo <- suppressWarnings(suppressMessages(ggpubr::gghistogram(complexity.df, x = "EII",
                                                                             add = "mean",
                                                                             fill = "#00AFBB",
                                                                             color = "black",
                                                                             rug = TRUE,
                                                                             ylab = "Frequency",
                                                                             xlab = "Ergodicity Information Index")+
                                                           ggplot2::geom_vline(xintercept = EII, color = "#00AFBB", linetype = "dotted")))

  ## Return Results:
  results <- vector("list")
  results$knifejack.ergoInfo <- complexity.estimates2
  results$p.value.twosided <- two.sided
  results$effect <- ifelse(p.greater<p.greater, "Greater", "Less")
  results$plot.dist <- plot.knifejackErgoInfo
  return(results)
}
#----
