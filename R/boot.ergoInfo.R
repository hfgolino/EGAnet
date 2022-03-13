#' Bootstrap Test for the Ergodicity Information Index
#'
#' @description Computes a parametric Bootstrap Test for the Ergodicity Information Index, comparing the
#' empirical Ergodicity Information index to values obtained in data generated using \code{N} parametric bootstraps of the correlation matrix estimated using the
#' \code{\link[EGAnet]{dynEGA}} function, for the population structure. The p-values in the bootstrap test can be calculated as \code{(sum(EII>=boot.EII)+1)/(iter+1)} and as
#' \code{(sum(EII<=boot.EII)+1)/(iter+1)}, where EII is the empirical Ergodicity Information Index, boot.EII is the values of the Ergodicity Information Index obtained
#' in the bootstraped samples, and \code{iter} is the number of random samples generated in the simulation. The two-sided p-value is computed as two times the lowest p-value. In the bootstrap Test for the Ergodicity Information Index,
#' the null hypothesis is that the empirical value of EII is equal to the values of EII obtained in multiple individuals with the same structure as the population structure estimated
#' via \code{\link[EGAnet]{dynEGA}}.
#' Small values of p indicate that is very unlikely to obtain an EII as large as the one obtained in the empirical sample if the null hypothesis is true (i.e. all individuals have the same structure as the population structure), thus there is convincing evidence that the empirical Ergodicity Information Index is
#' different than it could be expected if all individuals had a similar latent structure.
#'
#' @param dynEGA.pop A dynEGA or a dynEGA.pop.ind object.
#'
#' @param iter Numeric integer.
#' Number of random samples to generate in the Monte-Carlo simulation.
#' At least \code{500} is recommended
#'
#' @param EII Numeric.
#' Empirical Ergodicity Information Index obtained via the \code{\link[EGAnet]{ergoInfo}} function.
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
#' @param embed Integer.
#' Number of embedded dimensions (the number of observations to be used in the \code{\link[EGAnet]{Embed}} function). For example,
#' an \code{"embed = 5"} will use five consecutive observations to estimate a single derivative.
#' Default is \code{"embed = 5"}.
#'
#' @param tau Integer.
#' Number of observations to offset successive embeddings in the \code{\link[EGAnet]{Embed}} function. A tau of one uses adjacent observations.
#' Default is \code{"tau = 1"}.
#'
#' @param delta Integer.
#' The time between successive observations in the time series.
#' Default is \code{"delta = 1"}.
#'
#' @param derivatives Integer.
#' The order of the derivative to be used in the EGA procedure. Default to 1.

#' @param corr Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{cor_auto}}}
#' {Computes the correlation matrix using the \code{\link[qgraph]{cor_auto}} function from
#' \code{\link[qgraph]{qgraph}}}.
#'
#' \item{\strong{\code{pearson}}}
#' {Computes Pearson's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#'
#' \item{\strong{\code{spearman}}}
#' {Computes Spearman's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#' }
#'
#' @param model Character.
#' A string indicating the method to use. Defaults to \code{glasso}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter.
#' This is the default method}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#'
#' }
#'
#' @param model.args List.
#' A list of additional arguments for \code{\link[EGAnet]{EBICglasso.qgraph}}
#' or \code{\link[EGAnet]{TMFG}}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#'
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
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
#'
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
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
#' testing.ergoinfo <- boot.ergoInfo(dynEGA.pop = dyn1, iter = 10,EII = eii1,
#' embed = 5, tau = 1, delta = 1, derivatives = 1,
#' model = "glasso", ncores = 2, corr = "pearson")
#' }}
#'
#' @return Returns a list containing:
#'
#' \item{boot.ergoInfo}{The values of the Ergodicity Information Index obtained in the Monte-Carlo Simulation}
#'
#' \item{p.value.twosided}{The p-value of the Monte-Carlo test for the Ergodicity Information Index.
#' The null hypothesis is that the empirical Ergodicity Information index is equal to the expected value of the EII if the all individuals
#' had similar latent structures.}
#'
#' \item{effect}{Indicates wheter the empirical EII is greater or less then the Monte-Carlo obtained EII.}
#'
#' \item{plot.dist}{Histogram of the bootstrapped ergodicity information index}
#'
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#'
#' @export
# Bootstrap Test for the Ergodicity Information Index
# Updated 29.12.2021

boot.ergoInfo <- function(dynEGA.pop,
                          iter,
                          EII,
                          use,
                          embed,
                          tau,
                          delta,
                          derivatives,
                          model, model.args = list(),
                          algorithm = c("walktrap", "louvain"),
                          algorithm.args = list(),
                          corr, ncores, ...
){


  #### MISSING ARGUMENTS HANDLING ####

  if(missing(dynEGA.pop))
  {  # Warning
    warning(
      "The 'dynEGA.pop' argument is missing. Please, provide the name of the object
        created using the 'dynEGA' or the 'dynEGA.pop.ind' functions."
    )
  }


  if(missing(use))
  {use <- "edge.list"
  }else{use}

  if(missing(embed))
  {embed <- 5
  }else{embed}

  if(missing(tau))
  {tau <- 1
  }else{tau}

  if(missing(delta))
  {delta <- 1
  }else{delta}

  if(missing(derivatives))
  {derivatives <- 1
  }else{derivatives}

  if(missing(model))
  {model <- "glasso"
  }else{model}

  if(missing(corr))
  {corr <- "cor_auto"
  }else{corr}

  if(missing(algorithm))
  {algorithm <- "walktrap"
  }else{algorithm}

  if(missing(ncores))
  {ncores <- ceiling(parallel::detectCores() / 2)
  }else{ncores}



  # Initialize Data list
  data.sim <- vector("list", length = iter)

  #let user know data generation has started
  message("\nGenerating the Data...\n", appendLF = FALSE)

  # Get derivative estimates (and associated information)
  derivative_estimates <- dynEGA.pop$Derivatives$EstimatesDF
  N <- nrow(derivative_estimates)
  IDs <- derivative_estimates[,ncol(derivative_estimates)]
  unique.ids <- unique(IDs)
  time.points <- table(IDs)
  time.points <- time.points+(embed-1)
  # Initialize Data list
  data.sim <- vector("list", length = iter)
  for(i in 1:iter){
    data.sim[[i]] <- vector("list", length = length(unique.ids))
  }

  if(class(dynEGA.pop)=="dynEGA"){
    for(i in 1:iter){
      for(j in 1:length(unique.ids)){
        data.sim[[i]][[j]] <- MASS_mvrnorm(n = time.points[[j]], mu = rep(0, ncol(dynEGA.pop$dynEGA$cor.dat)), Sigma = as.matrix(Matrix::nearPD(solve(dynEGA.pop$network))$mat))
      }
    }
  } else if(class(dynEGA.pop)=="dynEGA.ind.pop"){
    for(i in 1:iter){
      for(j in 1:length(unique.ids)){
        data.sim[[i]][[j]] <- as.data.frame(MASS_mvrnorm(n = time.points[[j]], mu = rep(0, ncol(dynEGA.pop$dynEGA.pop$cor.data)), Sigma = as.matrix(Matrix::nearPD(solve(dynEGA.pop$dynEGA.pop$network))$mat)))
        data.sim[[i]][[j]]$ID <- rep(i, each = time.points[[j]])
      }
    }
  }

  data.sim.df <- lapply(
    seq_along(data.sim), function(i){
      as.data.frame(long_results(data.sim[[i]]))
    }
  )

  variab <- ncol(data.sim.df[[1]])-1

  #initialize correlation matrix list
  boot.data.ids <- vector("list", length = iter)
  list.results.sim <- boot.data.ids

  #let user know data generation has started
  message("\nEstimating the Population and Individual Structures...\n", appendLF = FALSE)
  
  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  # parallel::clusterExport(cl = cl,
  #                         varlist = c("list.results.sim",
  #                                     "complexity.estimates"),
  #                         envir=environment())
  #
  # ^^^ Only necessary when testing outside of package ^^^

  #Compute DynEGA in the population and in the individuals
  boot.data <- pbapply::pblapply(X = data.sim.df, cl = cl,
                                 FUN = EGAnet::dynEGA.ind.pop,
                                 n.embed = embed, tau = tau,
                                 delta = delta, id = variab+1, use.derivatives = derivatives,
                                 algorithm = algorithm, algorithm.args = algorithm.args,
                                 model = model, model.args = model.args, corr = corr)


  #let user know data generation has started
  message("Estimating the Ergodicity Information Index\n", appendLF = FALSE)

  complexity.estimates <- pbapply::pblapply(X = boot.data, cl = cl,
                                            FUN = ergoInfo, use = use)
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

  ## Compute the P-value of the bootstrap test:
  p.greater <- (sum(EII>=complexity.estimates2)+1)/(iter+1)
  p.lower <- (sum(EII<=complexity.estimates2)+1)/(iter+1)
  p.values <- c(p.greater, p.lower)
  two.sided <- 2*min(p.values)

  # Plot:
  complexity.df <- data.frame(EII = complexity.estimates2, ID = 1:length(complexity.estimates2))
  plot.bootErgoInfo <- suppressWarnings(suppressMessages(ggpubr::gghistogram(complexity.df, x = "EII",
                                                                             add = "mean",
                                                                             fill = "#00AFBB",
                                                                             color = "black",
                                                                             rug = TRUE)))

  ## Return Results:
  results <- vector("list")
  results$boot.ergoInfo <- complexity.estimates2
  results$p.value.twosided <- two.sided
  results$effect <- ifelse(p.greater<p.greater, "Greater", "Less")
  results$plot.dist <- plot.bootErgoInfo
  return(results)
}
#----
