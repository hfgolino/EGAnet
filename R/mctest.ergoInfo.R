#' Monte-Carlo Test for the Ergodicity Information Index
#' @description Computes a Monte-Carlo Test for the Ergodicity Information Index, comparing the
#' empirical Ergodicity Information index to values obtained in a Monte-Carlo simulation in which all individuals
#' have a similar latent structure. The p-values in the Monte-Carlo test can be calculated as \code{(sum(EII>=MC.EII)+1)/(iter+1)} and as
#' \code{(sum(EII<=MC.EII)+1)/(iter+1)}, where EII is the empirical Ergodicity Information Index, MC.EII is the values of the Ergodicity Information Index obtained
#' in the simulation, and \code{iter} is the number of random samples generated in the simulation. The two-sided p-value is computed as two times the lowest p-value. In the Monte-Carlo Test for the Ergodicity Information Index,
#' the null hypothesis is that the empirical value of EII is equal to the Monte-Carlo value of EII obtained in multiple individuals with a similar latent structure.
#' Small values of p indicate that is very unlikely to obtain an EII as large as the one obtained in the empirical sample if the null hypothesis is true, thus there is convincing evidence that the empirical Ergodicity Information Index is
#' different than it could be expected if all individuals had a similar latent structure, conditioned on the parameters used to simulate the data.
#'
#' @param iter Numeric integer.
#' Number of random samples to generate in the Monte-Carlo simulation.
#' At least \code{500} is recommended
#'
#' @param N Numeric integer.
#' Number of individuals to simulate data from, using the \code{\link[EGAnet]{simDFM}} function.
#'
#' @param EII Numeric.
#' Empirical Ergodicity Information Index obtained via the \code{\link[EGAnet]{ergoInfo}} function.
#'
#' @param variab Number of variables per factor.
#'
#' @param timep Number of time points.
#'
#' @param nfact Number of factors.
#'
#' @param error Value to be used to construct a diagonal matrix Q. This matrix is p x p covariance matrix Q that will
#' generate random errors following a multivariate normal distribution with mean zeros.
#' The value provided is squared before constructing Q.
#'
#' @param dfm A string indicating the dynamical factor model to use. Defaults to \code{"DAFS"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{DAFS}}}
#' {Simulates data using the direct autoregressive factor score model.
#' This is the default method}
#'
#' \item{\strong{\code{RandomWalk}}}
#' {Simulates data using a dynamic factor model with random walk factor scores.}
#'}
#'
#' @param loadings Magnitude of the loadings.
#'
#' @param autoreg Magnitude of the autoregression coefficients.
#' Default is \code{"autoreg = 0.8"}.
#'
#' @param crossreg Magnitude of the cross-regression coefficients.
#' Default is \code{"crossreg = 0.1"}.
#'
#' @param var.shock Magnitude of the random shock variance.
#' Default is \code{"var.shock = 0.18"}.
#'
#' @param cov.shock Magnitude of the random shock covariance
#' Default is \code{"cov.shock = 0.36"}.
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
#'
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
#' @param embed Integer.
#' Number of embedded dimensions (the number of observations to be used in the \code{\link[EGAnet]{Embed}} function). For example,
#' an \code{"embed = 5"} will use five observations to estimate a single derivative. Defaults to \code{embed = 5}.
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
#' or \code{\link[NetworkToolbox]{TMFG}}
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
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @examples
#'
#' \dontrun{
#' \donttest{
#' dyn1 <- dynEGA.ind.pop(data = sim.dynEGA, n.embed = 5, tau = 1,
#'                       delta = 1, id = 21, group = 22, use.derivatives = 1,
#'                     model = "glasso", ncores = 2, corr = "pearson")
#'
#' eii1 <- ergoInfo(data = dyn1)$EII
#'
#' dist.ergoinfo <- mctest.ergoInfo(iter = 10, N = 10, EII = eii1,
#' variab = 4,
#' timep = 100, nfact = 2, error = 0.05, dfm = "DAFS", loadings = 0.55, autoreg = 0.8,
#' crossreg = 0.1, var.shock = 0.18, cov.shock = 0.36, embed = 5, tau=1, delta=1, derivatives=1,
#' model = "glasso", ncores = 2, corr = "pearson")
#' }}
#'
#' @return Returns a list containing:
#'
#' \item{mc.ergoInfo}{The values of the Ergodicity Information Index obtained in the Monte-Carlo Simulation}
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
# Updated 21.10.2020


mctest.ergoInfo <- function(iter, N,
                          EII,
                          variab,
                          timep,
                          nfact,
                          error,
                          dfm,
                          loadings,
                          autoreg,
                          crossreg,
                          var.shock,
                          cov.shock,
                          embed, tau, delta, derivatives,
                          model, model.args = list(),
                          algorithm = c("walktrap", "louvain"),
                          algorithm.args = list(),
                          corr, ncores, ...
){


  #### MISSING ARGUMENTS HANDLING ####

  if(missing(dfm))
  {dfm <- "DAFS"
  }else{dfm}

  if(missing(autoreg))
  {autoreg <- 0.8
  }else{autoreg}

  if(missing(crossreg))
  {crossreg <- 0.1
  }else{crossreg}

  if(missing(var.shock))
  {var.shock <- 0.18
  }else{var.shock}

  if(missing(cov.shock))
  {cov.shock <- 0.36
  }else{cov.shock}

  if(missing(embed))
  {embed <- 5
  }else{embed}

  if(missing(tau))
  {tau <- 5
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
  for(i in 1:iter){
    data.sim[[i]] <- vector("list", length = N)
  }



  for(i in 1:iter){
    for(j in 1:N){
      #generate data
      data.sim[[i]][[j]] <- EGAnet::simDFM(variab = variab, timep = timep, nfact = nfact, error = error,
                                   dfm = dfm, loadings = loadings, autoreg = autoreg,
                                   crossreg = crossreg, var.shock = var.shock,
                                   cov.shock = cov.shock, burnin = 1000)$data
    }
  }

  data.sim.df <- vector("list", length = iter)
  for(i in 1:iter){
    data.sim.df[[i]] <- purrr::map_df(data.sim[[i]], ~as.data.frame(.))
    data.sim.df[[i]]$ID <- rep(1:N, each = timep)
  }


  #initialize correlation matrix list
  sim.dynEGA <- vector("list", length = iter)
  sim.dynEGA.ids <- vector("list", length = iter)
  list.results.sim <- vector("list", length = iter)
  complexity.estimates <- vector("list", length = iter)

  #let user know data generation has started
  message("\nEstimating the Population and Individual Structures...\n", appendLF = FALSE)

  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("sim.dynEGA", "list.results.sim",
                                      "complexity.estimates"),
                          envir=environment())

  #Compute DynEGA in the population and in the individuals
  sim.dynEGA <- pbapply::pblapply(X = data.sim.df, cl = cl,
                                  FUN = EGAnet::dynEGA.ind.pop,
                                  n.embed = embed, tau = tau,
                                  delta = delta, id = (variab*nfact)+1, use.derivatives = derivatives,
                                  algorithm = algorithm, algorithm.args = algorithm.args,
                                  model = model, model.args = model.args, corr = corr)


  #let user know data generation has started
  message("Estimating the Ergodicity Information Index\n", appendLF = FALSE)


  complexity.estimates <- pbapply::pblapply(X = sim.dynEGA, cl = cl,
                                            FUN = EGAnet::ergoInfo)
  parallel::stopCluster(cl)


  complexity.estimates2 <- vector("list")
  for(i in 1:length(complexity.estimates)){
    complexity.estimates2[[i]] <- complexity.estimates[[i]]$EII
  }

  complexity.estimates2 <- unlist(complexity.estimates2)

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
  results$mc.ergoInfo <- complexity.estimates2
  results$p.value.twosided <- two.sided
  results$effect <- ifelse(p.greater<p.greater, "Greater", "Less")
  results$plot.dist <- plot.bootErgoInfo
  return(results)
}
#----
