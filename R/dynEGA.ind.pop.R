#' Dynamic EGA
#' 
#' @description DynEGA estimates dynamic factors in multivariate time series (i.e. longitudinal data, panel data, intensive longitudinal data) at multiple
#' time scales, in different levels of analysis: individuals (intraindividual structure) and population (structure of the population).
#' Exploratory graph analysis is applied in the derivatives estimated using generalized local linear approximation (\code{\link[EGAnet]{glla}}). Instead of estimating factors by modeling how variables are covarying, as in traditional
#' EGA, dynEGA is a dynamic model that estimates the factor structure by modeling how variables are changing together.
#' GLLA is a filtering method for estimating derivatives from data that uses time delay embedding and a variant of Savitzky-Golay filtering to accomplish the task.
#'
#' @param data A data frame with the variables to be used in the analysis.
#' The data frame should be in a long format (i.e. observations for the
#' same individual (for example, individual 1) are placed in order,
#' from time 1 to time t, followed by the observations from individual 2, also 
#' ordered from time 1 to time t.)
#'
#' @param n.embed Integer.
#' Number of embedded dimensions (the number of observations to be used in the \code{\link[EGAnet]{Embed}} function). For example,
#' an \code{n.embed = 5} will use five consecutive observations to estimate a single derivative.
#'
#' @param tau Integer.
#' Number of observations to offset successive embeddings in the \code{\link[EGAnet]{Embed}} function. A tau of one uses adjacent observations.
#' Default is \code{tau = 1}.
#'
#' @param delta Integer.
#' The time between successive observations in the time series.
#' Default is \code{delta = 1}.
#'
#' @param id Numeric.
#' Number of the column identifying each individual.
#'
#'
#' @param use.derivatives Integer.
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
#' @param ... Additional arguments.
#' Used for deprecated arguments from previous versions of \code{\link{EGA}}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @examples
#' # Obtain data
#' sim.dynEGA <- sim.dynEGA # bypasses CRAN checks
#' 
#' \donttest{# Dynamic EGA individual and population structure
#' dyn.ega1 <- dynEGA.ind.pop(
#'   data = sim.dynEGA, n.embed = 5, tau = 1,
#'   delta = 1, id = 21, use.derivatives = 1, 
#'   ncores = 2, corr = "pearson"
#' )}
#' 
#' @export
#'
# Updated 14.06.2022
dynEGA.ind.pop <- function(data, n.embed, tau = 1, delta = 1,
                           id = NULL,
                           use.derivatives = 1,
                           model = c("glasso", "TMFG"), model.args = list(),
                           algorithm = c("walktrap", "louvain"), algorithm.args = list(),
                           corr = c("cor_auto", "pearson", "spearman"),
                           ncores, ...){
  
  # Get additional arguments
  add.args <- list(...)
  
  # Check if steps has been input as an argument
  if("steps" %in% names(add.args)){
    
    # Give deprecation warning
    warning(
      paste(
        "The 'steps' argument has been deprecated in all EGA functions.\n\nInstead use: algorithm.args = list(steps = ", add.args$steps, ")",
        sep = ""
      )
    )
    
    # Handle the number of steps appropriately
    algorithm.args$steps <- add.args$steps
  }
  
  
  # MISSING ARGUMENTS HANDLING
  if(missing(id))
  {stop("The 'id' argument is missing! \n The number of the column identifying each individual must be provided!")
  }else{id <- id}
  
  if(missing(corr))
  {corr <- "cor_auto"
  }else{corr <- match.arg(corr)}
  
  if(missing(ncores))
  {ncores <- ceiling(parallel::detectCores() / 2)
  }else{ncores}
  
  # # Setting the order:
  # 
  # order = 2
  # 
  # # ### Spliting by ID:
  # 
  # #number of cases
  # cases <- unique(data[,id])
  # 
  # #initialize data list
  # datalist <- vector("list", length = length(cases))
  # datalist <- split(data[,-c(id)],data[,id])
  
  # Remove group
  if("group" %in% tolower(colnames(data))){
    datalist <- data[,-which(tolower(colnames(data)) == "group")]
  }else{
    datalist <- data
  }
  
  # Perform dynEGA for individuals first
  ega_individuals <- dynEGA(
    data = datalist, n.embed = n.embed,
    tau = tau, delta = delta, level = "individual",
    id = id, use.derivatives = use.derivatives,
    model = model, model.args = model.args,
    algorithm = algorithm, algorithm.args = algorithm.args,
    corr = corr, ncores = ncores
  )
  
  # Stack derivatives
  derivatives_df <- ega_individuals$Derivatives$EstimatesDF
  
  # Obtain proper derivatives
  derivatives_df <- derivatives_df[,grep(
    paste("Ord", use.derivatives, sep = ""),
    colnames(derivatives_df)
  )]
  
  # Message user for population structure
  message("Level: Population...", appendLF = FALSE)
  
  # Estimate population structure
  ega_pop <- suppressWarnings(
    EGA(
      derivatives_df,
      model = model, model.args = model.args,
      algorithm = algorithm, algorithm.args = algorithm.args,
      corr = corr, plot.EGA = FALSE
    )
  )
  
  # Set up EGA population object
  ega_population <- list()
  ega_population$Derivatives <- ega_individuals$Derivatives
  ega_population$dynEGA <- ega_pop
  ega_population$dynEGA$Methods <- ega_individuals$dynEGA$Methods
  
  # Change class
  class(ega_population) <- "dynEGA"
  
  # ega_population <- dynEGA(
  #   data = datalist, n.embed = n.embed,
  #   tau = tau, delta = delta, level = "population",
  #   id = id, use.derivatives = use.derivatives,
  #   model = model, model.args = model.args,
  #   algorithm = algorithm, algorithm.args = algorithm.args,
  #   corr = corr, ncores = ncores
  # )
  
  # # Derivatives list
  # derivlist <- ega_population$Derivatives$Estimates
  # 
  # # Obtain `use.derivatives` to get individual data
  # data.individuals <- lapply(derivlist, function(x){
  #   x[,grep(
  #     paste("Ord", use.derivatives, sep = ""),
  #     colnames(x)
  #   )]
  # })
  # 
  # # ### Estimating the derivatives using GLLA:
  # # 
  # # #let user know derivatives estimation has started
  # # message("\nComputing derivatives using GLLA...\n", appendLF = FALSE)
  # # 
  # # #initialize derivatives list
  # # derivlist <- list()
  # # 
  # # #Parallel processing
  # # cl <- parallel::makeCluster(ncores)
  # # 
  # # #Export variables
  # # parallel::clusterExport(cl = cl,
  # #                         varlist = c("datalist", "derivlist", "cases"),
  # #                         envir=environment())
  # # 
  # # # GLLA Estimation:
  # # glla.multi <- function(data, n.embed = n.embed, tau = tau, delta = delta, order = order){
  # #   order.deriv <- paste0("Ord",seq(from = 0, to = order))
  # #   data.est <- vector("list")
  # #   for(i in 1:ncol(data)){
  # #     data.est[[i]] <- as.data.frame(EGAnet::glla(data[,i], n.embed = n.embed, tau = tau, delta = delta, order = order))
  # #     data.est[[i]] <- as.data.frame(data.est[[i]])
  # #   }
  # #   data.est2 <- vector("list")
  # #   for(i in 0:order+1){
  # #     data.est2[[i]] <- sapply(data.est, "[[", i)
  # #   }
  # #   
  # #   data.estimates <- data.frame(Reduce(cbind, data.est2))
  # #   colnames(data.estimates) <- paste(colnames(data), rep(order.deriv, each = ncol(data)), sep = ".")
  # #   return(data.estimates)
  # # }
  # # 
  # # #Compute derivatives per ID
  # # derivlist <- pbapply::pblapply(X = datalist, cl = cl,
  # #                                FUN = glla.multi,
  # #                                n.embed = n.embed, tau = tau, delta = delta, order = order)
  # # 
  # # ### Estimating the dimensionality structure using EGA:
  # # 
  # # message("Estimating the dimensionality structure using EGA...\n", appendLF = FALSE)
  # # 
  # # for(i in 1:length(cases)){
  # #   derivlist[[i]]$ID <- data[which(data[,id]==cases[i]),id][-c(1:(n.embed-1))]
  # # }
  # # 
  # # names(derivlist) <- paste0("ID", cases)
  # # # Population Level:
  # # message("Level: Population...\n", appendLF = FALSE)
  # # 
  # # data.all <- data.frame(Reduce(rbind, derivlist))
  # # 
  # # # EGA Part
  # # 
  # # if(use.derivatives == 0){
  # #   ega1 <- EGA.estimate(data = data.all[,1:ncol(data[,-c(id)])],
  # #                        model = model, model.args = model.args,
  # #                        algorithm = algorithm, algorithm.args = algorithm.args,
  # #                        corr = corr)}
  # # if(use.derivatives == 1){
  # #   ega1 <- EGA.estimate(data = data.all[,(ncol(data[,-c(id)])+1):(ncol(data[,-c(id)])*2)],
  # #                        model = model, model.args = model.args,
  # #                        algorithm = algorithm, algorithm.args = algorithm.args,
  # #                        corr = corr)}
  # # if(use.derivatives==2){
  # #   init <- (ncol(data[,-c(id)])*2)+1
  # #   cols <- seq(from = init, to = init+ncol(data[,-c(id)])-1)
  # #   ega1 <- EGA.estimate(data = data.all[,cols],
  # #                        model = model, model.args = model.args,
  # #                        algorithm = algorithm, algorithm.args = algorithm.args,
  # #                        corr = corr)}
  # # 
  # # parallel::stopCluster(cl)
  # # 
  # # # Level: Individual (intraindividual structure):
  # # message("Level: Individual (Intraindividual Structure)...\n", appendLF = FALSE)
  # # 
  # # data.all <- data.frame(Reduce(rbind, derivlist))
  # # 
  # # # Which derivatives to use:
  # # if(use.derivatives == 0){
  # #   colstouse <- colnames(data.all[,1:ncol(data[,-c(id)])])}
  # # if(use.derivatives == 1){
  # #   colstouse <- colnames(data.all[,(ncol(data[,-c(id)])+1):(ncol(data[,-c(id)])*2)])}
  # # if(use.derivatives==2){
  # #   init <- (ncol(data[,-c(id)])*2)+1
  # #   cols <- seq(from = init, to = init+ncol(data[,-c(id)])-1)
  # #   colstouse <- colnames(data.all[,cols])
  # # }
  # # 
  # # #initialize data list
  # # data.individuals <- vector("list", length = length(cases))
  # # data.individuals <- split(data.all[,colstouse],data.all$ID)
  # # names(data.individuals) <- paste0("ID", cases)
  # 
  # message("Level: Individual (Intraindividual Structure)...", appendLF = FALSE)
  # 
  # #Parallel processing
  # cl <- parallel::makeCluster(ncores)
  # 
  # #Export variables
  # # parallel::clusterExport(cl = cl,
  # #                         varlist = c("data.individuals", "cases"),
  # #                         envir=environment())
  # 
  # # EGA estimates per individual:
  # ega.list.individuals <- list()
  # 
  # ega.list.individuals <- pbapply::pblapply(X = data.individuals, cl = cl,
  #                                           FUN = EGA.estimate,
  #                                           model = model, model.args = model.args,
  #                                           algorithm = algorithm, algorithm.args = algorithm.args,
  #                                           corr = corr)
  # parallel::stopCluster(cl)
  # 
  # #let user know results have been computed
  # message("done", appendLF = TRUE)
  
  # Results:
  results <- vector("list")
  results$Derivatives <- ega_individuals$Derivatives # vector("list")
  # results$Derivatives$Estimates <- derivlist
  # results$Derivatives$EstimatesDF <- data.all
  results$dynEGA.pop <- ega_population # ega1
  results$dynEGA.ind <- ega_individuals # ega.list.individuals
  # if(use.derivatives == 0){
  #   results$data.all <- data.all[,1:ncol(data[,-c(id)])]}
  # if(use.derivatives == 1){
  #   results$data.all <- data.all[,(ncol(data[,-c(id)])+1):(ncol(data[,-c(id)])*2)]}
  # if(use.derivatives == 2){
  #   results$data.all <- data.all[,cols]}
  results$data.all <- derivatives_df
  results$data.individuals <- ega_individuals$Derivatives$Estimates # data.individuals
  results$model <- model
  class(results) <- "dynEGA.ind.pop"
  return(results)
}
