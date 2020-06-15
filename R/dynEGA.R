#' Dynamic Exploratory Graph Analysis
#'
#' @description Estimates dynamic factors in multivariate time series (i.e. longitudinal data, panel data, intensive longitudinal data) at multiple
#' time scales, in different levels of analysis: individuals (intraindividual structure), groups or population (structure of the population).
#' Exploratory graph analysis is applied in the derivatives estimated using generalized local linear approximation (\code{\link[EGAnet]{glla}}). Instead of estimating factors by modeling how variables are covarying, as in traditional
#' EGA, dynEGA is a dynamic model that estimates the factor structure by modeling how variables are changing together.
#' GLLA is a filtering method for estimating derivatives from data that uses time delay embedding and a variant of Savitzky-Golay filtering to accomplish the task.
#'
#' @param data A dataframe with the variables to be used in the analysis. The dataframe should be in a long format (i.e. observations for the same individual (for example, individual 1) are placed in order, from time 1 to time t, followed by the observations from individual 2, also ordered from time 1 to time t.)
#'
#' @param n.embed Integer.
#' Number of embedded dimensions (the number of observations to be used in the \code{\link[EGAnet]{Embed}} function). For example,
#' an \code{"n.embed = 5"} will use five consecutive observations to estimate a single derivative.
#'
#' @param tau Integer.
#' Number of observations to offset successive embeddings in the \code{\link[EGAnet]{Embed}} function. A tau of one uses adjacent observations.
#' Default is \code{"tau = 1"}.
#'
#' @param delta Integer.
#' The time between successive observations in the time series.
#' Default is \code{"delta = 1"}.
#'
#' @param level Character.
#' A string indicating the level of analysis. If the interest is
#' in modeling the intraindividual structure only (one dimensionality structure per individual), then \code{level} should be set to \code{"individual"}.
#' If the interest is in the structure of a group of individuals, then \code{level} should be set to \code{"group"}.
#' Finally, if the interest is in the population structure, then \code{level} should be set to \code{"population"}.
#'
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{individual}}}
#' {Estimates the dynamic factors per individual. This should be the prefered method is one is interested in
#' in the factor structure of individuals. An additional parameter (\code{"id"}) needs to be provided identifying each individual.}
#'
#' \item{\strong{\code{group}}}
#' {Estimates the dynamic factors for each group.
#' An additional parameter (\code{"group"}) needs to be provided identifying the group membership.}
#'
#' \item{\strong{\code{population}}}
#' {Estimates the dynamic factors of the population}
#'
#'}
#' @param id Numeric.
#' Number of the column identifying each individual.
#'
#' @param group Numeric or character.
#' Number of the column identifying group membership. Must be specified only if \code{level = "group"}.
#'
#' @param use.derivatives Integer.
#' The order of the derivative to be used in the EGA procedure. Default to 1.
#'
#' @param plot.EGA Logical.
#' If TRUE, returns a plot of the network and its estimated dimensions.
#' Defaults to TRUE
#'
#' @param cor Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
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
#' @param model A string indicating the network method to use (\code{\link{EGA.estimate}}).
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
#' @param algorithm A string indicating the community detection algorithm to use.
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
#' @param steps Number of steps to be used in \code{\link[igraph]{cluster_walktrap}} algorithm.
#' Defaults to 4.
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
#' @author Hudson Golino <hfg9s at virginia.edu>
#'
#' @examples
#'
#' \donttest{
#' # Population structure:
#'dyn.random <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#'delta = 1, id = 21, group = 22, use.derivatives = 1,
#'level = "population", model = "glasso", ncores = 2)
#'
#'
#'# Group structure:
#'dyn.group <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#'delta = 1, id = 21, group = 22, use.derivatives = 1,
#'level = "group", model = "glasso", ncores = 2)
#'
#'# Intraindividual structure:
#'
#'dyn.individual <- dynEGA(data = sim.dynEGA, n.embed = 5, tau = 1,
#'delta = 1, id = 21, group = 22, use.derivatives = 1,
#'level = "individual", model = "glasso", ncores = 2)
#'}
#' @references
#'
#' Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010)
#' Generalized local linear approximation of derivatives from time series. In S.-M. Chow, E. Ferrer, & F. Hsieh (Eds.),
#' \emph{The Notre Dame series on quantitative methodology. Statistical methods for modeling human dynamics: An interdisciplinary dialogue},
#' (p. 161-178). \emph{Routledge/Taylor & Francis Group}.
#' doi:\href{https://doi.org/10.1037/a0016622}{10.1037/a0016622}
#'
#' Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
#' Using derivative estimates to describe intraindividual variability at multiple time scales.
#' \emph{Psychological Methods}, \emph{14(4)}, 367-386.
#' doi:\href{https://doi.org/10.1037/a0016622}{10.1037/a0016622}
#'
#' Golino, H. F., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PloS one}, \emph{12(6)}, e0174035..
#' doi: \href{https://doi.org/10.1371/journal.pone.0174035}{journal.pone.0174035}

#' Savitzky, A., & Golay, M. J. (1964).
#' Smoothing and differentiation of data by simplified least squares procedures.
#' \emph{Analytical Chemistry}, \emph{36(8)}, 1627-1639.
#' doi:\href{https://doi.org/10.1021/ac60214a047}{10.1021/ac60214a047}
#'
#'
#' @importFrom stats cor rnorm runif na.omit
#'
#' @export
# dynEGA
# Updated 15.06.2020
#'
dynEGA <- function(data, n.embed, tau = 1, delta = 1,
                   level = c("individual", "group", "population"),
                   id = NULL, group = NULL,
                   use.derivatives = 1,
                   model = c("glasso", "TMFG"),
                   algorithm = c("walktrap", "louvain"),
                   plot.EGA = TRUE,
                   cor = c("cor_auto", "pearson", "spearman"),
                   steps = 4,
                   ncores){
  #### MISSING ARGUMENTS HANDLING ####
  if(missing(id))
  {stop("The 'id' argument is missing! \n The number of the column identifying each individual must be provided!")
  }else{id <- id}

  if(missing(level))
  {level <- "population"
  }else if(level == "group"){
    if(missing(group))
    {stop("Level set to 'group', but the 'group' argument is missing! \n The number of the column identifying each group must be provided!")
    }else{group <- group}
  }else{level <- match.arg(level)}

  if(missing(group))
  {group <- ncol(data)+1
  }else{group <- group}

  if(missing(cor))
  {cor <- "cor_auto"
  }else{cor <- match.arg(cor)}

  if(missing(ncores))
  {ncores <- ceiling(parallel::detectCores() / 2)
  }else{ncores}

  # Setting the order:

  order = 2

  ### Spliting by ID:

  #number of cases
  cases <- unique(data[,id])

  #initialize data list
  datalist <- vector("list", length = length(cases))
  datalist <- split(data[,-c(id, group)],data[,id])

  ### Estimating the derivatives using GLLA:

  #let user know derivatives estimation has started
  message("\nComputing derivatives using GLLA...\n", appendLF = FALSE)


  #initialize derivatives list
  derivlist <- list()

  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  parallel::clusterExport(cl = cl,
                          varlist = c("datalist", "derivlist", "cases"),
                          envir=environment())

  # GLLA Estimation:
  glla.multi <- function(data, n.embed = n.embed, tau = tau, delta = delta, order = order){
    order.deriv <- paste0("Ord",seq(from = 0, to = order))
    data.est <- vector("list")
    for(i in 1:ncol(data)){
      data.est[[i]] <- as.data.frame(EGAnet::glla(data[,i], n.embed = n.embed, tau = tau, delta = delta, order = order))
      data.est[[i]] <- as.data.frame(data.est[[i]])
    }
    data.est2 <- vector("list")
    for(i in 0:order+1){
      data.est2[[i]] <- sapply(data.est, "[[", i)
    }

    data.estimates <- data.frame(Reduce(cbind, data.est2))
    colnames(data.estimates) <- paste(colnames(data), rep(order.deriv, each = ncol(data)), sep = ".")
    return(data.estimates)
  }

  #Compute derivatives per ID
  derivlist <- pbapply::pblapply(X = datalist, cl = cl,
                                 FUN = glla.multi,
                                 n.embed = n.embed, tau = tau, delta = delta, order = order)

  ### Estimating the dimensionality structure using EGA:

  message("Estimating the dimensionality structure using EGA...\n", appendLF = FALSE)

  for(i in 1:length(cases)){
    derivlist[[i]]$ID <- data[which(data[,id]==cases[i]),id][-c(1:(n.embed-1))]
  }

  names(derivlist) <- paste0("ID", cases)
  # Population Level:
  if(level == "population"){
    message("Level: Population...\n", appendLF = FALSE)

    data.all <- data.frame(Reduce(rbind, derivlist))

    # EGA Part

    if(use.derivatives == 0){
      ega1 <- EGA.estimate(data = data.all[,1:ncol(data[,-c(id, group)])],
                           model = model, algorithm = algorithm,
                           steps = steps, cor = cor)}
    if(use.derivatives == 1){
      ega1 <- EGA.estimate(data = data.all[,(ncol(data[,-c(id, group)])+1):(ncol(data[,-c(id, group)])*2)],
                           model = model, algorithm = algorithm,
                           steps = steps, cor = cor)}
    if(use.derivatives==2){
      init <- (ncol(data[,-c(id, group)])*2)+1
      cols <- seq(from = init, to = init+ncol(data[,-c(id, group)])-1)
      ega1 <- EGA.estimate(data = data.all[,cols], model = model, algorithm = algorithm,
                           steps = steps, cor = cor)}
  }

  parallel::stopCluster(cl)

  # Level Group:
  if(level == "group"){
    message("Level: Group...\n", appendLF = FALSE)

    group.memb <- unique(data[,group])

    for(i in 1:length(cases)){
      derivlist[[i]]$Group <- data[which(data[,id]==cases[i]),group][-c(1:(n.embed-1))]
    }

    data.all <- data.frame(Reduce(rbind, derivlist))

    # Which derivatives to use:
    if(use.derivatives == 0){
      colstouse <- colnames(data.all[,1:ncol(data[,-c(id, group)])])}
    if(use.derivatives == 1){
      colstouse <- colnames(data.all[,(ncol(data[,-c(id, group)])+1):(ncol(data[,-c(id, group)])*2)])}
    if(use.derivatives==2){
      init <- (ncol(data[,-c(id, group)])*2)+1
      cols <- seq(from = init, to = init+ncol(data[,-c(id, group)])-1)
      colstouse <- colnames(data.all[,cols])
    }

    #initialize data list
    data.groups <- vector("list", length = length(group.memb))
    data.groups <- split(data.all[,colstouse],data.all$Group)
    names(data.groups) <- paste0("Group", group.memb)

    #Parallel processing
    cl <- parallel::makeCluster(ncores)

    #Export variables
    parallel::clusterExport(cl = cl,
                            varlist = c("data.groups", "group.memb"),
                            envir=environment())

    # EGA Part

    ega.list.groups <- list()

    #Compute derivatives per Group
    ega.list.groups <- pbapply::pblapply(X = data.groups, cl = cl,
                                         FUN = EGA.estimate,
                                         model = model, algorithm = algorithm,
                                         steps = steps, cor = cor)
    parallel::stopCluster(cl)
  }

  # Level: Individual (intraindividual structure):
  if(level == "individual"){
    message("Level: Individual (Intraindividual Structure)...\n", appendLF = FALSE)

    data.all <- data.frame(Reduce(rbind, derivlist))

    # Which derivatives to use:
    if(use.derivatives == 0){
      colstouse <- colnames(data.all[,1:ncol(data[,-c(id, group)])])}
    if(use.derivatives == 1){
      colstouse <- colnames(data.all[,(ncol(data[,-c(id, group)])+1):(ncol(data[,-c(id, group)])*2)])}
    if(use.derivatives==2){
      init <- (ncol(data[,-c(id, group)])*2)+1
      cols <- seq(from = init, to = init+ncol(data[,-c(id, group)])-1)
      colstouse <- colnames(data.all[,cols])
    }

    #initialize data list
    data.individuals <- vector("list", length = length(cases))
    data.individuals <- split(data.all[,colstouse],data.all$ID)
    names(data.individuals) <- paste0("ID", cases)


    #Parallel processing
    cl <- parallel::makeCluster(ncores)

    #Export variables
    parallel::clusterExport(cl = cl,
                            varlist = c("data.individuals", "cases"),
                            envir=environment())

    # EGA estimates per individual:
    ega.list.individuals <- list()

    ega.list.individuals <- pbapply::pblapply(X = data.individuals, cl = cl,
                                              FUN = EGA.estimate,
                                              model = model, algorithm = algorithm,
                                              steps = steps, cor = cor)
    parallel::stopCluster(cl)
  }

  #let user know results have been computed
  message("done", appendLF = TRUE)

  # Results:
  results <- vector("list")
  results$Derivatives <- vector("list")
  results$Derivatives$Estimates <- derivlist
  results$Derivatives$EstimatesDF <- data.all
  if(level == "population"){
    results$dynEGA <- ega1
    dim.variables <- data.frame(items = colnames(data[-c(id, group)]), dimension = ega1$wc)
    dim.variables <- dim.variables[order(dim.variables[, 2]),]
    results$dynEGA$dim.variables <- dim.variables
    class(results) <- "dynEGA"
  }else if(level == "group"){
    results$dynEGA <- ega.list.groups
    results$Derivatives$Estimates.Groups <- data.groups
    class(results) <- "dynEGA.Groups"
    dim.variables <- list()
    for(i in 1:length(group.memb)){
    dim.variables[[i]] <- data.frame(items = colnames(data[-c(id, group)]), dimension = ega.list.groups[[i]]$wc)
    dim.variables[[i]] <- dim.variables[[i]][order(dim.variables[[i]][, 2]),]
    results$dynEGA[[i]]$dim.variables <- dim.variables[[i]]}
    }else if(level == "individual"){
      results$dynEGA <- ega.list.individuals
      class(results) <- "dynEGA.Individuals"
      dim.variables <- list()
      for(i in 1:length(cases)){
        dim.variables[[i]] <- data.frame(items = colnames(data[-c(id, group)]), dimension = ega.list.individuals[[i]]$wc)
        dim.variables[[i]] <- dim.variables[[i]][order(dim.variables[[i]][, 2]),]
        results$dynEGA[[i]]$dim.variables <- dim.variables[[i]]}
    }
    return(results)
}
#----
