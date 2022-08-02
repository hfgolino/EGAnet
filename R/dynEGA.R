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
#' @param corr Type of correlation matrix to compute. The default uses \code{"pearson"}.
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
#' A string indicating the method to use.
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
#' Defaults to \code{"walktrap"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#' 
#' \item{\strong{\code{leiden}}}
#' {Computes the Leiden algorithm using \code{\link[igraph]{cluster_leiden}}.
#' Defaults to \code{objective_function = "modularity"}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Louvain algorithm using \code{\link[igraph]{cluster_louvain}}}
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
#' \donttest{# Population structure
#' dyn.random <- dynEGA(
#'   data = sim.dynEGA, n.embed = 5, tau = 1,
#'   delta = 1, id = 21, group = 22, use.derivatives = 1,
#'   level = "population", ncores = 2, corr = "pearson"
#' )
#'
#' # Plot population structure
#' plot(dyn.random)
#'
#' # Group structure
#' dyn.group <- dynEGA(
#'   data = sim.dynEGA, n.embed = 5, tau = 1,
#'   delta = 1, id = 21, group = 22, use.derivatives = 1,
#'   level = "group", ncores = 2, corr = "pearson"
#' )
#'
#' # Plot group structure
#' plot(dyn.group, ncol = 2, nrow = 1)
#'
#' # Intraindividual structure
#' dyn.individual <- dynEGA(
#'   data = sim.dynEGA, n.embed = 5, tau = 1,
#'   delta = 1, id = 21, group = 22, use.derivatives = 1,
#'   level = "individual", ncores = 2, corr = "pearson"
#' )
#'
#' # Plot individual structure (participant 1)
#' plot(dyn.individual, id = 1)}
#'
#' @references
#' Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010)
#' Generalized local linear approximation of derivatives from time series. In S.-M. Chow, E. Ferrer, & F. Hsieh (Eds.),
#' \emph{The Notre Dame series on quantitative methodology. Statistical methods for modeling human dynamics: An interdisciplinary dialogue},
#' (p. 161-178). \emph{Routledge/Taylor & Francis Group}.
#'
#' Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
#' Using derivative estimates to describe intraindividual variability at multiple time scales.
#' \emph{Psychological Methods}, \emph{14(4)}, 367-386.
#'
#' Golino, H., Christensen, A. P., Moulder, R. G., Kim, S., & Boker, S. M. (2021).
#' Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
#' \emph{Psychometrika}.
#'
#' Savitzky, A., & Golay, M. J. (1964).
#' Smoothing and differentiation of data by simplified least squares procedures.
#' \emph{Analytical Chemistry}, \emph{36(8)}, 1627-1639.
#'
#' @importFrom stats cor rnorm runif na.omit
#'
#' @export
# dynEGA
# Updated 18.07.2022
dynEGA <- function(data, n.embed, tau = 1, delta = 1,
                   level = c("individual", "group", "population"),
                   id = NULL, group = NULL,
                   use.derivatives = 1,
                   model = c("glasso", "TMFG"), model.args = list(),
                   algorithm = c("walktrap", "leiden", "louvain"), algorithm.args = list(),
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

  # Check if cor has been input as an argument
  if("cor" %in% names(add.args)){

    # Give deprecation warning
    warning(
      paste(
        "The 'cor' argument has been deprecated in dynEGA.\n\nInstead use: corr = ", add.args$cor, ")",
        sep = ""
      )
    )

    # Handle the number of steps appropriately
    corr <- add.args$cor
  }

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

  if(missing(model)){
    model <- "glasso"
  }else{
    model <- match.arg(model)
  }

  if(missing(algorithm)){
    algorithm <- "walktrap"
  }else{
    algorithm <- match.arg(algorithm)
  }

  if(missing(corr))
  {corr <- "pearson"
  }else{corr <- match.arg(corr)}

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
  if("group" %in% tolower(colnames(data))){
    datalist <- split(data[,-c(id, group)], data[,id])
  }else{
    datalist <- split(data[,-c(id)], data[,id])
  }

  ### Estimating the derivatives using GLLA:

  # Let user know derivatives estimation has started
  message("\nComputing derivatives using GLLA...\n", appendLF = FALSE)

  #Parallel processing
  cl <- parallel::makeCluster(ncores)

  #Export variables
  # parallel::clusterExport(cl = cl,
  #                         varlist = c("datalist"),
  #                         envir=environment())

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

  # Compute derivatives per ID
  derivlist <- pbapply::pblapply(X = datalist, cl = cl,
                                 FUN = glla.multi,
                                 n.embed = n.embed, tau = tau, delta = delta, order = order)

  # Stop cluster
  parallel::stopCluster(cl)

  ### Estimating the dimensionality structure using EGA:

  message("Estimating the dimensionality structure using EGA...\n", appendLF = FALSE)

  for(i in 1:length(cases)){
    derivlist[[i]]$ID <- data[which(data[,id]==cases[i]),id][-c(1:(n.embed-1))]
  }

  names(derivlist) <- paste0("ID", cases)

  # Population Level:
  if(level == "population"){

    # Message user about level
    message("Level: Population...", appendLF = FALSE)

    # Stack derivatives
    data.all <- data.frame(Reduce(rbind, derivlist))

    # Obtain derivative indices
    derivative_index <- switch(
      as.character(use.derivatives),
      "0" = grep("Ord0", colnames(data.all)),
      "1" = grep("Ord1", colnames(data.all)),
      "2" = grep("Ord2", colnames(data.all))
    )

    # Estimate using EGA
    ega1 <- suppressWarnings(
      EGA(data = data.all[, derivative_index],
          model = model, model.args = model.args,
          algorithm = algorithm, algorithm.args = algorithm.args,
          corr = corr, plot.EGA = FALSE)
    )

    # Message user that results are done
    message("done", appendLF = TRUE)
  }

  # Level Group:
  if(level == "group"){

    # Message user about level
    message("Level: Group...", appendLF = FALSE)

    # Obtain group memberships
    group.memb <- unique(data[,group])

    # Assign derivatives to groups
    for(i in 1:length(cases)){
      derivlist[[i]]$Group <- data[which(data[,id]==cases[i]),group][-c(1:(n.embed-1))]
    }

    # Stack group memberships
    data.all <- data.frame(Reduce(rbind, derivlist))

    # Obtain derivative indices
    derivative_index <- switch(
      as.character(use.derivatives),
      "0" = grep("Ord0", colnames(data.all)),
      "1" = grep("Ord1", colnames(data.all)),
      "2" = grep("Ord2", colnames(data.all))
    )

    # Initialize data list
    data.groups <- vector("list", length = length(group.memb))
    data.groups <- split(data.all[,derivative_index], data.all$Group)
    names(data.groups) <- paste0("Group", group.memb)

    # Parallel processing
    cl <- parallel::makeCluster(ncores)

    # Export variables
    # parallel::clusterExport(cl = cl,
    #                         varlist = c("data.groups", "group.memb"),
    #                         envir=environment())

    # Compute derivatives per Group
    ega.list.groups <- pbapply::pblapply(X = data.groups, cl = cl,
                                         FUN = EGA,
                                         model = model, model.args = model.args,
                                         algorithm = algorithm, algorithm.args = algorithm.args,
                                         corr = corr, plot.EGA = FALSE)

    # Stop cluster
    parallel::stopCluster(cl)

  }

  # Level: Individual (intraindividual structure):
  if(level == "individual"){

    # Message user about level
    message("Level: Individual (Intraindividual Structure)...\n", appendLF = FALSE)

    # Stack group memberships
    data.all <- data.frame(Reduce(rbind, derivlist))

    # Obtain derivative indices
    derivative_index <- switch(
      as.character(use.derivatives),
      "0" = grep("Ord0", colnames(data.all)),
      "1" = grep("Ord1", colnames(data.all)),
      "2" = grep("Ord2", colnames(data.all))
    )

    # Initialize data list
    data.individuals <- vector("list", length = length(cases))
    data.individuals <- split(data.all[,derivative_index], data.all$ID)
    names(data.individuals) <- paste0("ID", cases)

    # Get number of variables
    initial.nvar <- unlist(lapply(data.individuals, ncol))

    # Remove variables from participants with no variance
    # in their derivatives
    data.individuals_var <- lapply(data.individuals, function(x){
      indices <- which(apply(x, 2, sd) == 0)
      if(length(indices) != 0){
        x[,-indices]
      }else{x}
    })
    # ^^ creates new object to keep reduced columns
    # separate to avoid errors when creating `dim.variables` later

    # Get number of variables
    final.nvar <- unlist(lapply(data.individuals_var, ncol))

    # Get warnings
    warning.idx <- which(initial.nvar != final.nvar)
    if(length(warning.idx) != 0){

      for(i in 1:length(warning.idx)){
        warning(
          paste(
            names(warning.idx)[i],
            "had variables with no variance. Some variables will be disconnected in their network."
          )
        )
      }

    }

    #Parallel processing
    cl <- parallel::makeCluster(ncores)

    #Export variables
    # parallel::clusterExport(cl = cl,
    #                         varlist = c("data.individuals_var"),
    #                         envir=environment())

    # EGA estimates per individual:
    # op <- pbapply::pboptions(type = "none")
    ega.list.individuals <- pbapply::pblapply(
      X = data.individuals_var, cl = cl,
      FUN = EGA,
      model = model, model.args = model.args,
      algorithm = algorithm, algorithm.args = algorithm.args,
      corr = corr, plot.EGA = FALSE
    )
    # pbapply::pboptions(op)

    # Stop cluster
    parallel::stopCluster(cl)

    # Add back disconnected variables
    if(length(warning.idx) != 0){
      for(i in 1:length(warning.idx)){

        # Obtain target individual
        target_individual <- ega.list.individuals[[names(warning.idx)[i]]]

        # Obtain network
        network <- target_individual$network

        # Obtain columns to use
        colstouse <- colnames(data.all)[derivative_index]

        # Identify zero variance variables
        zero_var <- setdiff(colstouse, colnames(network))

        # Add zero variance variables to network
        new_network <- matrix(
          0,
          ncol = length(colstouse),
          nrow = length(colstouse)
        )
        colnames(new_network) <- c(colnames(network), zero_var)
        row.names(new_network) <- colnames(new_network)

        # Insert existing variables into network
        new_network[1:nrow(network), 1:ncol(network)] <- network

        # Update wc
        wc <- target_individual$wc
        add_wc <- rep(NA, length(zero_var))
        names(add_wc) <- zero_var
        new_wc <- c(wc, add_wc)

        # Update correlation matrix
        correlation <- target_individual$correlation

        # Add zero variance variables to correlation
        new_correlation <- matrix(
          NA,
          ncol = length(colstouse),
          nrow = length(colstouse)
        )
        colnames(new_correlation) <- c(colnames(correlation), zero_var)
        row.names(new_correlation) <- colnames(new_correlation)

        # Insert existing variables into correlation
        new_correlation[1:nrow(correlation), 1:ncol(correlation)] <- correlation

        # Update order
        new_network <- new_network[colstouse, colstouse]
        new_wc <- new_wc[colstouse]
        new_correlation <- new_correlation[colstouse, colstouse]

        # Put back into object
        target_individual$network <- new_network
        target_individual$wc <- new_wc
        target_individual$correlation <- new_correlation
        ega.list.individuals[[names(warning.idx)[i]]] <- target_individual

      }
    }

  }

  # Results:
  results <- vector("list")
  results$Derivatives <- vector("list")
  results$Derivatives$Estimates <- derivlist
  results$Derivatives$EstimatesDF <- data.all


  # Set up methods
  methods <- list(); methods$glla <- list(); methods$EGA <- list()

  # GLLA methods
  methods$glla$n.embed <- n.embed
  methods$glla$tau <- tau
  methods$glla$delta <- delta
  methods$glla$derivatives <- use.derivatives

  # EGA methods
  methods$EGA$model <- model
  methods$EGA$model.args <- model.args
  methods$EGA$algorithm <- algorithm
  methods$EGA$algorithm.args <- algorithm.args
  methods$EGA$corr <- corr

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
      ega.list.individuals <- ega.list.individuals[!unlist(
        lapply(ega.list.individuals, is.null)
      )]
      results$dynEGA <- ega.list.individuals
      class(results) <- "dynEGA.Individuals"
      dim.variables <- list()
      for(i in 1:length(ega.list.individuals)){
        dim.variables[[i]] <- data.frame(items = colnames(data.individuals[[i]]), dimension = ega.list.individuals[[i]]$wc)
        dim.variables[[i]] <- dim.variables[[i]][order(dim.variables[[i]][, 2]),]
        results$dynEGA[[i]]$dim.variables <- dim.variables[[i]]
      }
  }

  # Attach methods to dynEGA output
  results$dynEGA$Methods <- methods

  return(results)
}
#----
