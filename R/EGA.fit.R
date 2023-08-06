#' @title \code{\link[EGAnet]{EGA}} Optimal Model Fit using the Total Entropy Fit Index  (\code{\link[EGAnet]{tefi}})
#'
#' @description Estimates the best fitting model using \code{\link[EGAnet]{EGA}}.
#' The number of steps in the \code{\link[igraph]{cluster_walktrap}} detection
#' algorithm is varied and unique community solutions are compared using
#' \code{\link[EGAnet]{tefi}}.
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param n Numeric (length = 1).
#' Sample size if \code{data} is a correlation matrix
#'
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"auto"} --- }
#' {Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)}
#' 
#' \item{\code{"pearson"} --- }
#' {Pearson's correlation is computed for all variables regardless of
#' categories}
#' 
#' \item{\code{"spearman"} --- }
#' {Spearman's rank-order correlation is computed for all variables
#' regardless of categories}
#' 
#' }
#' 
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#' 
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"} --- }
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"} --- }
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param model Character (length = 1).
#' Defaults to \code{"glasso"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"BGGM"} --- }
#' {Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{\link[BGGM]{estimate}} for more details}
#' 
#' \item{\code{"glasso"} --- }
#' {Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details}
#' 
#' \item{\code{"TMFG"} --- }
#' {Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details}
#' 
#' }
#' 
#' @param algorithm Character or \code{\link{igraph}} \code{cluster_*} function.
#' Three options are listed below but all are available
#' (see \code{\link[EGAnet]{community.detection}} for other options):
#' 
#' \itemize{
#'
#' \item{\code{"leiden"} --- }
#' {See \code{\link[igraph]{cluster_leiden}} for more details.
#' \emph{Note}: The Leiden algorithm will default to the
#' Constant Potts Model objective function
#' (\code{objective_function = "CPM"}). Set
#' \code{objective_function = "modularity"} to use
#' modularity instead (see examples). By default, searches along
#' resolutions from 0 to \code{max(abs(network))} or the maximum
#' absolute edge weight in the network in 0.01 increments
#' (\code{resolution_parameter = seq.int(0, max(abs(network)), 0.01)}). 
#' For modularity, searches along resolutions from 0 to 2 in 0.05 increments
#' (\code{resolution_parameter = seq.int(0, 2, 0.05)}) by default.
#' Use the argument \code{resolution_parameter} to change the search parameters
#' (see examples)}
#' 
#' \item{\code{"louvain"} --- }
#' {See \code{\link[EGAnet]{community.consensus}} for more details.
#' By default, searches along resolutions from 0 to 2 in 0.05 increments
#' (\code{resolution_parameter = seq.int(0, 2, 0.05)}). Use the argument \code{resolution_parameter}
#' to change the search parameters (see examples)}
#' 
#' \item{\code{"walktrap"} --- }
#' {This algorithm is the default. See \code{\link[igraph]{cluster_walktrap}} for more details.
#' By default, searches along 3 to 8 steps (\code{steps = 3:8}). Use the argument \code{steps}
#' to change the search parameters (see examples)}
#' 
#' }
#' 
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
#' 
#' @param verbose Boolean.
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}}, \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}}, \code{\link[EGAnet]{community.consensus}}, 
#' and \code{\link[EGAnet]{EGA.estimate}}
#'
#' @return Returns a list containing:
#' 
#' \item{EGA}{\code{\link[EGAnet]{EGA}} results of the best fitting solution}
#' 
#' \item{EntropyFit}{\code{\link[EGAnet]{tefi}} fit values for each solution}
#' 
#' \item{Lowest.EntropyFit}{The best fitting solution based on \code{\link[EGAnet]{tefi}}}
#' 
#' \item{parameter.space}{Parameter values used in search space}
#' 
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Estimate optimal EGA with Walktrap
#' fit.walktrap <- EGA.fit(
#'   data = wmt, algorithm = "walktrap",
#'   steps = 3:8, # default
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#' 
#' # Estimate optimal EGA with Leiden and CPM
#' fit.leiden <- EGA.fit(
#'   data = wmt, algorithm = "leiden",
#'   objective_function = "CPM", # default
#'   # resolution_parameter = seq.int(0, max(abs(network)), 0.01),
#'   # For CPM, the default max resolution parameter
#'   # is set to the largest absolute edge in the network
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#' 
#' # Estimate optimal EGA with Leiden and modularity
#' fit.leiden <- EGA.fit(
#'   data = wmt, algorithm = "leiden",
#'   objective_function = "modularity",
#'   resolution_parameter = seq.int(0, 2, 0.05),
#'   # default for modularity
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#' 
#' \dontrun{
#' # Estimate optimal EGA with Louvain
#' fit.louvain <- EGA.fit(
#'   data = wmt, algorithm = "louvain",
#'   resolution_parameter = seq.int(0, 2, 0.05), # default
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )}
#'
#' @references
#' \strong{Entropy fit measures} \cr
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Neito, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (in press).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#' 
#' \strong{Simulation for EGA.fit} \cr
#' Jamison, L., Christensen, A. P., & Golino, H. (under review).
#' Optimizing Walktrap's community detection in networks using the Total Entropy Fit Index.
#' \emph{PsyArXiv}.
#' 
#' \strong{Leiden algorithm} \cr
#' Traag, V. A., Waltman, L., & Van Eck, N. J. (2019).
#' From Louvain to Leiden: guaranteeing well-connected communities.
#' \emph{Scientific Reports}, \emph{9}(1), 1-12.
#' 
#' \strong{Louvain algorithm} \cr
#' Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008). 
#' Fast unfolding of communities in large networks. 
#' \emph{Journal of Statistical Mechanics: Theory and Experiment}, \emph{2008}(10), P10008.
#' 
#' \strong{Walktrap algorithm} \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#'
#' @export
#' 
# EGA fit ----
# Updated 06.08.2023
EGA.fit <- function(
    data, n = NULL,
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    plot.EGA = TRUE, verbose = FALSE,
    ...
)
{
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", c("auto", "cor_auto", "pearson", "spearman"))
  corr <- swiftelse(corr == "cor_auto", "auto", corr) # deprecate `cor_auto`
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", EGA.fit)
  
  # Argument errors
  EGA.fit_errors(data, n, plot.EGA, verbose)

  # Obtain ellipse arguments
  ellipse <- list(...)
  
  # Get fit function
  fit_FUN <- switch(
    algorithm, # other algorithms error in `set_default`
    "walktrap" = walktrap_fit,
    "louvain" = louvain_fit,
    "leiden" = leiden_fit
  )
  
  # Get result
  ## `EGA.estimate` will handle the data processing (e.g., correlations)
  fit_result <- fit_FUN(
    data = data, n = n, corr = corr,
    na.data = na.data, model = model, 
    verbose = verbose, ellipse = ellipse
  )
  
  # Set objective function (will be NULL for Louvain and Walktrap)
  objective_function <- fit_result$objective_function
  
  # Obtain only unique solutions
  search_unique <- unique_solutions(fit_result$search_matrix)
  
  # Obtain absolute correlation matrix
  if(model != "bggm"){
    correlation_matrix <- fit_result$ega_result$cor.data
  }else{
    
    # Needs correlation matrix for BGGM
    correlation_matrix <- obtain_sample_correlations(
      data = data, n = 1, # avoids throwing error
      corr = corr, na.data = na.data,
      verbose = verbose, ...
    )$correlation_matrix
    
  }
  
  # Determine best fitting solution
  fit_values <- apply(
    search_unique, 1, function(membership){
      tefi(correlation_matrix, membership)
    }$VN.Entropy.Fit
  )
  
  # Determine best fit index
  best_index <- which.min(fit_values)
  
  # Obtain best solution
  best_solution <- search_unique[best_index,]
  
  # Add methods to membership attributes
  attr(best_solution, "methods") <- list(
    algorithm = obtain_algorithm_name(algorithm),
    # `obtain_algorithm_name` is with `community.detection`
    objective_function = objective_function
  )
  
  # Set class (proper `EGA.estimate` printing)
  class(best_solution) <- "EGA.community"
  
  # Set up EGA with appropriate memberships
  fit_result$ega_result$wc <- best_solution
  fit_result$ega_result$n.dim <- unique_length(best_solution)
  
  # Set up best fit results
  best_fit <- list(
    EGA = fit_result$ega_result,
    EntropyFit = fit_values,
    Lowest.EntropyFit = fit_values[best_index]
  )
  
  # Ensure consistent naming with `EGA`
  if("cor.data" %in% names(best_fit$EGA)){
    names(best_fit$EGA)[
      names(best_fit$EGA) == "cor.data"
    ] <- "correlation"
  }
  
  # Add parameters
  best_fit$parameter.space <- fit_result$parameters
  
  # No need for attributes (all necessary information S3 is available)
  
  # Set class
  class(best_fit) <- "EGA.fit"
  
  # Check for plot
  if(plot.EGA && sum(best_fit$EGA$network != 0)){
    
    # Set up plot
    best_fit$Plot.EGA <- plot(best_fit, ...)
    
    # Actually send the plot
    silent_plot(best_fit$Plot.EGA)
    
    
  }
  
  # Return results
  return(best_fit)
  
}

#' @noRd
# Errors ----
# Updated 30.07.2023
EGA.fit_errors <- function(data, n, plot.EGA, verbose)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame"))
  
  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1)
    typeof_error(n, "numeric")
  }
  
  # 'plot.EGA' errors
  length_error(plot.EGA, 1)
  typeof_error(plot.EGA, "logical")
  
  # 'verbose' errors
  length_error(verbose, 1)
  typeof_error(verbose, "logical")
  
}

# Bug checking ----
## Basic input
# data = wmt2[,7:24]; n = NULL; corr = "auto"
# na.data = "pairwise"; model = "glasso"; algorithm = "louvain"
# plot.EGA = TRUE; verbose = FALSE
# ellipse = list()

#' @exportS3Method 
# S3 Print Method ----
# Updated 27.06.2023
print.EGA.fit <- function(x, ...)
{
  
  # Print network estimation
  print(x$EGA$network)
  
  # Add break space
  cat("\n----\n\n")
  
  # Obtain membership
  membership <- x$EGA$wc
  
  # Determine number of communities
  communities <- unique_length(membership) 
  
  # Obtain attributes
  community_attributes <- attr(membership, "methods")
  
  # Obtain algorithm name (if available)
  algorithm <- community_attributes$algorithm
  
  # Check for Leiden
  if(algorithm == "Leiden"){
    
    # Obtain objective function
    objective_function <- community_attributes$objective_function
    
    # Set up algorithm name
    objective_name <- swiftelse(
      is.null(objective_function),
      "CPM", objective_function
    )
    
    # Expand "CPM"
    objective_name <- swiftelse(
      objective_name == "CPM",
      "Constant Potts Model", "Modularity"
    )
    
    # Finalize algorithm name
    algorithm <- paste(algorithm, "with", objective_name)
    
  }
  
  # Check for parameter addition
  algorithm <- paste0(
    algorithm,
    paste0(
      " (", swiftelse(
        algorithm == "Walktrap",
        "Steps = ",
        "Resolution = "
      ), names(x$Lowest.EntropyFit),
      ")"
    )
  )
  
  # Set up methods
  cat(paste0("Algorithm: "), algorithm)
  
  # Add breakspace
  cat("\n\n")
  
  # Print communities
  cat(paste0("Number of communities: "), communities)
  cat("\n\n") # Add breakspace
  
  # Remove class and attribute for clean print
  membership <- remove_attributes(membership)
  
  # Print membership
  print(membership)
  
  # Add breakspace
  cat("\n\n----\n\n")
  
  # Add TEFI value
  cat(paste0("TEFI: ", round(x$Lowest.EntropyFit, 3)))
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 05.07.2023
summary.EGA.fit <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Plot Method ----
# Updated 23.06.2023
plot.EGA.fit <- function(x, ...)
{
  
  # Return plot
  single_plot(
    network = x$EGA$network,
    wc = x$EGA$wc,
    ...
  )
  
}

#' @noRd
# Determine unique solutions ----
# Updated 07.07.2023 
unique_solutions <- function(search_matrix)
{

  # Obtain counts to eliminate duplicates
  search_unique <- unique(search_matrix, MARGIN = 1)
  
  # Obtain original rows
  original_rows <- dim(search_unique)[1]
  
  # Remove solutions with singleton
  search_unique <- na.omit(search_unique)

  # Singleton rows
  singleton_rows <- dim(search_unique)[1]
  
  # Remove unidimensional solutions
  search_unique <- search_unique[
    apply(search_unique, 1, function(x){unique_length(x) != 1}),, drop = FALSE
  ]
  
  # Unidimensional rows
  unidimensional_rows <- dim(search_unique)[1]
  
  # Check for errors to send
  if(singleton_rows == 0){
    stop("All solutions were determined to include at least one singleton community", call. = FALSE)
  }else if(singleton_rows == 1 & unidimensional_rows == 0){
    stop("All solutions were determined to be unidimensional and/or at least one singleton community", call. = FALSE)
  }else if(unidimensional_rows == 0){
    stop("All solutions were determined to be unidimensional", call. = FALSE)
  }
  
  # Return unique solutions
  return(search_unique)

}

#' @noRd
# Fit for Walktrap ----
# Updated 07.07.2023
walktrap_fit <- function(
    data, n, corr, na.data, model,
    verbose, ellipse
)
{
  
  # Check for parameter search space
  if(!"steps" %in% names(ellipse)){
    steps <- 3:8 # default
  }else{
    
    # Set steps
    steps <- ellipse$steps
    
    # Remove objective function from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "steps"]
    
  }
  
  # Get the length of the steps
  step_length <- length(steps)
  
  # Perform EGA with first parameter
  ega_result <- do.call(
    what = EGA.estimate,
    args = c(
      list( # Necessary call
        data = data, n = n, corr = corr,
        na.data = na.data, model = model,
        algorithm = "walktrap", verbose = verbose,
        steps = steps[1]
      ),
      ellipse # pass on ellipse
    )
  )
  
  # Set up search matrix
  search_matrix <- matrix(
    nrow = step_length,
    ncol = length(ega_result$wc),
    dimnames = list(steps, names(ega_result$wc))
  )
  
  # Add first parameter
  search_matrix[1,] <- ega_result$wc
  
  # The network won't change so apply the community detection
  # algorithm over the rest of the parameters
  for(i in 2:step_length){
    search_matrix[i,] <- do.call(
      what = community.detection,
      args = c(
        list( # Necessary call
          network = ega_result$network,
          algorithm = "walktrap",
          steps = steps[i]
        ),
        ellipse # pass on ellipse
      )
    )
  }
  
  # Return results
  return(
    list(
      ega_result = ega_result,
      search_matrix = search_matrix,
      parameters = steps
    )
  )
  
}

#' @noRd
# Fit for Louvain ----
# Updated 07.07.2023
louvain_fit <- function(
    data, n, corr, na.data, model,
    verbose, ellipse
)
{
  
  # Check for parameter search space
  if(!"resolution_parameter" %in% names(ellipse)){
    resolution_parameter <- seq.int(0, 2, 0.05) # default
  }else{
    
    # Set resolution parameter
    resolution_parameter <- ellipse$resolution_parameter
    
    # Remove resolution parameter from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "resolution_parameter"]
    
  }
  
  # Get the length of the resolution parameter
  resolution_parameter_length <- length(resolution_parameter)
  
  # Perform EGA with first parameter
  ega_result <- do.call(
    what = EGA.estimate,
    args = c(
      list( # Necessary call
        data = data, n = n, corr = corr,
        na.data = na.data, model = model,
        algorithm = "louvain", verbose = verbose,
        resolution_parameter = resolution_parameter[1]
      ),
      ellipse # pass on ellipse
    )
  )
  
  # Set up search matrix
  search_matrix <- matrix(
    nrow = resolution_parameter_length,
    ncol = length(ega_result$wc),
    dimnames = list(resolution_parameter, names(ega_result$wc))
  )
  
  # Add first parameter
  search_matrix[1,] <- ega_result$wc
  
  # The network won't change so apply the community detection
  # algorithm over the rest of the parameters
  for(i in 2:resolution_parameter_length){
    search_matrix[i,] <- do.call(
      what = community.consensus,
      args = c(
        list( # Necessary call
          network = ega_result$network,
          resolution = resolution_parameter[i]
        ),
        ellipse # pass on ellipse
      )
    )
  }
  
  # Return results
  return(
    list(
      ega_result = ega_result,
      search_matrix = search_matrix,
      parameters = resolution_parameter
    )
  )
  
}

#' @noRd
# Fit for Leiden ----
# Updated 20.07.2023
leiden_fit <- function(
    data, n, corr, na.data, model,
    verbose, ellipse
)
{
  
  # Check for objective function
  if(!"objective_function" %in% names(ellipse)){
    
    # Default for {EGAnet}
    objective_function <- "modularity"
    
    # Send warning
    warning(
      paste0(
        "{EGAnet} uses \"modularity\" as the default objective function in the Leiden algorithm. ",
        "In contrast, {igraph} uses \"CPM\". Set `objective_function = \"CPM\"` to use the Constant Potts ",
        "Model in {EGAnet}"
      )
    )
    
  }else{
    
    # Set objective function
    objective_function <- ellipse$objective_function
    
    # Remove objective function from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "objective_function"]
    
  }
  
  # Perform EGA with first parameter
  ega_result <- do.call(
    what = EGA.estimate,
    args = c(
      list( # Necessary call
        data = data, n = n, corr = corr,
        na.data = na.data, model = model,
        algorithm = "leiden", verbose = verbose,
        resolution_parameter = 0,
        # start with resolution parameter at zero
        # zero is guaranteed to be unidimensional
        objective_function = objective_function
      ),
      ellipse # pass on ellipse
    )
  )
  
  # Check for parameter search space
  if(!"resolution_parameter" %in% names(ellipse)){
    
    # Switch based on objective function
    if(objective_function == "CPM"){
      resolution_parameter <- seq.int(0, max(abs(ega_result$network)), 0.01) # default 
    }else if(objective_function == "modularity"){
      resolution_parameter <- seq.int(0, 2, 0.05) # default
    }
    
  }else{
    
    # Set resolution parameter
    resolution_parameter <- ellipse$resolution_parameter
    
    # Remove resolution parameter from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "resolution_parameter"]
    
  }
  
  # Get the length of the resolution parameter
  resolution_parameter_length <- length(resolution_parameter)
  
  # Set up search matrix
  search_matrix <- matrix(
    nrow = resolution_parameter_length,
    ncol = length(ega_result$wc),
    dimnames = list(resolution_parameter, names(ega_result$wc))
  )
  
  # Add first parameter
  search_matrix[1,] <- ega_result$wc
  
  # The network won't change so apply the community detection
  # algorithm over the rest of the parameters
  for(i in 2:resolution_parameter_length){
    search_matrix[i,] <- do.call(
      what = community.detection,
      args = c(
        list( # Necessary call
          network = ega_result$network,
          algorithm = "leiden",
          resolution_parameter = resolution_parameter[i],
          objective_function = objective_function
        ),
        ellipse # pass on ellipse
      )
    )
  }
  
  # Return results
  return(
    list(
      ega_result = ega_result,
      search_matrix = search_matrix,
      parameters = resolution_parameter,
      objective_function = objective_function
    )
  )
  
}

