#' \code{\link[EGAnet]{EGA}} Optimal Model Fit using the Total Entropy Fit Index  (\code{\link[EGAnet]{tefi}})
#'
#' @description Estimates the best fitting model using \code{\link[EGAnet]{EGA}}.
#' The number of steps in the \code{\link[igraph]{cluster_walktrap}} detection
#' algorithm is varied and unique community solutions are compared using
#' \code{\link[EGAnet]{tefi}}.
#'
#' @param data Numeric matrix or data frame.
#' Either data representing \emph{only} the variables of interest or
#' a correlation matrix. Data that are not numeric will be
#' removed from the dataset
#'
#' @param n Numeric (length = 1).
#' Sample size if \code{data} is a correlation matrix
#'
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"} to automatically compute
#' appropriate correlations using \code{\link[EGAnet]{auto.correlate}}.
#' \code{"pearson"} and \code{"spearman"} are provide for completeness.
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#' 
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"}}
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"}}
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param model Character (length = 1).
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"BGGM"}}
#' {Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{\link[BGGM]{estimate}} for more details}
#' 
#' \item{\code{"glasso"}}
#' {Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details}
#' 
#' \item{\code{"TMFG"}}
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
#' \item{\code{"leiden"}}
#' {See \code{\link[igraph]{cluster_leiden}} for more details.
#' \emph{Note}: The Leiden algorithm will default to the
#' Constant Potts Model objective function
#' (\code{objective_function = "CPM"}). Set
#' \code{objective_function = "modularity"} to use
#' modularity instead (see examples). By default, searches along
#' resolutions from 0 to 2 in 0.001 increments
#' (\code{resolution_parameter = seq.int(0, 2, 0.001)}). Use the argument \code{resolution_parameter}
#' to change the search parameters}
#' 
#' \item{\code{"louvain"}}
#' {See \code{\link[igraph]{cluster_louvain}} for more details.
#' By default, searches along resolutions from 0 to 2 in 0.001 increments
#' (\code{resolution_parameter = seq.int(0, 2, 0.001)}). Use the argument \code{resolution_parameter}
#' to change the search parameters}
#' 
#' \item{\code{"walktrap"}}
#' {This algorithm is the default. See \code{\link[EGAnet]{cluster_walktrap}} for more details.
#' By default, searches along 3 to 8 steps (\code{steps = 3:8}). Use the argument \code{steps}
#' to change the search parameters}
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
#' \code{\link[EGAnet]{community.detection}}, and \code{\link[EGAnet]{EGA.estimate}}
#'
#' @return Returns a list containing:
#'
#' \item{EGA}{The \code{\link[EGAnet]{EGA}} output for the best fitting model}
#'
#' \item{steps}{The number of steps used in the best fitting model from
#' the \code{\link[igraph]{cluster_walktrap}} algorithm}
#' 
#' \item{resolution_parameter}{The resolution parameter used in the best fitting model from
#' the \code{\link[igraph]{cluster_leiden}} algorithm}
#'
#' \item{EntropyFit}{The \code{\link[EGAnet]{tefi}} Index for the unique solutions given the range of steps
#' (vector names represent the number of steps)}
#'
#' \item{Lowest.EntropyFit}{The lowest value for the \code{\link[EGAnet]{tefi}} Index}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Estimate optimal EGA with Walktrap
#' fit.walktrap <- EGA.fit(
#'   data = wmt, algorithm = "walktrap",
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#' 
#' # Estimate optimal EGA with Louvain
#' fit.louvain <- EGA.fit(
#'   data = wmt, algorithm = "louvain",
#'   resolution_parameter = seq.int(0, 2, 0.50),
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#' 
#' # Estimate optimal EGA with Leiden and modularity
#' fit.leiden <- EGA.fit(
#'   data = wmt, algorithm = "leiden",
#'   objective_function = "modularity",
#'   resolution_parameter = seq.int(0, 2, 0.50),
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#'
#' @references
#' # Entropy fit measures \cr
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Neito, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (in press).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#' 
#' # Simulation for EGA.fit \cr
#' Jamison, L., Christensen, A. P., & Golino, H. (under review).
#' Optimizing Walktrap's community detection in networks using the Total Entropy Fit Index.
#' \emph{PsyArXiv}.
#' 
#' # Leiden algorithm \cr
#' Traag, V. A., Waltman, L., & Van Eck, N. J. (2019).
#' From Louvain to Leiden: guaranteeing well-connected communities.
#' \emph{Scientific Reports}, \emph{9}(1), 1-12.
#' 
#' # Walktrap algorithm \cr
#' Pons, P., & Latapy, M. (2006).
#' Computing communities in large networks using random walks.
#' \emph{Journal of Graph Algorithms and Applications}, \emph{10}, 191-218.
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# EGA fit
# Updated 26.06.2023
EGA.fit <- function (
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
  corr <- ifelse(corr == "cor_auto", "auto", corr) # deprecate `cor_auto`
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", EGA.fit)

  # Obtain ellipse arguments
  ellipse <- list(...)
  
  # `EGA.estimate` will handle the data processing (e.g., correlations)
  ## Branch for Walktrap or Leiden/Louvain
  if(algorithm == "walktrap"){
    
    # Set objective function to NULL
    objective_function <- NULL
    
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
          algorithm = algorithm, verbose = verbose,
          steps = steps[1]
        ),
        ellipse # pass on ellipse
      )
    )
    
    # Set up search matrix
    search_matrix <- matrix(
      0, nrow = step_length,
      ncol = length(ega_result$wc)
    )
    
    # Add names
    colnames(search_matrix) <- names(ega_result$wc)
    row.names(search_matrix) <- steps
    
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
            algorithm = algorithm,
            steps = steps[i]
          ),
          ellipse # pass on ellipse
        )
      )
    }
    
  }else{
    
    # Check for Leiden algorithm
    if(algorithm == "leiden"){

      # If Leiden, then check for objective function
      if(!"objective_function" %in% names(ellipse)){
        objective_function <- "CPM"
        # default for {igraph} and `community.detection`
      }else{
        
        # Set objective function
        objective_function <- ellipse$objective_function
        
        # Remove objective function from ellipse to
        # make other arguments available in `do.call`
        ellipse <- ellipse[names(ellipse) != "objective_function"]
        
      }

    }else{ # Set to NULL to avoid conflict with Louvain
      objective_function <- NULL
    }
    
    # Check for parameter search space
    if(!"resolution_parameter" %in% names(ellipse)){
      resolution_parameter <- seq.int(0, 2, 0.001) # default
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
          algorithm = algorithm, verbose = verbose,
          resolution_parameter = resolution_parameter[1],
          objective_function = objective_function
        ),
        ellipse # pass on ellipse
      )
    )
    
    # Set up search matrix
    search_matrix <- matrix(
      0, nrow = resolution_parameter_length,
      ncol = length(ega_result$wc),
      dimnames = list(
        resolution_parameter, names(ega_result$wc)
      )
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
            algorithm = algorithm,
            resolution_parameter = resolution_parameter[i],
            objective_function = objective_function
          ),
          ellipse # pass on ellipse
        )
      )
    }
    
  }
  
  # Obtain only unique solutions
  search_unique <- unique_solutions(search_matrix)
  
  # Obtain absolute correlation matrix
  if(model != "bggm"){
    correlation_matrix <- ega_result$cor.data
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
  
  # Obtain signed argument
  if(!"signed" %in% names(ellipse)){
    signed <- FALSE # default
  }else{signed <- ellipse$signed}
  
  # Add methods to membership attributes
  attr(best_solution, "methods") <- list(
    algorithm = obtain_algorithm_name(algorithm),
    # `obtain_algorithm_name` is with `community.detection`
    signed = signed, objective_function = objective_function
  )
  
  # Set class (proper `EGA.estimate` printing)
  class(best_solution) <- "EGA.community"
  
  # Set up EGA with appropriate memberships
  ega_result$wc <- best_solution
  ega_result$n.dim <- unique_length(best_solution)
  
  # Set up best fit results
  best_fit <- list(
    EGA = ega_result,
    EntropyFit = fit_values,
    Lowest.EntropyFit = fit_values[best_index]
  )
  
  # Add parameters
  if(algorithm == "walktrap"){
    best_fit$parameter.space <- steps
  }else{
    best_fit$parameter.space <- resolution_parameter
  }
  
  # No need for attributes (all necessary information S3 is available)
  
  # Set class
  class(best_fit) <- "EGA.fit"
  
  # Check for plot
  if(isTRUE(plot.EGA)){
    
    # Set up plot
    best_fit$Plot.EGA <- plot(best_fit)
    
    # Actually send the plot
    plot(best_fit$Plot.EGA)
    
  }
      
  return(best_fit)
}

# Bug checking ----
## Basic input
# data = wmt2[,7:24]; n = NULL; corr = "auto"
# na.data = "pairwise"; model = "glasso"; algorithm = "leiden"
# plot.EGA = TRUE; verbose = FALSE
# ellipse = list(objective_function = "modularity")

#' @exportS3Method 
# S3 Print Method ----
# Updated 23.06.2023
print.EGA.fit <- function(x, ...)
{
  
  # Print network estimation
  print(x$EGA$network)
  
  # Add break space
  cat("\n----\n\n")
  
  # Obtain membership
  membership <- x$EGA$wc
  
  # Determine number of communities
  communities <- length(na.omit(unique(membership)))
  
  # Obtain algorithm name (if available)
  algorithm <- attr(membership, "methods")$algorithm
  
  # Check for signed
  algorithm_name <- ifelse(
    attr(membership, "methods")$signed,
    paste("Signed", algorithm),
    algorithm
  )
  
  # Check for Leiden
  if(algorithm == "Leiden"){
    
    # Obtain objective function
    objective_function <- attr(membership, "methods")$objective_function
    
    # Set up algorithm name
    objective_name <- ifelse(
      is.null(objective_function),
      "CPM", objective_function
    )
    
    # Expand "CPM"
    objective_name <- ifelse(
      objective_name == "CPM",
      "Constant Potts Model", "Modularity"
    )
    
    # Finalize algorithm name
    algorithm_name <- paste(
      algorithm, "with", objective_name
    )
    
  }
  
  # Check for parameter addition
  algorithm_name <- paste0(
    algorithm_name,
    paste0(
      " (", ifelse(
        algorithm == "Walktrap",
        "Steps = ",
        "Resolution = "
      ), names(x$Lowest.EntropyFit),
      ")"
    )
  )
  
  # Set up methods
  cat(paste0("Algorithm: "), algorithm_name)
  
  # Add breakspace
  cat("\n\n")
  
  # Add TEFI value
  cat(paste0("TEFI: ", round(x$Lowest.EntropyFit, 3)))
  
  # Add breakspace
  cat("\n\n")
  
  # Print communities
  cat(paste0("Number of communities: "), communities)
  cat("\n\n") # Add breakspace
  
  # Remove attribute for clean print
  attr(membership, which = "class") <- NULL
  attr(membership, which = "methods") <- NULL
  
  # Print membership
  print(membership)
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 23.06.2023
summary.EGA.fit <- function(object, ...)
{
  
  # Print network estimation
  print(object$EGA$network)
  
  # Add break space
  cat("\n----\n\n")
  
  # Obtain membership
  membership <- object$EGA$wc
  
  # Determine number of communities
  communities <- length(na.omit(unique(membership)))
  
  # Obtain algorithm name (if available)
  algorithm <- attr(membership, "methods")$algorithm
  
  # Check for signed
  algorithm_name <- ifelse(
    attr(membership, "methods")$signed,
    paste("Signed", algorithm),
    algorithm
  )
  
  # Check for Leiden
  if(algorithm == "Leiden"){
    
    # Obtain objective function
    objective_function <- attr(membership, "methods")$objective_function
    
    # Set up algorithm name
    objective_name <- ifelse(
      is.null(objective_function),
      "CPM", objective_function
    )
    
    # Expand "CPM"
    objective_name <- ifelse(
      objective_name == "CPM",
      "Constant Potts Model", "Modularity"
    )
    
    # Finalize algorithm name
    algorithm_name <- paste(
      algorithm, "with", objective_name
    )
    
  }
  
  # Check for parameter addition
  algorithm_name <- paste0(
    algorithm_name,
    paste0(
      " (", ifelse(
        algorithm == "Walktrap",
        "Steps = ",
        "Resolution = "
      ), names(object$Lowest.EntropyFit),
      ")"
    )
  )
  
  # Set up methods
  cat(paste0("Algorithm: "), algorithm_name)
  
  # Add breakspace
  cat("\n\n")
  
  # Add TEFI value
  cat(paste0("TEFI: ", round(object$Lowest.EntropyFit, 3)))
  
  # Add breakspace
  cat("\n\n")
  
  # Print communities
  cat(paste0("Number of communities: "), communities)
  cat("\n\n") # Add breakspace
  
  # Remove attribute for clean print
  attr(membership, which = "class") <- NULL
  attr(membership, which = "methods") <- NULL
  
  # Print membership
  print(membership)
  
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
# Updated 26.06.2023 
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
    stop("All solutions were determined to include at least one singleton community")
  }else if(singleton_rows == 1 & unidimensional_rows == 0){
    stop("All solutions were determined to be unidimensional and/or at least one singleton community")
  }else if(unidimensional_rows == 0){
    stop("All solutions were determined to be unidimensional")
  }
  
  # Return unique solutions
  return(search_unique)

}




