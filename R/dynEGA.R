#' @title Dynamic Exploratory Graph Analysis
#'
#' @description Estimates dynamic communities in multivariate time series 
#' (e.g., panel data, longitudinal data, intensive longitudinal data) at multiple
#' time scales and at different levels of analysis: 
#' individuals (intraindividual structure), groups, and population (interindividual structure)
#' 
#' @param data Matrix or data frame.
#' Participants and variable should be in long format such that
#' row \emph{t} represents observations for all variables at time point 
#' \emph{t} for a participant. The next row, \emph{t + 1}, represents
#' the next measurement occasion for that same participant. The next
#' participant's data should immediately follow, in the same pattern,
#' after the previous participant
#' 
#' \code{data} should have an ID variable labeled \code{"ID"}; otherwise, it is
#' assumed that the data represent the population
#' 
#' For groups, \code{data} should have a Group variable labeled \code{"Group"};
#' otherwise, it is assumed that there are no groups in \code{data}
#' 
#' Arguments \code{id} and \code{group} can be specified to tell the function
#' which column in \code{data} it should use as the ID and Group variable, respectively
#' 
#' A measurement occasion variable is not necessary and should be \emph{removed}
#' from the data before proceeding with the analysis
#' 
#' @param id Numeric or character (length = 1).
#' Number or name of the column identifying each individual.
#' Defaults to \code{NULL}
#'
#' @param group Numeric or character (length = 1).
#' Number of the column identifying group membership.
#' Defaults to \code{NULL}
#'
#' @param n.embed Numeric (length = 1).
#' Defaults to \code{5}.
#' Number of embedded dimensions (the number of observations to 
#' be used in the \code{\link[EGAnet]{Embed}} function). For example,
#' an \code{"n.embed = 5"} will use five consecutive observations 
#' to estimate a single derivative
#'
#' @param tau Numeric (length = 1).
#' Defaults to \code{1}.
#' Number of observations to offset successive embeddings in 
#' the \code{\link[EGAnet]{Embed}} function.
#' Generally recommended to leave "as is"
#'
#' @param delta Numeric (length = 1).
#' Defaults to \code{1}.
#' The time between successive observations in the time series (i.e, lag).
#' Generally recommended to leave "as is"
#' 
#' @param use.derivatives Numeric (length = 1).
#' Defaults to \code{1}.
#' The order of the derivative to be used in the analysis.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{0} --- }
#' {No derivatives; consistent with moving average}
#' 
#' \item{\code{1} --- }
#' {First-order derivatives; interpretated as "velocity" or 
#' rate of change over time}
#' 
#' \item{\code{2} --- }
#' {Second-order derivatives; interpreted as "acceleration" or
#' rate of the rate of change over time}
#' 
#' }
#' 
#' Generally recommended to leave "as is"
#'
#' @param level Character vector (up to length of 3).
#' A character vector indicating which level(s) to estimate:
#'
#' \itemize{
#'
#' \item{\code{"individual"} --- }
#' {Estimates \code{\link[EGAnet]{EGA}} for each individual in \code{data}
#' (intraindividual structure; requires an \code{"ID"} column, see \code{data})}
#' 
#' \item{\code{"group"} --- }
#' {Estimates \code{\link[EGAnet]{EGA}} for each group in \code{data}
#' (group structure; requires a \code{"Group"} column, see \code{data})}
#' 
#' \item{\code{"population"} --- }
#' {Estimates \code{\link[EGAnet]{EGA}} across all \code{data}
#' (interindividual structure)}
#'
#' }
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
#' \item{\code{"cor_auto"} --- }
#' {Uses \code{\link[qgraph]{cor_auto}} to compute correlations. Arguments
#' can be passed along to the function}
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
#' @param algorithm Character or 
#' \code{\link{igraph}} \code{cluster_*} function (length = 1).
#' Defaults to \code{"walktrap"}.
#' Three options are listed below but all are available
#' (see \code{\link[EGAnet]{community.detection}} for other options):
#' 
#' \itemize{
#'
#' \item{\code{"leiden"} --- }
#' {See \code{\link[igraph]{cluster_leiden}} for more details}
#' 
#' \item{\code{"louvain"} --- }
#' {By default, \code{"louvain"} will implement the Louvain algorithm using 
#' the consensus clustering method (see \code{\link[EGAnet]{community.consensus}} 
#' for more information). This function will implement
#' \code{consensus.method = "most_common"} and \code{consensus.iter = 1000} 
#' unless specified otherwise}
#' 
#' \item{\code{"walktrap"} --- }
#' {See \code{\link[igraph]{cluster_walktrap}} for more details}
#' 
#' }
#'
#' @param uni.method Character (length = 1).
#' What unidimensionality method should be used? 
#' Defaults to \code{"louvain"}.
#' Available options:
#' 
#' \itemize{
#'
#' \item{\code{"expand"} --- }
#' {Expands the correlation matrix with four variables correlated 0.50.
#' If number of dimension returns 2 or less in check, then the data 
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This method was used in the Golino et al.'s (2020)
#' \emph{Psychological Methods} simulation}
#'
#' \item{\code{"LE"} --- }
#' {Applies the Leading Eigenvector algorithm
#' (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvector solution is used; otherwise, regular EGA
#' is used. This method was used in the Christensen et al.'s (2023) 
#' \emph{Behavior Research Methods} simulation}
#' 
#' \item{\code{"louvain"} --- }
#' {Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix. If the number of dimensions is 1, 
#' then the Louvain solution is used; otherwise, regular EGA is used. 
#' This method was validated Christensen's (2022) \emph{PsyArXiv} simulation.
#' Consensus clustering can be used by specifying either
#' \code{"consensus.method"} or \code{"consensus.iter"}}
#' 
#' }
#' 
#' @param ncores Numeric (length = 1).
#' Number of cores to use in computing results.
#' Defaults to \code{ceiling(parallel::detectCores() / 2)} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' If you're unsure how many cores your computer has,
#' then type: \code{parallel::detectCores()}
#'
#' @param verbose Boolean (length = 1).
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not display progress
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}}, and
#' \code{\link[EGAnet]{EGA}}
#' 
#' @return A list containing:
#' 
#' \item{Derivatives}{A list containing:
#' 
#' \itemize{
#' 
#' \item{\code{Estimates} --- }
#' {A list the length of the unique IDs containing
#' data frames of zero- to second-order derivatives for each ID in \code{data}}
#' 
#' \item{\code{EstimatesDF} --- }
#' {A data frame of derivatives across all IDs containing
#' columns of the zero- to second-order derivatives as well as \code{id} and
#' \code{group} variables (\code{group} is automatically set to \code{1}
#' for all if no \code{group} is provided)}
#' 
#' }
#' 
#' }
#' 
#' \item{dynEGA}{A list containing:
#' 
#' \itemize{
#' 
#' \item{\code{population} --- }
#' {If \code{level} includes \code{"populaton"}, then
#' the \code{\link[EGAnet]{EGA}} results for the entire sample}
#' 
#' \item{\code{group} --- }
#' {If \code{level} includes \code{"group"}, then
#' a list containing the \code{\link[EGAnet]{EGA}} results for each \code{group}}
#' 
#' \item{\code{individual} --- }
#' {If \code{level} includes \code{"individual"}, then
#' a list containing the \code{\link[EGAnet]{EGA}} results for each \code{id}}
#' 
#' }
#' 
#' }
#' 
#' @details Derivatives for each variable's time series for each participant are
#' estimated using generalized local linear approximation (see \code{\link[EGAnet]{glla}}).
#' \code{\link[EGAnet]{EGA}} is then applied to these derivatives to model how variables
#' are changing together over time. Variables that change together over time are detected
#' as communities
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Population structure
#' simulated_population <- dynEGA(
#'   data = sim.dynEGA, level = "population"
#'   # uses simulated data in package
#'   # useful to understand how data should be structured
#' )
#'
#' # Group structure
#' simulated_group <- dynEGA(
#'   data = sim.dynEGA, level = "group"
#'   # uses simulated data in package
#'   # useful to understand how data should be structured
#' )
#' 
#' \dontrun{
#' # Individual structure
#' simulated_individual <- dynEGA(
#'   data = sim.dynEGA, level = "individual",
#'   ncores = 2, # use more for quicker results
#'   verbose = TRUE # progress bar
#' )
#' 
#' # Population, group, and individual structure
#' simulated_all <- dynEGA(
#'   data = sim.dynEGA, 
#'   level = c("individual", "group", "population"),
#'   ncores = 2, # use more for quicker results
#'   verbose = TRUE # progress bar
#' )
#' 
#' # Plot population
#' plot(simulated_all$dynEGA$population)
#' 
#' # Plot groups
#' plot(simulated_all$dynEGA$group)
#' 
#' # Plot individual
#' plot(simulated_all$dynEGA$individual, id = 1)
#' 
#' # Step through all plots
#' # Unless `id` is specified, 4 random IDs
#' # will be drawn from individuals
#' plot(simulated_all)}
#'
#' @references
#' \strong{Generalized local linear approximation} \cr
#' Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010)
#' Generalized local linear approximation of derivatives from time series. In S.-M. Chow, E. Ferrer, & F. Hsieh (Eds.),
#' \emph{The Notre Dame series on quantitative methodology. Statistical methods for modeling human dynamics: An interdisciplinary dialogue},
#' (p. 161-178). \emph{Routledge/Taylor & Francis Group}.
#'
#' Deboeck, P. R., Montpetit, M. A., Bergeman, C. S., & Boker, S. M. (2009)
#' Using derivative estimates to describe intraindividual variability at multiple time scales.
#' \emph{Psychological Methods}, \emph{14(4)}, 367-386.
#' 
#' \strong{Original dynamic EGA implementation} \cr
#' Golino, H., Christensen, A. P., Moulder, R. G., Kim, S., & Boker, S. M. (2021).
#' Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
#' \emph{Psychometrika}.
#'
#' \strong{Time delay embedding procedure} \cr
#' Savitzky, A., & Golay, M. J. (1964).
#' Smoothing and differentiation of data by simplified least squares procedures.
#' \emph{Analytical Chemistry}, \emph{36(8)}, 1627-1639.
#' 
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#'
#' @export
#' 
# dynEGA ----
# Updated 04.09.2023
dynEGA <- function(
    # `dynEGA` arguments
    data,  id = NULL, group = NULL,
    n.embed = 5, tau = 1, delta = 1, use.derivatives = 1,
    level = c("individual", "group", "population"),
    # `EGA` arguments
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    ncores, verbose = TRUE, ...
){
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", dynEGA)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}
  
  # Handle level (allow multiple to be estimated at once)
  if(missing(level)){
    level <- "population" # default to full sample
  }else{level <- match.arg(level, several.ok = TRUE)}
  
  # Argument errors (return data in case of tibble)
  data <- dynEGA_errors(
    data, id, group, n.embed, tau, delta,
    use.derivatives, ncores, verbose
  )
  
  # Get dimensions of the data
  dimensions <- dim(data)
  
  # Get "ID" and "Group" attributes of data if they exist
  # "id" and "group" columns will be added as "ID" and "Group"
  # attributes to the data
  # These columns will be removed from the data!
  data <- get_attributes(data, dimensions, id, group, level)
  
  # Get variable names
  variable_names <- dimnames(data)[[2]]
  
  # Split data into lists based on ID
  individual_data <- split(data, attributes(data)$ID)
  
  # Set up to compute GLLA
  # Avoids computation of weights participant x variable times
  # Leads to about 6x faster computation
  L <- glla_setup(n.embed, tau, delta, order = 2)
  
  # Get derivatives for each participant
  derivative_list <- individual_derivatives(
    individual_data, variable_names, n.embed, tau, L,
    individual_attributes = attributes(data)
  )
  
  # Set up return list
  results <- list(
    Derivatives = list(
      Estimates = derivative_list,
      EstimatesDF = data.frame(
        do.call(rbind, derivative_list),
        id = ulapply(derivative_list, attr, "ID"),
        group = ulapply(derivative_list, attr, "Group")
      )
    )
  )
  
  # Get derivatives to use
  derivative_index <- grep(
    switch(
      as.character(use.derivatives),
      "0" = "Ord0", "1" = "Ord1", "2" = "Ord2"
    ),
    dimnames(results$Derivatives$EstimatesDF)[[2]]
  )
  
  # Compute Dynamic EGA at each level
  ## Population
  if("population" %in% level){
    
    # Estimate population EGA
    results$dynEGA$population <- EGA(
      results$Derivatives$EstimatesDF[,derivative_index],
      corr = corr, na.data = na.data, model = model,
      algorithm = algorithm, uni.method = uni.method,
      plot.EGA = FALSE, verbose = FALSE, ...
    )
    
    # Add class
    class(results$dynEGA$population) <- "dynEGA.Population"
    
  }
  
  ## Group
  if("group" %in% level){
    
    # Split derivatives by Group
    group_data <- lapply(
      split(
        derivative_list, ulapply(
          derivative_list, function(x){
            unique(attributes(x)$Group)
          }
        )
      ), do.call, what = rbind
    )
    
    # Estimate group EGA
    results$dynEGA$group <- lapply(
      group_data, function(x){
        EGA(
          x[,derivative_index],
          corr = corr, na.data = na.data, model = model,
          algorithm = algorithm, uni.method = uni.method,
          plot.EGA = FALSE, verbose = FALSE, ...
        )
      }
    )
    
    # Add class
    class(results$dynEGA$group) <- "dynEGA.Group"
    
  }
  
  ## Individual
  if("individual" %in% level){
    
    # Get proper derivatives and usable data
    # Add zero variance variables as attributes
    usable_derivatives <- lapply(
      derivative_list, function(x){
        
        # First, get proper derivatives
        proper_derivatives <- x[,derivative_index]
        
        # Next, check for zero variance variables
        zero_variance <- lvapply(
          as.data.frame(proper_derivatives),
          function(x){sd(x, na.rm = TRUE) == 0}
        )
        
        # Return data with all non-zero variance variables
        # and with zero variance variables as an attribute
        return(
          structure(
            proper_derivatives[,!zero_variance],
            ID = unique(attr(x, "ID")),
            zero_variance = zero_variance
          )
        )
        
      }
    )
    
    # Estimate individual EGA
    individual_results <- parallel_process(
      iterations = length(usable_derivatives),
      datalist = usable_derivatives,
      EGA, # Use `EGA`
      corr = corr, na.data = na.data, model = model,
      algorithm = algorithm, uni.method = uni.method,
      plot.EGA = FALSE, verbose = FALSE, ...,
      ncores = ncores, progress = verbose
    )
    
    # Process individual results
    # Handle zero variance variables within each individual
    results$dynEGA$individual <- process_individual_dynEGA(
      usable_derivatives, individual_results
    )
    
    # Add class
    class(results$dynEGA$individual) <- "dynEGA.Individual"
    
  }
  
  # Add attributes to overall results
  ## EGA attributes will already be attached to `results$dynEGA`
  attr(results, "glla") <- list(
    n.embed = n.embed, tau = tau, delta = delta,
    use.derivatives = use.derivatives
  )
  
  # Add overall class
  class(results) <- "dynEGA"
  
  # Return results
  return(results)
  
}

# Bug Checking ----
# data = sim.dynEGA; n.embed = 5; tau = 1; delta = 1
# level = c("individual", "group", "population")
# id = NULL; group = NULL; use.derivatives = 1
# corr = "auto"; na.data = "pairwise"
# model = "glasso"; algorithm = "walktrap"
# uni.method = "louvain"; ncores = 8
# verbose = FALSE; ellipse = list()

#' @noRd
# Errors ----
# Updated 19.08.2023
dynEGA_errors <- function(
    data, id, group, n.embed, tau, delta,
    use.derivatives, ncores, verbose
)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "dynEGA")
  
  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }
  
  # 'id' errors
  if(!is.null(id)){
    length_error(id, 1, "dynEGA")
    typeof_error(id, c("numeric", "character"), "dynEGA")
  }
  # 'group' errors
  if(!is.null(group)){
    length_error(group, 1, "dynEGA")
    typeof_error(group, c("numeric", "character"), "dynEGA")
  }
  
  # 'n.embed' errors
  length_error(n.embed, 1, "dynEGA")
  typeof_error(n.embed, "numeric", "dynEGA")
  range_error(n.embed, c(3, Inf), "dynEGA")
  
  # 'tau' errors
  length_error(tau, 1, "dynEGA")
  typeof_error(tau, "numeric", "dynEGA")
  range_error(tau, c(1, 3), "dynEGA") # don't allow more than 3
  
  # 'delta' errors
  length_error(delta, 1, "dynEGA")
  typeof_error(delta, "numeric", "dynEGA")
  range_error(delta, c(1, 3), "dynEGA") # don't allow more than 3
  
  # 'use.derivatives' errors
  length_error(use.derivatives, 1, "dynEGA")
  typeof_error(use.derivatives, "numeric", "dynEGA")
  range_error(use.derivatives, c(0, 2), "dynEGA") # only 0, 1, or 2
  
  # 'ncores' errors
  length_error(ncores, 1)
  typeof_error(ncores, "numeric", "dynEGA")
  range_error(ncores, c(1, parallel::detectCores()), "dynEGA")
  
  # 'verbose' errors
  length_error(verbose, 1, "dynEGA")
  typeof_error(verbose, "logical", "dynEGA")
  
  # Return data in case of tibble
  return(data)
  
}

#' @exportS3Method 
# S3 Print Method (General) ----
# Updated 15.08.2023
print.dynEGA <- function(x, ...)
{
  
  # Print general dynEGA
  cat(
    styletext(
      text = styletext(
        text =  "Dynamic Exploratory Graph Analysis\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Get methods attributes
  glla_attributes <- attr(x, "glla")
  
  # Print GLLA information at the top _only_
  cat("Number of Embeddings: ", glla_attributes$n.embed)
  cat("\nEmbedding Offset (tau): ", glla_attributes$tau)
  cat("\nTime between Observations (delta/lag): ", glla_attributes$delta)
  cat(
    "\nDerivatives: ", switch(
      as.character(glla_attributes$use.derivatives),
      "0" = "Moving Average (0)",
      "1" = "First-order (1)",
      "2" = "Second-order (2)"
    )
  )
  
  # Add full breakspace for each level
  cat(
    paste0(
      "\n\n",
      paste0(rep("-", options("width")), collapse = ""),
      "\n\n"
    )
  )
  
  # Get `EGA` objects
  ega_objects <- get_EGA_object(x)
  
  # Determine NULLs
  null_objects <- !lvapply(ega_objects, is.null)
  
  # Print population first
  if(null_objects["population"]){
    print(ega_objects$population)
  }
  
  # Print group second
  if(null_objects["group"]){
    
    # Check for breakspace
    if(null_objects["population"]){
      cat(
        paste0(
          "\n\n",
          paste0(rep("-", options("width")), collapse = ""),
          "\n\n"
        )
      )
    }
    
    # Print group
    print(ega_objects$group)

  }
  
  # Print individuals third
  if(null_objects["individual"]){
    
    # Check for breakspace
    if(null_objects["population"] | null_objects["group"]){
      cat(
        paste0(
          "\n\n",
          paste0(rep("-", options("width")), collapse = ""),
          "\n\n"
        )
      )
    }
    
    # Print individual
    print(ega_objects$individual)
    
  }
  
}

#' @exportS3Method 
# S3 Print Method (Population) ----
# Updated 07.07.2023
print.dynEGA.Population <- function(x, ...)
{
  
  # Print population
  cat(
    styletext(
      text = styletext(
        text =  "Population\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Switch class to "EGA"
  class(x) <- "EGA"
  
  # Borrow print from `EGA`
  print(x)
  
}

#' @exportS3Method 
# S3 Print Method (Group) ----
# Updated 07.07.2023
print.dynEGA.Group <- function(x, ...)
{
  
  # Print group
  cat(
    styletext(
      text = styletext(
        text =  "Group\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Get number of groups
  group_length <- length(x)
  
  # Loop over and print
  for(group in seq_len(group_length)){
    
    # Print group name
    cat("Group: ", names(x)[group], "\n\n")
    
    # Switch class to "EGA"
    class(x[[group]]) <- "EGA"
    
    # Borrow print from `EGA`
    print(x[[group]])
    
    # Add breakspace if not last group
    if(group != group_length){
      cat("\n\n------------\n\n")
    }
    
  }
  
}

#' @exportS3Method 
# S3 Print Method (Individual) ----
# Updated 11.10.2023
print.dynEGA.Individual <- function(x, ...)
{
  
  # Print individual
  cat(
    styletext(
      text = styletext(
        text =  "Individual\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Follow `bootEGA` print strategy
  
  # Print network information
  send_network_methods(x[[1]]$network, boot = TRUE)
  
  # Add line break
  cat("\n")
  
  # Get unidimensional attributes
  unidimensional_attributes <- attr(x[[1]], "unidimensional")
  
  # Obtain unidimensional method
  unidimensional_method <- switch(
    unidimensional_attributes$uni.method,
    "expand" = "Expand",
    "le" = "Leading Eigenvector",
    "louvain" = "Louvain"
  )
  
  # Set up unidimensional print
  if(
    unidimensional_method == "Louvain" &
    "consensus.iter" %in% names(unidimensional_attributes$consensus)
  ){
    
    # Set up consensus attributes
    consensus_attributes <- unidimensional_attributes$consensus
    
    # Obtain consensus name
    consensus_name <- switch(
      consensus_attributes$consensus.method,
      "highest_modularity" = "Highest Modularity",
      "iterative" = "Iterative",
      "most_common" = "Most Common",
      "lowest_tefi" = "Lowest TEFI"
    )
    
    # Update unidimensional method text
    unidimensional_method <- paste0(
      unidimensional_method, " (", consensus_name,
      " for ", consensus_attributes$consensus.iter,
      " iterations)"
    )
    
  }
  
  # Print unidimensional
  cat("Unidimensional Method: ", unidimensional_method)
  
  # Add breakspace
  cat("\n\n----\n\n")
  
  # Get communities
  communities <- nvapply(x, function(x){unique_length(x$wc)})
  
  # Set up like `bootEGA`
  
  # Print number of cases
  cat("Number of cases: ", length(x), "\n")
  
  # Print median
  cat(
    paste0(
      "\nMedian dimensions: ", median(communities, na.rm = TRUE)
    )
  )
  
  # Add space
  cat("\n")
  
  # Get frequencies
  frequencies <- fast_table(communities)
  
  # Get frequency table
  frequency_df <- as.data.frame(
    t(data.frame(
      names(frequencies),
      frequencies
    ))
  )
  
  # Adjust dimension names (`dimnames` doesn't work)
  colnames(frequency_df) <- NULL
  row.names(frequency_df) <- c("", "Frequency: ")
  # Finally, print
  print(frequency_df)
  
}

#' @exportS3Method 
# S3 Summary Method (General) ----
# Updated 07.07.2023
summary.dynEGA <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Summary Method (Population) ----
# Updated 07.07.2023
summary.dynEGA.Population <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Summary Method (Group) ----
# Updated 07.07.2023
summary.dynEGA.Group <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Summary Method (Individual) ----
# Updated 07.07.2023
summary.dynEGA.Individual <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Plot Method (General) ----
# Updated 31.07.2023
plot.dynEGA <- function(x, base = 1, id = NULL, ...)
{
  
  # Determine non-NULLs
  non_null_objects <- !lvapply(get_EGA_object(x), is.null)
  
  # Get number of NULLs
  null_total <- sum(!non_null_objects)
  
  # Plot population first
  if(non_null_objects["population"]){
    
    # Send plot
    plot(x$dynEGA$population, ...)
    
    # Print message only if other `dynEGA` objects
    if(null_total != 2){
      message("dynEGA Population")
    }
    
  }
  
  # Print group second
  if(non_null_objects["group"]){
    
    # Check for breakspace
    if(non_null_objects["population"]){
      
      # Allow user to proceed at their own pace
      sink <- readline("Press <ENTER> for 'Group' plot")
      
    }
    
    # Plot group
    plot(x$dynEGA$group, base = base, ...)
    
    # Print message only if other `dynEGA` objects
    if(null_total != 2){
      message("dynEGA Group")
    }
    
  }
  
  # Print individuals third
  if(non_null_objects["individual"]){
    
    # Check for breakspace
    if(non_null_objects["population"] | non_null_objects["group"]){
      # Allow use to proceed at their own pace
      sink <- readline("Press <ENTER> for 'Individual' plot")
    }
    
    # Print individual
    plot(x$dynEGA$individual, id = id, ...)
    
    # Print message only if other `dynEGA` objects
    if(null_total != 2){
      message("dynEGA Individual")
    }
    
  }
  
}

#' @exportS3Method 
# S3 Plot Method (Population) ----
# Updated 07.07.2023
plot.dynEGA.Population <- function(x, ...)
{
  silent_plot(
    single_plot(
      network = x$network,
      wc = x$wc,
      ...
    )
  )
}

#' @exportS3Method 
# S3 Plot Method (Group) ----
# Updated 07.07.2023
plot.dynEGA.Group <- function(x, base = 1, ...)
{
  
  # Get names
  group_names <- names(x)
  
  # Convert base group name to number (if necessary)
  if(is.character(base)){
    base <- which(group_names == base)
  }
  
  # Order group names
  group_names <- c(group_names[base], group_names[-base])
  
  # Remove base from groups
  base_object <- x[[base]]
  
  # Extract other groups
  other_objects <- x[-base]
  
  # Get sequence length of other objects
  sequence_length <- seq_len(length(other_objects))
  
  # Get base plot
  base_plot <- silent_load(
    single_plot(
      network = base_object$network,
      wc = base_object$wc,
      arguments = TRUE,
      ...
    )
  )
  
  # Set up comparison plots
  comparison_plots <- lapply(
    sequence_length, function(i){
      compare_plots(
        comparison_network = other_objects[[i]]$network,
        comparison_wc = other_objects[[i]]$wc,
        plot_ARGS = base_plot$ARGS
      )
    }
  )
  
  # Converge into single list
  plotlist <- c(
    base = list(base_plot$network_plot),
    comparison_plots[sequence_length]
  )
  
  # Set up for comparison
  silent_plot(
    ggpubr::ggarrange(
      plotlist = plotlist,
      labels = group_names,
      legend = "bottom",
      ...
    )
  )
  
}

#' @exportS3Method 
# S3 Plot Method (Individual) ----
# Updated 30.07.2023
plot.dynEGA.Individual <- function(x, base = 1, id = NULL, ...)
{
  
  # Get ellipse
  ellipse <- list(...)
  
  # Check for ID
  if(!is.null(id)){
    
    # Convert IDs to numeric (if necessary)
    if(is.character(id)){
      id <- which(names(x) == id)
    }
    
    # Check for multiple IDs
    if(length(id) != 1){ # Perform similar operation to groups
      
      # Get names
      ID_names <- names(x)
      
      # Order group names
      ID_names <- c(ID_names[id[base]], ID_names[id[-base]])
      
      # Use first ID as base
      base_object <- x[[id[base]]]
      
      # Extract other IDs
      other_objects <- x[id[-base]]
      
      # Get sequence length of other objects
      sequence_length <- seq_len(length(other_objects))
      
      # Get base plot
      base_plot <- silent_load(
        do.call(
          what = basic_plot_setup,
          args = c(
            list(
              network = base_object$network,
              wc = base_object$wc,
              arguments = TRUE
            ), ellipse
          )
        )
      )
  
      # Set up comparison plots
      comparison_plots <- lapply(
        sequence_length, function(i){
          compare_plots(
            comparison_network = other_objects[[i]]$network,
            comparison_wc = other_objects[[i]]$wc,
            plot_ARGS = base_plot$ARGS
          )
        }
      )
      
      # Converge into single list
      plotlist <- c(
        base = list(base_plot$network_plot),
        comparison_plots[sequence_length]
      )
      
      # `ggarrange` does not like non-plot arguments for its ellipse
      ellipse <- ellipse[
        names(ellipse) %in% names(formals(ggpubr::ggarrange))
      ]
      
      # Set up for comparison
      silent_plot(
        do.call(
          what = ggpubr::ggarrange,
          args = c(
            list(
              plotlist = plotlist,
              labels = ID_names,
              legend = "bottom"
            ), ellipse
          )
        )
      )
      
    }else{
      
      # If only one ID, then plot it
      silent_plot(
        single_plot(
          network = x[[id]]$network,
          wc = x[[id]]$wc,
          ...
        )
      )
      
    }
    
  }else{ # No ID provided, then randomly some plots
    
    # Get number of individuals
    ID_length <- length(x)
    
    # Ensure there are at least four; otherwise, all of them
    ID_numbers <- min(ID_length, 4)
    
    # Randomly sample from IDs
    random_IDs <- shuffle(seq_len(ID_length), size = ID_numbers)
    
    # Get ID names
    ID_names <- names(x)[random_IDs]
    
    # Remove base from individuals
    base_object <- x[[random_IDs[1]]]
    
    # Extract other groups
    other_objects <- x[random_IDs[-1]]
    
    # Get sequence length of other objects
    sequence_length <- seq_len(length(other_objects))
    
    # Get base plot
    base_plot <- silent_load(
      do.call(
        what = basic_plot_setup,
        args = c(
          list(
            network = base_object$network,
            wc = base_object$wc,
            arguments = TRUE
          ), ellipse
        )
      )
    )
    
    # Set up comparison plots
    comparison_plots <- lapply(
      sequence_length, function(i){
        compare_plots(
          comparison_network = other_objects[[i]]$network,
          comparison_wc = other_objects[[i]]$wc,
          plot_ARGS = base_plot$ARGS
        )
      }
    )
    
    # Converge into single list
    plotlist <- c(
      base = list(base_plot$network_plot),
      comparison_plots[sequence_length]
    )
    
    # `ggarrange` does not like non-plot arguments for its ellipse
    ellipse <- ellipse[
      names(ellipse) %in% names(formals(ggpubr::ggarrange))
    ]
    
    # Set up for comparison
    silent_plot(
      do.call(
        what = ggpubr::ggarrange,
        args = c(
          list(
            plotlist = plotlist,
            labels = ID_names,
            legend = "bottom"
          ), ellipse
        )
      )
    )
    
  }
  
}

#' @noRd
# Get ID from data ----
# Updated 28.09.2023
get_ID <- function(data, id, level, variable_names, dimensions)
{
  
  # First, check for 'id' in variable names (revert)
  if("id" %in% variable_names){
    
    # Get 'id' from data if already listed as a column
    id <- which(variable_names == "id")
    
    # Get 'id'
    ID <- data[,id]
    
    # Remove 'id' from data
    data <- data[,-id]
    
  }else if(is.null(id)){
    
    # Set id to be `1` (same as "population")
    ID <- rep(1, dim(data)[1])
    
    # Send warning
    if("individual" %in% level){
      warning(
        "'level' included \"individual\" but no 'id' was provided. Setting all 'id' to `1` or same as level = \"population\"",
        call. = FALSE
      )
    }
    
  }else{ # Otherwise, check for proper IDs
    
    # Check first for length
    id_length <- length(id)
    
    # If not a column name or number, then vector
    if(id_length != 1){ 
      
      # Check for errors
      object_error(id, "vector", "dynEGA")
      length_error(id, dimensions[1], "dynEGA")
      
      # Make the vector into ID
      ID <- id
      
    }else{
      
      # Regardless of column name or number,
      # it can be extracted directly from the data
      ID <- data[,id]
      
      # If 'id' is column name, then turn it to number
      if(is.character(id)){
        id <- which(variable_names == tolower(id))
      }
      
      # Remove 'id' from data
      data <- data[,-id]
      
    }
    
  }
  
  # Return list
  return(list(data = data, ID = ID))
  
}

#' @noRd
# Get Group from data ----
# Updated 28.09.2023
get_Group <- function(data, group, level, variable_names, dimensions)
{
  
  # First, check for 'group' in variable names (revert)
  if("group" %in% variable_names){
    
    # Get 'group' from data if already listed as a column
    group <- which(variable_names == "group")
    
    # Get 'group'
    Group <- data[,group]
    
    # Remove 'group' from data
    data <- data[,-group]
    
  }else if(is.null(group)){
    
    # Set group to be `1` (same as "population")
    Group <- rep(1, dim(data)[1])
    
    # Send warning
    if("group" %in% level){
      warning(
        "'level' included \"group\" but no 'id' was provided. Setting all 'group' to `1` or same as level = \"population\"",
        call. = FALSE
      )
    }
    
  }else{ # Otherwise, check for proper Groups
    
    # Check first for length
    group_length <- length(group)
    
    # If not a column name or number, then vector
    if(group_length != 1){ 
      
      # Check for errors
      object_error(group, "vector", "dynEGA")
      length_error(group, dimensions[1], "dynEGA")
      
      # Make the vector into ID
      Group <- group
      
    }else{
      
      # Regardless of column name or number,
      # it can be extracted directly from the data
      Group <- data[,group]
      
      # If 'group' is column name, then turn it to number
      if(is.character(group)){
        group <- which(variable_names == tolower(group))
      }
      
      # Remove 'group' from data
      data <- data[,-group]
      
    }
    
  }
  
  # Return list
  return(list(data = data, Group = Group))
  
}

#' @noRd
# Get "ID" and "Group" attributes ----
# Updated 07.07.2023
get_attributes <- function(data, dimensions, id, group, level)
{
  
  # Check for proper ID
  # If "individual" in 'level', then assigns "ID" as attribute
  # If 'id' is in the data, then it is removed!
  ID_result <- get_ID(data, id, level, tolower(dimnames(data)[[2]]), dimensions)
  
  # Check for proper Group
  # If "group" in 'level', then assigns "Group" as attribute
  # If 'group' is in the data, then it is removed!
  Group_result <- get_Group(
    ID_result$data, group, level,
    tolower(dimnames(ID_result$data)[[2]]),
    dimensions
  )
  
  # Get data from Group result
  data <- Group_result$data
  
  # Return data with attributes assigned
  return(
    structure(
      data, ID = ID_result$ID, Group = Group_result$Group
    )
  )
  
}

#' @noRd
# Variable derivatives ----
# Updated 07.07.2023
variable_derivatives <- function(variable, n.embed, tau, L)
{
  
  # Get derivatives
  derivatives <- Embed(variable, n.embed, tau) %*% L
  
  # Add column names
  dimnames(derivatives)[[2]] <- paste0("Ord", 0:2)
  
  # Return derivatives
  return(derivatives)
  
}

#' @noRd
# Individual derivatives ----
# Updated 09.07.2023
individual_derivatives <- function(
    individual_data, variable_names,
    n.embed, tau, L,
    individual_attributes
)
{
  
  # Get derivatives for each participant
  participant_derivatives <- lapply(
    seq_along(individual_data), function(index){
      
      # Apply over variables
      derivatives <- do.call(
        cbind, lapply(
          as.data.frame(individual_data[[index]]),
          variable_derivatives, n.embed, tau, L
        )
      )
      
      # Add names
      dimnames(derivatives)[[2]] <- paste(
        rep(variable_names, each = 3),
        dimnames(derivatives)[[2]], sep = "."
      )
      
      # Get ID
      ID <- names(individual_data)[index]
      
      # Get Group
      Group <- unique(
        individual_attributes$Group[individual_attributes$ID == ID]
      )
      
      # Get length of time series derivatives\
      ts_length <- dim(derivatives)[1]
      
      # Return derivatives with updated attributes
      return(
        structure(
          derivatives,
          ID = rep(ID, ts_length),
          Group = rep(Group, ts_length)
        )
      )
      
    }
  )
  
  # Add IDs to derivatives
  names(participant_derivatives) <- unique(individual_attributes$ID)
  
  # Return derivatives
  return(participant_derivatives)
  
}

#' @noRd
# Organize individual `dynEGA` results ----
# Updated 09.07.2023
process_individual_dynEGA <- function(usable_derivatives, individual_results)
{
  
  # Regardless, add IDs to individual results
  names(individual_results) <- cvapply(usable_derivatives, attr, "ID")
  
  # Determine whether there are any IDs with zero variance derivatives
  zero_variance_ID <- lvapply(usable_derivatives, function(x){any(attr(x, "zero_variance"))})
  
  # Check for any zero variance IDs
  if(any(zero_variance_ID)){
    
    # Send warning about which IDs
    warning(
      paste(
        paste0("IDs: ", paste0(names(individual_results)[zero_variance_ID], collapse = ", ")),
        "\nhad derivatives with zero variance (no change in values across time). These IDs will have disconnected nodes in their network and missing community memberships"
      ), call. = FALSE
    )
    
    # Pre-compute node names
    node_names <- names(attr(usable_derivatives[[1]], "zero_variance"))
    
    # As well as number of nodes
    nodes <- length(node_names)
    
    # Adjust networks, memberships, and dimension outputs
    individual_results[zero_variance_ID] <- lapply(
      which(zero_variance_ID), function(ID){
        
        # Get zero variance attribute
        zero_variance <- attr(usable_derivatives[[ID]], "zero_variance")
        
        # Initialize new network
        new_network <- matrix(
          NA, nrow = nodes, ncol = nodes,
          dimnames = list(node_names, node_names)
        )
        
        # Initialize new correlations with same attributes
        new_correlations <- new_network
        
        # Add network into new network
        new_network[!zero_variance, !zero_variance] <- 
          individual_results[[ID]]$network
        
        # Restore attributes
        attributes(new_network) <- c(
          attributes(new_network),
          attributes(individual_results[[ID]]$network)[c("methods", "class")]
        )
        
        # Return new network
        individual_results[[ID]]$network <- new_network
        
        # Initialize new memberships
        new_memberships <- rep(NA, nodes)
        
        # Add names
        names(new_memberships) <- node_names
        
        # Add memberships into new memberships
        new_memberships[!zero_variance] <- individual_results[[ID]]$wc
        
        # Restore attributes
        attributes(new_memberships) <- c(
          attributes(new_memberships),
          attributes(individual_results[[ID]]$wc)[c("methods", "class")]
        )
        
        # Return new memberships
        individual_results[[ID]]$wc <- new_memberships
        
        # Add correlation into new correlations
        new_correlations[!zero_variance, !zero_variance] <- 
          individual_results[[ID]]$correlation
        
        # Return new correlations (no attributes to pass)
        individual_results[[ID]]$correlation <- new_correlations
        
        # Initialize new dimension output
        dim.variables <- fast.data.frame(
          c(node_names, as.vector(new_memberships)),
          nrow = nodes, ncol = 2,
          colnames = c("items", "dimension")
        )
        
        # Dimension variables data frame by dimension
        individual_results[[ID]]$dim.variables <- dim.variables[
          order(dim.variables$dimension),
        ]
        
        # At long last... return
        return(individual_results[[ID]])
        
      }
    )
    
  }
  
  # Return individual results
  return(individual_results)
  
}
