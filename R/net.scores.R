#' Network Scores
#'
#' @description This function computes network scores computed based on
#' each node's \code{strength} within each
#' community (i.e., factor) in the network (see \code{\link[EGAnet]{net.loads}}).
#' These values are used as network "factor loadings" for the weights of each item.
#' Notably, network analysis allows nodes to contribution to more than one community.
#' These loadings are considered in the network scores. In addition,
#' if the construct is a hierarchy (e.g., personality questionnaire;
#' items in facet scales in a trait domain), then an overall
#' score can be computed (see argument \code{global}). An important difference
#' is that the network scores account for cross-loadings in their
#' estimation of scores
#'
#' @param data Matrix or data frame.
#' Must be a dataset
#'
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object.
#' An adjacency matrix of network data
#'
#' @param wc Numeric.
#' A vector of community assignments.
#' Not necessary if an \code{\link[EGAnet]{EGA}} object
#' is input for argument \code{A}
#'
#' @param rotation Character.
#' A rotation to use, like factor loadings, to obtain
#' a simple structure. For a list of rotations,
#' see \link{GPArotation}
#' 
#' @param method Character.
#' Factor scoring method to use. For a list of scoring methods,
#' see \code{\link[psych]{factor.scores}}.
#' Defaults to \code{"network"}
#' 
#' "network" scores can also be computed as a formative method
#' using \code{"component"} or \code{"network"}. The core difference
#' between these two methods is that \code{"network"} will use the
#' weight contributions to communities using each variable's
#' standard deviation in combination with their loadings
#'
#' @param impute Character.
#' In the presence of missing data, imputation can be implemented. Currently,
#' three options are available:
#' 
#' @param ... Additional arguments.
#' Arguments to be passed onto \code{\link[EGAnet]{net.loads}}
#'
#' \itemize{
#'
#' \item{\strong{\code{none}}}
#' {No imputation is performed. This is the default.}
#'
#' \item{\strong{\code{mean}}}
#' {The "mean" value of the columns are used to replace the missing data.}
#'
#' \item{\strong{\code{median}}}
#' {The "median" value of the columns are used to replace the missing data.}
#'
#' }
#'
#' @return Returns a list containing:
#'
#' \item{unstd.scores}{The unstandardized network scores for each participant
#' and community (including the overall score)}
#'
#' \item{std.scores}{The standardized network scores for each participant
#' and community (including the overall score)}
#'
#' \item{commCor}{Partial correlations between the specified or identified communities}
#'
#' \item{loads}{Standardized network loadings for each item in each dimension
#' (computed using \code{\link[EGAnet]{net.loads}})}
#'
#' @details For more details, type \code{vignette("Network_Scores")}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )}
#'
#' # Network scores
#' net.scores(data = wmt, A = ega.wmt)
#' 
#' \dontrun{
#' # Produce Methods section
#' methods.section(
#'   ega.wmt,
#'   stats = "net.scores"
#' )}
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}, \emph{53}, 1563-1580.
#' 
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}, 1095-1108.
#'
#' Golino, H., Christensen, A. P., Moulder, R., Kim, S., & Boker, S. M. (2021).
#' Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
#' \emph{Psychometrika}.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Scores
# Updated: 13.04.2023
# Add rotation: 20.10.2022
net.scores <- function (
    data, A, wc, rotation = "oblimin",
    method = "network", impute = "none",
    ...
)
{
  
  # Missing arguments handling
  if(missing(data)){
    stop("Argument 'data' is required for analysis")
  }
  
  # Set network and memberships
  if(is(A, "EGA")){
    wc <- A$wc
    A <- A$network
  }else if(is(A, "dynEGA")){
    wc <- A$dynEGA$wc
    A <- A$dynEGA$network
  }else if(missing(A)){
    stop("Adjacency matrix is required for analysis")
  }else if(missing(wc)){
    wc <- rep(1, ncol(data))
  }
  
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Perform imputation
  if(impute != "none"){
    data <- imputation(data = data, impute = impute)
  }
  
  # Compute network loadings
  loadings <- net.loads(A = A, wc = wc, ...)
  
  # Check for method to compute scores
  score_results <- compute_scores(
    loadings_object = loadings,
    data = data, method = method,
    wc = wc
  )
  
  # Set up results list
  results <- list(
    scores = list(
      std.scores = as.data.frame(scale(score_results$unrotated)),
      rot.scores = as.data.frame(scale(score_results$rotated))
    ),
    loadings = loadings
  )

  # Class
  class(results) <- "NetScores"
  
  return(results)
}

# Imputation ----
#' @noRd
# Function to handle imputation
imputation <- function(data, impute)
{
  
  # Determine missing data
  missing_data <- which(is.na(data), arr.ind = TRUE)
  
  # Check for mean or median imputation
  if(impute == "mean"){
    
    # Compute means
    variable_impute <- colMeans(data, na.rm = TRUE)
    
  }else if(impute == "median"){
    
    # Compute medians
    variable_impute <- apply(data, 2, median, na.rm = TRUE)
    
  }
  
  # Determine unique columns
  unique_columns <- unique(missing_data[,"col"])
  
  # Loop over unique columns
  for(column in unique_columns){
    
    # Obtain target column
    target_column <- missing_data[,"col"] == column
    
    # Obtain target rows
    target_rows <- missing_data[,"row"][target_column]
    
    # Populate data
    data[target_rows, column] <- variable_impute[column]
    
  }
  
  # Return data
  return(data)
  
}

# Network scores computation ----
#' @noRd
# Function to compute network scores
network_scores <- function(loads, data, wc)
{
  
  # Initialize matrix for network scores
  scores <- matrix(
    NA, nrow = nrow(data),
    ncol = ncol(loads)
  )
  
  # Reorder data to match loadings
  data <- data[,row.names(loads)]
  
  # # Determine signs and whether coding is necessary
  # signs <- numeric(nrow(loads))
  # names(signs) <- row.names(loads)
  # 
  # # Obtain unique memberships
  # unique_wc <- na.omit(unique(wc))
  # 
  # # Loop over to obtain direct signs
  # for(current_wc in unique_wc){
  # 
  #   # Obtain signs for current community
  #   target_loadings <- loads[wc == current_wc, as.character(current_wc)]
  # 
  #   # Obtain signs
  #   signs[names(target_loadings)] <- sign(target_loadings)
  # 
  # }
  # 
  # # Reverse data (if necessary)
  # if(any(signs == -1)){
  # 
  #   # Flip loadings
  #   loads <- loads * signs
  # 
  #   # Loop over all variables
  #   for(i in 1:ncol(data)){
  # 
  #     # Only do something with negative signs
  #     if(signs[i] == -1){
  # 
  #       # Obtain number of categories
  #       target_categories <- length(na.omit(unique(data[,i])))
  # 
  #       # Determine if data are categorical
  #       if(target_categories <= 7){
  # 
  #         # Reverse code data
  #         data[,i] <- (max(data[,i], na.rm = TRUE) + 1) - data[,i]
  # 
  #       }else{
  # 
  #         # Flip signs of data
  #         data[,i] <- -data[,i]
  # 
  #       }
  # 
  #     }
  # 
  #   }
  # 
  # }
  
  
  # Loop over communities
  for(i in 1:ncol(loads)){
    
    # Obtain target loadings
    target_loadings <- loads[,i]
    
    # Identify which loadings which are not zero
    non_zero <- target_loadings != 0
    
    # Obtain data for non-zero loadings
    non_zero_loadings <- target_loadings[non_zero]
    non_zero_data <- data[,non_zero]
    
    # Obtain standard deviations
    standard_devs <- apply(non_zero_data, 2, sd, na.rm = TRUE)
    
    # Obtain relative weight
    relative <- non_zero_loadings / standard_devs
    relative_weight <- relative / sum(relative, na.rm = TRUE)
    
    # Multiple by data
    score <- as.vector( # Ensure vector
      colSums(t(non_zero_data) * relative_weight, na.rm = FALSE)
    )
    
    # Add to matrix
    scores[,i] <- score
  
  }
  
  # Add column names
  colnames(scores) <- colnames(loads)
  
  # Return scores
  return(scores)
  
}

# Scores computation ----
#' @noRd
# Wrapper to compute scores
compute_scores <- function(loadings_object, data, method, wc)
{
  
  # Set methods
  methods <- c(
    "Thurstone", "tenBerge", "Anderson",
    "Bartlett", "Harman", "components",
    "network"
  )
  
  # Check for method
  if(!method %in% methods){
    
    # Set up methods
    methods_text <- paste0("\"", methods, "\"", collapse = ", ")
    
    # Send error
    stop(
      paste0(
        "Method \"", method, "\" in not available. \n",
        "Current options include: ",
        methods_text
      )
    )
    
  }
  
  # Method must exist, so continue
  if(method == "network"){
    
    # Compute unrotated scores
    unrotated <- network_scores(
      loads = loadings_object$std,
      data = data, wc = wc
    )
    
    # Compute rotated scores
    rotated <- network_scores(
      loads = loadings_object$rotated$loadings,
      data = data, wc = wc
    )
    
  }else{
    
    # Compute unrotated scores
    unrotated <- psych::factor.scores(
      x = data,
      f = loadings_object$std,
      method = method
    )$scores
    
    # Compute rotated scores
    rotated <- psych::factor.scores(
      x = data,
      f = loadings_object$rotated$loadings,
      Phi = loadings_object$rotated$Phi,
      method = method
    )$scores
    
  }
  
  # Set up results list
  results <- list(
    unrotated = unrotated,
    rotated = rotated
  )
  
  # Return scores
  return(results)
  
}
