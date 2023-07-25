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
#' a simple structure.
#' Defaults to \code{\link[GPArotation]{geominQ}}.
#' For a list of rotations, see \code{\link{GPArotation}}
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
# Network Scores ----
# Updated 15.07.2023
net.scores <- function (
    data, A, wc, rotation = NULL,
    loading.method = c("BRM", "experimental"),
    scoring.method = c(
      "Anderson", "Bartlett", "components",
      "Harman", "network", "tenBerge", "Thurstone"
    ), impute = c("mean", "median", "none"),
    ...
)
{
  
  # Error on:
  # 1. missing data
  # 2. symmetric matrix input as data
  if(missing(data)){ # Check for missing data
    stop(
      "Input for 'data' is missing. To compute scores, the original data are necessary",
      call. = FALSE
    )
  }else if(is_symmetric(data)){ # Check for symmetric matrix
      stop(
        "Input for 'data' was detected as symmetric. Correlation matrices cannot be used. The original data are required to estimate scores.",
        call. = FALSE
      )
  }
  
  # Ensure data is a matrix
  data <- as.matrix(data) 
  
  # Get ellipse arguments
  ellipse <- list(...)
  
  # Check for defunct "method"
  if("method" %in% names(ellipse)){
    scoring.method <- ellipse$method
  }
  
  # Check for missing arguments (argument, default, function)
  loading.method <- set_default(loading.method, "brm", net.loads)
  scoring.method <- set_default(scoring.method, "network", net.scores)
  impute <- set_default(impute, "none", net.scores)
  
  # Perform imputation
  if(impute != "none"){
    data <- imputation(data, impute)
  }
  
  # Compute network loadings (will handle `EGA` objects)
  loadings <- net.loads(
    A = A, wc = wc, loading.method = loading.method, 
    rotation = rotation, ...
  )
  
  # Return results
  return(
    list(
      scores = compute_scores(loadings, data, scoring.method),
      loadings = loadings
    )
  )
  
}

#' @noRd
# Imputation ----
# Updated 14.07.2023
imputation <- function(data, impute)
{
  
  # Get imputation function
  impute_values <- switch(
    impute,
    "mean" = colMeans(data, na.rm = TRUE),
    "median" = nvapply(as.data.frame(data), median, na.rm = TRUE)
  )
  
  # Loop over unique columns
  for(column in unique(missing_data[,"col"])){
    
    # Obtain target rows
    target_rows <- missing_data[,"row"][missing_data[,"col"] == column]
    
    # Populate data
    data[target_rows, column] <- variable_impute[column]
    
  }
  
  # Return data
  return(data)
  
}

#' @noRd
# Zero-out cross-loadings ----
# Consistent with hierarchical CFA
# Updated 25.07.2023
zero_out <- function(loadings, wc){
  
  # Get names of memberships
  wc_names <- names(wc)
  
  # Get loadings names
  loadings_names <- dimnames(loadings)[[2]]
  
  # Loop over unique memberships
  for(membership in unique(wc)){

    # Set cross-loadings to zero
    loadings[
      wc_names[wc == membership],
      loadings_names != membership
    ] <- 0
    
  }
  
  # Return loadings
  return(loadings)
  
}

#' @noRd
# Scores computation ----
# Updated 25.07.2023
compute_scores <- function(loadings, data, scoring.method)
{
  
  # Method must exist, so continue
  if(scoring.method == "network"){
    
    # Compute unrotated scores
    unrotated <- network_scores(
      data = data,
      loads = zero_out(
        loadings$std, attr(loadings, "membership")$wc
      )
    )
    
    # Compute rotated scores (if available)
    if(!is.null(loadings$rotated)){
      rotated <- network_scores(
        data = data,
        loads = zero_out(
          loadings$rotated$loadings,
          attr(loadings, "membership")$wc
        )
      )
    }else{rotated <- NULL}
    
  }else{
    
    # Switch lowercase method back to appropriate case
    scoring.method <- switch(
      scoring.method,
      "anderson" = "Anderson", 
      "bartlett" = "Bartlett",
      "components" = "components",
      "harman" = "Harman", 
      "tenberge" = "tenBerge",
      "thurstone" = "Thurstone"
    )
    
    # Compute unrotated scores
    unrotated <- psych::factor.scores(
      x = data,
      f = loadings$std,
      method = scoring.method
    )$scores
    
    # Compute rotated scores (if available)
    if(!is.null(loadings$rotated)){
      rotated <- psych::factor.scores(
        x = data,
        f = loadings$rotated$loadings,
        Phi = loadings$rotated$Phi,
        method = scoring.method
      )$scores
    }else{rotated <- NULL}
    
  }
  
  # Return scores
  return(
    list(
      std.scores = unrotated,
      rot.scores = rotated
    )
  )
  
}

#' @noRd
# Network scores computation ----
# Updated 24.07.2023
network_scores <- function(data, loads)
{
  return(scale(data) %*% loads[dimnames(data)[[2]],])
}
