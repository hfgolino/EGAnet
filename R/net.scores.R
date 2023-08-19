#' @title Network Scores
#'
#' @description This function computes network scores computed based on
#' each node's \code{strength} within each community in the network 
#' (see \code{\link[EGAnet]{net.loads}}). These values are used as "network loadings" 
#' for the weights of each variable.
#' 
#' Network scores are computed as a formative composite rather than a reflective factor.
#' This composite representation is consistent with no latent factors that psychometric
#' network theory proposes.
#' 
#' Scores can be computed as a "simple" structure, which is equivalent to a weighted
#' sum scores or as a "full" structure, which is equivalent to an EFA approach.
#' Conservatively, the "simple" structure approach is recommended until further
#' validation
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param A Network matrix, data frame, or \code{\link[EGAnet]{EGA}} object
#'
#' @param wc Numeric or character vector (length = \code{ncol(A)}).
#' A vector of community assignments.
#' If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
#' then \code{wc} is automatically detected
#'
#' @param loading.method Character (length = 1).
#' Sets network loading calculation based on implementation
#' described in \code{"BRM"} (Christensen & Golino, 2021) or
#' an \code{"experimental"} implementation.
#' Defaults to \code{"BRM"}
#' 
#' @param rotation Character.
#' A rotation to use to obtain a simpler structure. 
#' For a list of rotations, see \code{\link[GPArotation]{rotations}} for options.
#' Defaults to \code{NULL} or no rotation.
#' By setting a rotation, \code{scores} estimation will be
#' based on the rotated loadings rather than unrotated loadings
#' 
#' @param scores Character (length = 1).
#' How should scores be estimated?
#' Defaults to \code{"network"} for network scores.
#' Set to other scoring methods which will be computed using
#' \code{\link[psych]{factor.scores}} (see link for arguments
#' and explanations for other methods)
#' 
#' @param loading.structure Character (length = 1).
#' Whether simple structure or the saturated loading matrix
#' should be used when computing scores.
#' Defaults to \code{"simple"}
#' 
#' \code{"simple"} structure more closely mirrors sum scores and CFA; 
#' \code{"full"} structure more closely mirrors EFA
#' 
#' Simple structure is the more "conservative" (established) approach
#' and is therefore the default. Treat \code{"full"} as experimental
#' as proper vetting and validation has not been established
#'
#' @param impute Character (length = 1).
#' If there are any missing data, then imputation can be implemented. 
#' Available options:
#'
#' \itemize{
#'
#' \item{\code{"none"} --- }
#' {Default. No imputation is performed}
#'
#' \item{\code{"mean"} --- }
#' {The mean value of each variable is used to replace missing data
#' for that variable}
#'
#' \item{\code{"median"} --- }
#' {The median value of each variable is used to replace missing data
#' for that variable}
#'
#' }
#' 
#' @param ... Additional arguments to be passed on to 
#' \code{\link[EGAnet]{net.loads}} and
#' \code{\link[psych]{factor.scores}}
#'
#' @return Returns a list containing:
#'
#' \item{scores}{A list containing the standardized (\code{std.scores})
#' rotated (\code{rot.scores}) scores. If \code{rotation = NULL}, then
#' \code{rot.scores} will be \code{NULL}}
#' 
#' \item{loadings}{Output from \code{\link[EGAnet]{net.loads}}}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' # Network scores
#' net.scores(data = wmt, A = ega.wmt)
#'
#' @references
#' \strong{Original implementation and simulation for loadings} \cr
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}, \emph{53}, 1563-1580.
#' 
#' \strong{Preliminary simulation for scores} \cr
#' Golino, H., Christensen, A. P., Moulder, R., Kim, S., & Boker, S. M. (2021).
#' Modeling latent topics in social media using Dynamic Exploratory Graph Analysis: The case of the right-wing and left-wing trolls in the 2016 US elections.
#' \emph{Psychometrika}.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Scores ----
# Updated 19.08.2023
net.scores <- function (
    data, A, wc, 
    loading.method = c("BRM", "experimental"),
    rotation = NULL,
    scores = c(
      "Anderson", "Bartlett", "components",
      "Harman", "network", "tenBerge", "Thurstone"
    ),
    loading.structure = c("simple", "full"),
    impute = c("mean", "median", "none"),
    ...
)
{
  
  # All argument errors are handled here or in `net.loads`
  
  # Error on:
  # 1. missing data
  # 2. symmetric matrix input as data
  if(missing(data)){ # Check for missing data
    .handleSimpleError(
      h = stop,
      msg = "Input for 'data' is missing. To compute scores, the original data are necessary", 
      call = "net.scores"
    )
  }else if(is_symmetric(data)){ # Check for symmetric matrix
    .handleSimpleError(
      h = stop,
      msg = "Input for 'data' was detected as symmetric. Correlation matrices cannot be used. The original data are required to estimate scores.",
      call = "net.scores"
    )
  }
  
  # Ensure data is a matrix
  data <- as.matrix(usable_data(data, verbose = TRUE)) 
  
  # Get ellipse arguments
  ellipse <- list(...)
  
  # Check for defunct "method"
  if("method" %in% names(ellipse)){
    scores <- ellipse$method
  }
  
  # Check for missing arguments (argument, default, function)
  loading.method <- set_default(loading.method, "brm", net.loads)
  scores <- set_default(scores, "network", net.scores)
  loading.structure <- set_default(loading.structure, "simple", net.scores)
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
      scores = compute_scores(
        loadings, data, scores, loading.structure
      ),
      loadings = loadings
    )
  )
  
}

#' @noRd
# Imputation ----
# Updated 04.08.2023
imputation <- function(data, impute)
{
  
  # Get imputation function
  impute_values <- switch(
    impute,
    "mean" = colMeans(data, na.rm = TRUE),
    "median" = nvapply(as.data.frame(data), median, na.rm = TRUE)
  )
  
  # Get missing data
  missing_data <- which(is.na(data), arr.ind = TRUE)
  
  # Loop over unique columns
  for(column in unique(missing_data[,"col"])){
    
    # Obtain target rows
    target_rows <- missing_data[,"row"][missing_data[,"col"] == column]
    
    # Populate data
    data[target_rows, column] <- impute_values[column]
    
  }
  
  # Return data
  return(data)
  
}

#' @noRd
# Zero-out cross-loadings ----
# Consistent with hierarchical CFA
# Updated 28.07.2023
zero_out <- function(loadings, wc, loading.structure){
  
  # Check for loading structure
  if(loading.structure == "simple"){
    
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
    
  }
  
  # Return loadings
  return(loadings)
  
}

#' @noRd
# Scores computation ----
# Updated 04.08.2023
compute_scores <- function(
    loadings, data, scores, loading.structure
)
{
  
  # Method must exist, so continue
  if(scores == "network"){
    
    # Compute unrotated scores
    unrotated <- network_scores(
      data = data,
      loads = zero_out(
        loadings$std, attr(loadings, "membership")$wc,
        loading.structure
      )
    )
    
    # Compute rotated scores (if available)
    if(!is.null(loadings$rotated)){
      rotated <- network_scores(
        data = data,
        loads = zero_out(
          loadings$rotated$loadings,
          attr(loadings, "membership")$wc,
          loading.structure
        )
      )
    }else{rotated <- NULL}
    
  }else{
    
    # Switch lowercase method back to appropriate case
    scores <- switch(
      scores,
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
      method = scores
    )$scores
    
    # Compute rotated scores (if available)
    if(!is.null(loadings$rotated)){
      rotated <- psych::factor.scores(
        x = data,
        f = loadings$rotated$loadings,
        Phi = loadings$rotated$Phi,
        method = scores
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
