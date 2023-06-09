#' Automatic correlations
#'
#' This wrapper is similar to \code{\link[qgraph]{cor_auto}}. There
#' are some minor adjustments that make this function simpler and to
#' function within \code{\link{EGAnet}} as desired. \code{NA} values are not treated
#' as categories (this behavior differs from \code{\link[qgraph]{cor_auto}})
#'
#' @param data Matrix or data frame.
#' Should consist only of variables that are desired to be correlated
#' 
#' @param method Character (length = 1).
#' The standard correlation method to be used.
#' Defaults to \code{"pearson"}.
#' Using \code{"pearson"} will compute polychoric, tetrachoric, polyserial,
#' and biserial correlations for categorical and categorical/continuous correlations
#' by default. To obtain \code{"pearson"} correlations regardless, use \code{\link{cor}}.
#' Other options of \code{"kendall"} and \code{"spearman"} are provided for
#' completeness and use \code{\link{cor}} 
#'
#' @param ordinal.categories Numeric (length = 1).
#' \emph{Up to} the number of categories \emph{before} a variable is considered continuous.
#' Defaults to \code{7} categories before \code{8} is considered continuous
#' 
#' @param forcePD Boolean (length = 1).
#' Whether positive definite matrix should be enforced.
#' Defaults to \code{TRUE}
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
#' @param empty.method Character (length = 1).
#' Method for empty cell correction in \code{\link[EGAnet]{polychoric.matrix}}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"none"}}
#' {Adds no value (\code{empty.value = "none"})
#' to the empirical joint frequency table between two variables}
#' 
#' \item{\code{"zero"}}
#' {Adds \code{empty.value} to the cells with zero
#' in the joint frequency table between two variables}
#' 
#' \item{\code{"all"}}
#' {Adds \code{empty.value} to all
#' in the joint frequency table between two variables}
#' 
#' }
#' 
#' @param empty.value Character (length = 1).
#' Value to add to the joint frequency table cells in \code{\link[EGAnet]{polychoric.matrix}}.
#' Accepts numeric values between 0 and 1 or
#' specific methods:
#' 
#' \itemize{
#' 
#' \item{\code{"none"}}
#' {Adds no value (\code{0}) to the empirical joint
#' frequency table between two variables}
#' 
#' \item{\code{"point_five"}}
#' {Adds \code{0.5} to the cells defined by \code{empty.method}}
#' 
#' \item{\code{"one_over"}}
#' {Adds \code{1 / n} where \code{n} equals the number of cells
#' based on \code{empty.method}. For \code{empty.method = "zero"},
#' }
#' 
#' }
#' 
#' @param verbose Boolean (length = 1).
#' Whether messages should be printed.
#' Defaults to \code{FALSE}
#' 
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' # Obtain correlations
#' wmt_corr <- auto.correlate(wmt)
#'
#' @export
#'
# Automatic correlations ----
# Updated 09.06.2023
auto.correlate <- function(
    data, # Matrix or data frame
    method = c("kendall", "pearson", "spearman"), # allow changes to standard correlations
    ordinal.categories = 7, # consider ordinal up to 7 categories
    forcePD = TRUE, # ensure result is positive definite
    na.data = c("pairwise", "listwise"), # use available or complete values
    empty.method = c("none", "zero", "all"), # zero frequencies in categorical correlations
    empty.value = c("none", "point_five", "one_over"), # value to use in zero cells
    verbose = FALSE # don't print messages
)
{
  
  # Missing arguments
  ## Standard correlation method
  if(missing(method)){
    method <- "pearson"
  }else{method <- tolower(match.arg(method))}
  ## Missing data
  if(missing(na.data)){
    na.data <- "pairwise"
  }else{na.data <- tolower(match.arg(na.data))}
  ## Empty cell method
  if(missing(empty.method)){
    empty.method <- "none"
  }else{empty.method <- tolower(match.arg(empty.method))}
  ## Empty cell value
  if(missing(empty.value)){
    empty.value <- "none"
  }
  
  # Convert to matrix
  data <- as.matrix(data)
  
  # Set names
  data <- ensure_dimension_names(data)
  
  # Assumes some preprocessing has taken place
  # to only select for appropriate variables
  
  # Determine whether categorical correlations are necessary
  if(method != "pearson"){
    
    # Obtain correlation matrix
    correlation_matrix <- cor(
      x = data, use = na.data,
      method = method
    )
    
  }else{ # Proceed with determination of categorical correlations
    
    # Obtain the number of categories for each variables
    # See `helpers-general.R`
    categories <- data_categories(data)
  
    # Determine categorical variables
    categorical_variables <- which(categories <= ordinal.categories)
    
    # Determine number of categorical variables
    categorical_number <- length(categorical_variables)
    
    # Determine whether there are any categorical variables
    if(categorical_number != 0){
      
      # Determine continuous variables
      continuous_variables <- setdiff(
        1:ncol(data), categorical_variables
      )
      
      # Set up correlation matrix
      correlation_matrix <- matrix(
        0, nrow = ncol(data),
        ncol = ncol(data)
      )
      
      # Set names
      colnames(correlation_matrix) <- colnames(data)
      row.names(correlation_matrix) <- colnames(data)
      
      # Determine whether there are more than one categorical variables
      if(categorical_number > 1){
        
        # Compute categorical correlations (only correlation for categorical data)
        # Add correlations to correlation matrix
        correlation_matrix[
          categorical_variables, categorical_variables # ensure proper indexing
        ] <- polychoric.matrix(
          data = data[,categorical_variables], na.data = na.data,
          empty.method = empty.method, empty.value = empty.value
        )
        
      }
      
      # Determine whether there are more than one continuous variables
      if(length(continuous_variables) > 1){
        
        # Compute continuous correlations (only correlation for continuous data)
        # Add correlations to correlation matrix
        correlation_matrix[
          continuous_variables, continuous_variables # ensure proper indexing
        ] <- cor(
          x = data[,continuous_variables],
          use = na.data, method = method
        )
        
      }
      
      # Determine whether there are mixed variables
      if(categorical_number != ncol(data)){ # Check for mixed variables
        
        # Loop over categorical indices
        for(i in categorical_variables){
          
          # Polyserial correlations based on {polycor}
          mixed_correlations <- polyserial.vector(
            categorical_variable = data[,i],
            continuous_variables = data[,continuous_variables],
            na.data = na.data
          )
          
          # Compute categorical/continuous correlations
          correlation_matrix[i, continuous_variables] <- mixed_correlations
          
          # Fill other side
          correlation_matrix[continuous_variables, i] <- mixed_correlations
          
        }
        
      }
    
    }else{
      
      # Compute Pearson's correlations
      correlation_matrix <- cor(
        x = data, use = na.data,
        method = method
      )
      
    }

  }
  
  # Determine whether matrix is positive definite
  if(isTRUE(forcePD) & !all(eigen(correlation_matrix)$values > 0)){
    
    # Send warning to user (if `verbose`)
    if(isTRUE(verbose)){
      warning("Correlation matrix is not positive definite. Finding nearest positive definite matrix using `Matrix::nearPD`")
    }
    
    # Regardless, make matrix positive definite
    correlation_matrix <- as.matrix(
      Matrix::nearPD(
        correlation_matrix, corr = TRUE,
        ensureSymmetry = TRUE, keepDiag = TRUE
      )$mat
    )
    
  }
  
  # Return correlation matrix
  return(correlation_matrix)

}

# Bug checking ----
# ## Different categories
# set.seed(1234)
# data = latentFactoR::simulate_factors(
#   factors = 4, variables = 4,
#   loadings = 0.55, cross_loadings = 0.05,
#   correlations = 0.30, sample_size = 1000,
#   variable_categories = c(
#     rep(2, 4), rep(5, 4),
#     rep(7, 4), rep(Inf, 4)
#   )
# )$data;
# ordinal.categories = 7;
# method = "pearson"; forcePD = TRUE;
# na.data = "pairwise"; empty.method = "none";
# empty.value = "none"; verbose = FALSE;
# 
# # Compare against {qgraph}'s `cor_auto`
# qgraph_correlations <- qgraph::cor_auto(data)
# EGAnet_correlations <- auto.correlate(data)
# 
# # Difference
# max(abs(EGAnet_correlations - qgraph_correlations))
# # Biggest difference is between polyserial (7 categories with continuous)
# 
# ## Add missing data
# data[sample(1:length(data), 1000)] <- NA
# # Compare against {qgraph}'s `cor_auto`
# qgraph_correlations <- qgraph::cor_auto(
#   data, 
#   ordinalLevelMax = 8
#   # Needs to have 8 levels to account for missing data!!
# )
# EGAnet_correlations <- auto.correlate(data)
# 
# # Difference
# max(abs(EGAnet_correlations - qgraph_correlations))
# # Biggest difference is between polyserial (7 categories with continuous)
# 
## Zero cell counts
# data <- cbind(
#   c(0, 1, 2, 3, 4, 1, 2, 3, 1),
#   c(0, 1, 2, 2, 1, 2, 2, 4, 2)
# ); ordinal.categories = 7;
# method = "pearson"; forcePD = TRUE;
# na.data = "pairwise"; empty.method = "none";
# empty.value = "none"; verbose = FALSE;
#
# The above bug checks for categorical data
# have been verified to match the output
# of Turbofuns::PolychoricRM to a maximum
# difference less than or equal to 1.0e-06
# (or one step beyond floating point)

# Compute polyserial correlation ----
#' For a single categorical variable, compute correlations
#' with \emph{n} continuous variables
#' 
#' Uses two-step approximation from {polycor}'s \code{polyserial}
#' 
#' @noRd
# Updated 09.06.2023
polyserial.vector <- function(
    categorical_variable, continuous_variables,
    na.data = c("pairwise", "listwise")
)
{
  
  # Ensure matrices (see `helpers-general.R` for more details)
  categorical_variable <- force_matrix(categorical_variable)
  continuous_variables <- force_matrix(continuous_variables)
  
  # Determine cases based on `na.data` argument
  if(na.data == "pairwise"){
   
    # Pairwise cases
    categorical_cases <- apply(
      continuous_variables, 2, function(x){
        
        # Combine into single matrix
        combined_cases <- cbind(categorical_variable, x)
        
        # Remove cases and compute rows
        nrow(na.omit(combined_cases))
        
      }
    )
     
  }else if(na.data == "listwise"){
    
    # Combine into single matrix
    combined_cases <- cbind(categorical_variable, continuous_variables)
    
    # Remove cases and compute rows
    categorical_cases <- rep(
      nrow(na.omit(combined_cases)), ncol(continuous_variables)
    )
    
  }

  # Scale continuous data
  scaled_continuous <- scale(continuous_variables)
  
  # Compute thresholds
  thresholds <- obtain_thresholds(categorical_variable)
  
  # Compute correlation
  correlation <- sqrt((categorical_cases - 1) / categorical_cases) * 
    sd(categorical_variable, na.rm = TRUE) * 
    cor(categorical_variable, scaled_continuous, use = na.data) /
    sum(dnorm(thresholds))
  
  # Return correlation
  return(correlation)
  
  
}

# Compute thresholds ----
#' @noRd
# Updated 09.06.2023
obtain_thresholds <- function(categorical_variable)
{
  
  # Obtain table
  frequency <- table(categorical_variable)
  
  # Obtain cumulative sums
  cumulative_sum <- cumsum(frequency)
  
  # Obtain thresholds
  thresholds <- qnorm(
    cumulative_sum[-length(cumulative_sum)] / 
      cumulative_sum[length(cumulative_sum)]
  )
  
  # Return thresholds
  return(thresholds)
  
}












