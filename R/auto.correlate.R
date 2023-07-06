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
#' @param corr Character (length = 1).
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
#' @param ...
#' Not actually used but makes it either for general functionality
#' in the package
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
# Updated 03.07.2023
auto.correlate <- function(
    data, # Matrix or data frame
    corr = c("kendall", "pearson", "spearman"), # allow changes to standard correlations
    ordinal.categories = 7, # consider ordinal up to 7 categories
    forcePD = TRUE, # ensure result is positive definite
    na.data = c("pairwise", "listwise"), # use available or complete values
    empty.method = c("none", "zero", "all"), # zero frequencies in categorical correlations
    empty.value = c("none", "point_five", "one_over"), # value to use in zero cells
    verbose = FALSE, # don't print messages
    ... # not actually used
)
{
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "pearson", auto.correlate)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  empty.method <- set_default(empty.method, "none", auto.correlate)
  empty.value <- set_default(empty.value, "none", auto.correlate)
  
  # Convert to matrix
  data <- as.matrix(data)
  
  # Set names
  data <- ensure_dimension_names(data)
  
  # Get variable names
  variable_names <- dimnames(data)[[2]]
  
  # Set dimensions
  dimensions <- dim(data)
  
  # Assumes some preprocessing has taken place
  # to only select for appropriate variables
  
  # Determine whether categorical correlations are necessary
  if(corr != "pearson"){
    
    # Obtain correlation matrix
    correlation_matrix <- cor(
      x = data, use = na.data,
      method = corr
    )
    
  }else{ # Proceed with determination of categorical correlations
    
    # Obtain the number of categories for each variables
    categories <- data_categories(data)
  
    # Determine categorical variables
    categorical_variables <- categories <= ordinal.categories
    
    # Determine number of categorical variables
    categorical_number <- sum(categorical_variables)
    
    # Determine whether there are any categorical variables
    if(categorical_number != 0){
      
      # Determine continuous variables
      continuous_variables <- !categorical_variables
      
      # Set up correlation matrix
      correlation_matrix <- matrix(
        nrow = dimensions[2], ncol = dimensions[2],
        dimnames = list(variable_names, variable_names)
      )
      
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
      if(sum(continuous_variables) > 1){
        
        # Compute continuous correlations (only correlation for continuous data)
        # Add correlations to correlation matrix
        correlation_matrix[
          continuous_variables, continuous_variables # ensure proper indexing
        ] <- cor(
          x = data[,continuous_variables],
          use = na.data, method = corr
        )
        
      }
      
      # Determine whether there are mixed variables
      if(categorical_number != dimensions[2]){ # Check for mixed variables
        
        # Obtain continuous data (keep as matrix)
        continuous_data <- data[,continuous_variables, drop = FALSE]
        
        # Loop over categorical indices
        for(i in which(categorical_variables)){
        
          # Fill matrix
          correlation_matrix[continuous_variables, i] <-
            correlation_matrix[i, continuous_variables] <-
            polyserial.vector( # computes polyserial vector
              categorical_variable = data[,i], # drops to vector
              continuous_variables = continuous_data,
              na.data = na.data
            )
          
        }
        
      }
    
    }else{
      
      # Compute Pearson's correlations
      correlation_matrix <- cor(
        x = data, use = na.data,
        method = corr
      )
      
    }

  }
  
  # Determine whether matrix is positive definite
  if(isTRUE(forcePD) & !is_positive_definite(correlation_matrix)){
    
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
# corr = "pearson"; forcePD = TRUE;
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
# corr = "pearson"; forcePD = TRUE;
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
# Updated 03.07.2023
polyserial.vector <- function(
    categorical_variable, continuous_variables,
    na.data = c("pairwise", "listwise")
)
{
  
  # Determine cases based on `na.data` argument
  if(na.data == "pairwise"){
   
    # Pairwise cases
    categorical_cases <- colSums(
      !is.na(categorical_variable) * 
      !is.na(continuous_variables)
    )
    
  }else if(na.data == "listwise"){
    
    # Complete cases
    complete_cases <- complete.cases(cbind(categorical_variable, continuous_variables))
    
    # Repeat cases for number of continuous variables
    categorical_cases <- rep(sum(complete_cases), dim(continuous_variables)[2])
    
  }

  # Compute correlation
  return(
    sqrt((categorical_cases - 1) / categorical_cases) * 
      sd(categorical_variable, na.rm = TRUE) * 
      # Correlations with scaled continuous variables
      cor(categorical_variable, scale(continuous_variables), use = na.data) /
      # Compute sum of thresholds
      sum(dnorm(obtain_thresholds(categorical_variable)), na.rm = TRUE)
  )
  
  
}

# Compute thresholds ----
#' @noRd
# Updated 03.07.2023
obtain_thresholds <- function(categorical_variable)
{
  
  # Obtain cumulative sums from frequency table
  cumulative_sum <- cumsum(table(categorical_variable))
  
  # Obtain cumulative length
  cumsum_length <- length(cumulative_sum)
  
  # Obtain thresholds
  return(
    qnorm(
      cumulative_sum[-cumsum_length] / 
      cumulative_sum[cumsum_length]
    )
  )
  
}












