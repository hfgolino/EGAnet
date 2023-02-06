#%%%%%%%%%%%%%%%%%%%%%%%%
# FUNCTION FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Standard EGA arguments
# Updated 02.02.2023
ega_arguments <- function(arguments)
{
  
  # Set default arguments
  default_arguments <- list(
    data = NULL, n = NULL,
    corr = "cor_auto", uni.method = "louvain",
    model = "glasso", model.args = list(),
    algorithm = "walktrap", algorithm.args = list(),
    consensus.method = "most_common_tefi", consensus.iter = 100,
    plot.EGA = TRUE, plot.args = list()
  )
  
  # Match arguments with input
  matched_arguments <- match(names(arguments), names(default_arguments))
  
  # Check for whether any are not missing
  if(any(!is.na(matched_arguments))){
    
    # Remove missing
    replace_arguments <- na.omit(matched_arguments)
    
    # Replace arguments
    default_arguments[matched_arguments] <- arguments
    
  }
  
  # Return arguments
  return(default_arguments)
  
}

#' @noRd
# Standard {lavaan} arguments
# Updated 02.02.2023
lavaan_arguments <- function(arguments)
{
  
  # Set default arguments
  default_arguments <- list(
    model = NULL, data = NULL, ordered = NULL,
    sampling.weights = NULL, sample.cov = NULL,
    sample.mean = NULL, sample.th = NULL,
    sample.nobs = NULL, group = NULL, cluster = NULL,
    constraints = "", WLS.V = NULL, NACOV = NULL,
    std.lv = TRUE
  )
  
  # Match arguments with input
  matched_arguments <- match(names(arguments), names(default_arguments))
  
  # Check for whether any are not missing
  if(any(!is.na(matched_arguments))){
    
    # Remove missing
    replace_arguments <- na.omit(matched_arguments)
    
    # Replace arguments
    default_arguments[matched_arguments] <- arguments
    
  }
  
  # Return arguments
  return(default_arguments)
  
}

#' @noRd
# Make unidimensional CFA model
# Updated 02.02.2023
make_unidimensional_cfa <- function(variable_names)
{
 
  # Check for number of variables
  if(length(variable_names) == 2){
    
    # Make latent factor model
    # Fix all loadings to 1
    model <- paste0(
      "LF =~ ",
      paste0(
        "1*", variable_names, collapse = " + "
      )
    )
    
    # Fix residual variances to same value
    model <- paste(
      model,
      paste0(
        variable_names,
        "~~theta*",
        variable_names,
        collapse = " \n "
      ),
      sep = " \n "
    )
    
  }else{
    
    # Make latent factor model
    model <- paste0(
      "LF =~ ",
      paste(
        variable_names, collapse = " + "
      )
    )
    
  }
  
  # Return model
  return(model)
  
}

#' @noRd
# Determine estimator arguments
# Updated 02.02.2023
estimator_arguments <- function(lavaan_ARGS)
{
  
  # Obtain categories
  categories <- data_categories(
    data = lavaan_ARGS$data
  )
  
  # Check for categories
  if(any(categories < 6)){
    
    # Set arguments
    lavaan_ARGS$estimator <- "WLSMV"
    lavaan_ARGS$missing <- "pairwise"
    lavaan_ARGS$ordered <- names(categories)[
      categories < 6
    ]
    
  }else{
    
    # Set arguments
    lavaan_ARGS$estimator <- "MLR"
    lavaan_ARGS$missing <- "fiml"
    lavaan_ARGS$ordered <- FALSE
    
  }
  
  # Return arguments
  return(lavaan_ARGS)
  
  
}




