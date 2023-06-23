#%%%%%%%%%%%%%%%%%%%%%%%%#
#### {EGAnet} Helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%%#

# These helpers are internal functions that are repeatedly
# used across the package. They are implemented here once
# so that changes only need to be made here and not in
# multiple places across the package.
#
# The goal is to minimize redundant calls to improve reliability of code.
# Sections are separated by how they are associated with their utility
# in other functions.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# REPRODUCIBLE RANDOM GENERATION ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# For these random generation functions, the goal
# is reproducibility, which runs counter to "random".
# However, for bootstrap results to be reproduced
# there is a need to generate "random" samples that
# can be reproduced so that the user arrives at the
# same results *without* influencing R's random number
# generation or changing a user-defined seed previous
# to, while, or after using {EGAnet} functions.
#
# The `reproducible_sample` and `reproducible_seed`
# functions were written in C++ with the assistance
# of GPT-4 

#' @noRd
# Random normal wrapper ----
# About 4.5x faster than R's `rnorm`
# Updated 17.06.2023
rnorm_ziggurat <- function(n, mean = 0, sd = 1)
{
  return(mean + sd * RcppZiggurat::zrnorm(n))
}

#' @noRd
# Generate multivariate normal data ----
# Removes MASS package dependency (from version 7.3.54)
# Function is streamlined to avoid chained `if` calls in original function
# Updated 17.06.2023
MASS_mvrnorm <- function(
    n = 1, mu, Sigma,
    tol = 1e-06, empirical = FALSE, EISPACK = FALSE
)
{
  
  # Get number of variables
  p <- length(mu)
  
  # Assumes `mu` == `Sigma` dimensions
  # if (!all(dim(Sigma) == c(p, p)))
  
  # EISPACK is always FALSE
  # if (EISPACK)
  
  # Obtain eigenvectors and values
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values # Obtain values
  
  # SKIP positive definite check (should be from `auto.correlate`)
  # if (!all(ev >= -tol * abs(ev[1L])))
  
  # Generate data
  X <- matrix(rnorm_ziggurat(n = p * n), n)
  # X <- matrix(rnorm(p * n), n) # replaces `rnorm` so we can set seed
  
  # For empirical data
  if(isTRUE(empirical)){
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  
  # Obtain X
  X <- mu + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  
  # Skip checks
  # nm <- names(mu)
  # if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
  #   nm <- dn[[1L]]
  # dimnames(X) <- list(nm, NULL)
  
  # Return results
  if(n == 1){
    return(drop(X))
  }else{
    return(t(X))
  }
  
}

#' @noRd
# Generate multivariate normal data (quicker) ----
# Pre-computes `p`
# Avoids repeated calls to `eigen` in `MASS_mvrnorm`
# Pre-computes `eS$vectors %*% diag(sqrt(pmax(ev, 0)), p)`
# Updated 17.06.2023
MASS_mvrnorm_quick <- function(
    n = 1, mu, p, coV,
    tol = 1e-06, empirical = FALSE, EISPACK = FALSE
)
{
  
  # Get number of variables
  # p <- length(mu)
  
  # Assumes `mu` == `Sigma` dimensions
  # if (!all(dim(Sigma) == c(p, p)))
  
  # EISPACK is always FALSE
  # if (EISPACK)
  
  # Obtain eigenvectors and values
  # eS <- eigen(Sigma, symmetric = TRUE)
  # ev <- eS$values # Obtain values
  
  # SKIP positive definite check (should be from `auto.correlate`)
  # if (!all(ev >= -tol * abs(ev[1L])))
  
  # Generate data
  X <- matrix(rnorm_ziggurat(n = p * n), n)
  # X <- matrix(rnorm(p * n), n) # replaces `rnorm` so we can set seed
  
  # For empirical data
  if(isTRUE(empirical)){
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  
  # Obtain X
  # X <- mu + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  X <- mu + coV %*% t(X)
  
  # Skip checks
  # nm <- names(mu)
  # if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
  #   nm <- dn[[1L]]
  # dimnames(X) <- list(nm, NULL)
  
  # Return results
  if(n == 1){
    return(drop(X))
  }else{
    return(t(X))
  }
  
}

#' @noRd
# Generate reproducible parametric bootstrap ----
# (multivariate normal samples)
# A wrapper over `rnorm_ziggurat` and `MASS_mvrnorm`
# Updated 18.06.2023
reproducible_parametric <- function(
    samples, cases, mu, Sigma, seed
)
{
  
  # First, set seed for random normal generation
  RcppZiggurat::zsetseed(seed)
  # Seed does *not* affect R's seed and RNG
  
  # Next, perform pre-computations
  p <- length(mu) # avoids repeated calls to `length`
  eS <- eigen(Sigma, symmetric = TRUE) # avoids repeated calls to `eigen`
  coV <- eS$vectors %*% diag(sqrt(pmax(eS$values, 0)), p)
  # avoids repeated matrix computation
  
  # Then, generate samples
  return(
    lapply(
      seq_len(samples), MASS_mvrnorm_quick,
      n = cases, mu = mu,
      p = p, coV = coV
    )
  )
  
}

#' @noRd
# Generate reproducible seed ----
# Seed generation that does *not* affect R's RNG
# Get seeds from zero up to 32-bit maximum
# `0` is reserved for random seed
# Updated 17.06.2023
reproducible_seed <- function(n, seed = NULL)
{
  return(r_sample_seeds(n, ifelse(is.null(seed), 0, seed)))
}

#' @noRd
# Generate reproducible shuffling ----
# About 6x faster than R's `sample`
# Updated 17.06.2023
reproducible_sample <- function(
    x, size, replace, seed
){
  
  # Check for with replacement
  if(isTRUE(replace)){ # With replacement
    
    # Get indices
    ## Implemented in C++
    ## Internal function
    shuffled_indices <- r_sample_with_replacement(size, seed)
    
  }else{ # Without replacement
    
    # Get indices
    ## Implemented in C++
    ## Internal function
    shuffled_indices <- r_sample_without_replacement(seq_len(size), seed)
    
  }
  
  # Set shuffled data
  shuffled <- x[shuffled_indices]
  
  # Return shuffled data
  return(shuffled)
  
}

#' @noRd
# Generate reproducible resampling bootstrap ----
# (multivariate normal samples)
# A wrapper over `reproducible_seed` and `reproducible_sample`
# Updated 18.06.2023
reproducible_resampling <- function(
    data, samples, cases, seed
)
{
  
  # First, generate as many seeds as there are samples
  seeds <- reproducible_seed(n = samples, seed = seed)
  
  # Then, generate samples
  ## More direct approach than calling `reproducible_sample`
  return(
    lapply(
      seq_len(samples),
      function(i){
        data[
          r_sample_with_replacement(n = cases, seed = seeds[i]),
        ]
      }
    )
  )
  
}

#' @noRd
# Generate reproducible bootstrap data ----
# Wrapper for `reproducible_parametric` and
# `reproducible_resampling`
# Updated 18.06.2023
reproducible_bootstrap <- function(
    data, samples, cases, mu, Sigma, seed,
    type = c("parametric", "resampling")
)
{
  
  # Based on bootstrap type, generate data
  if(type == "parametric"){
    
    # Obtain bootstrap samples
    bootstrap_samples <- reproducible_parametric(
      samples = samples, cases = cases,
      mu = mu, Sigma = Sigma, seed = seed
    )
    
    
  }else if(type == "resampling"){
    
    # Obtain bootstrap samples
    bootstrap_samples <- reproducible_resampling(
      data = data, samples = samples,
      cases = cases, seed = seed
    )
    
  }
  
  # Return bootstrap samples
  return(bootstrap_samples)
  
}

#%%%%%%%%%%%%%%%%%%%%%%
# FASTER SEQUENCES ----
#%%%%%%%%%%%%%%%%%%%%%%

# I haven't looked into why `seq_len(n)` is faster
# than `1:n` or `dim(data)[1]` is faster than
# `nrow(data)` (same for columns) but both are
# faster than the latter commonly applied applications
#
# GPT-4 reports that `nrow` and `ncol` call `dim`

#' @noRd
# Faster row sequences ----
# Updated 13.06.2023
nrow_sequence <- function(data)
{
  return(seq_len(dim(data)[1]))
}


#' @noRd
# Faster column sequences ----
# Updated 13.06.2023
ncol_sequence <- function(data)
{
  return(seq_len(dim(data)[2]))
}

# # Evidence:
# #
# # Load benchmark
# library(microbenchmark)
# 
# # Create large matrix
# large_matrix <- matrix(0, nrow = 10000, ncol = 10000)
# 
# # Run row test (around 30 nanoseconds)
# microbenchmark(
#   1:nrow(large_matrix),
#   nrow_sequence(large_matrix),
#   times = 10000L,
#   control = list(warmup = 1000L)
# )
# 
# # Run column test (around 30 nanoseconds)
# microbenchmark(
#   1:ncol(large_matrix),
#   ncol_sequence(large_matrix),
#   times = 10000L,
#   control = list(warmup = 1000L)
# )
#
# Does it matter? Probably not...
#
# ¯\_(ツ)_/¯
#
# An actual real reason to use `seq_len`: handles edge cases
# when there is a length = 0

#%%%%%%%%%%%%%%%%%%%%%%%
# GENERAL FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Determine number of digits in a number ----
# Updated 27.05.2023
digits <- function(number)
{
  # Obtain the lowest value of log base 10 and add 1
  return(floor(log10(number)) + 1)
  
}

#' @noRd
# Force vector ----
# (usually for data frame rows/columns)
# Updated 27.05.2023
force_vector <- function(desired_vector)
{
  # Convert to matrix then vector
  new_vector <- as.vector(as.matrix(desired_vector))
  
  # Keep names (if any)
  if(!is.null(colnames(desired_vector))){
    names(new_vector) <- colnames(desired_vector)
  }
  
  return(new_vector)
  
}

#' @noRd
# Force matrix ----
# (usually for vectors)
# Updated 12.06.2023
force_matrix <- function(desired_matrix, dimension = c("col", "row"))
{
  
  # Check for missing dimension argument
  if(missing(dimension)){
    dimension <- "col"
  }else{dimension <- tolower(match.arg(dimension))}
  
  # Check for matrix form already
  if(is(desired_matrix, "matrix")){
    return(desired_matrix)
  }else if(is(desired_matrix, "data.frame")){
    return(as.matrix(desired_matrix))
  }else{# Convert vector to matrix
    
    # Set up as a single column or row
    if(dimension == "col"){
      return(matrix(desired_matrix, ncol = 1))
    }else{
      return(matrix(desired_matrix, nrow = 1))
    }
    
  }
  
}

#' @noRd
# Force numeric ----
# (usually for vectors)
# Updated 09.06.2023
force_numeric <- function(desired_numeric)
{
  
  # Check first for feasible coercion
  if(canCoerce(desired_numeric, "numeric")){
    
    # Check for edge cases
    if(is.factor(desired_numeric) | is.ordered(desired_numeric)){ # factor/ordered
      
      # Make character
      desired_numeric <- as.character(desired_numeric)
      
    }else if(is.complex(desired_numeric)){ # complex (force `NA``)
      
      # Find imaginary numbers (start as character)
      imaginary_character <- as.character(desired_numeric)
      
      # Determine which are imaginary
      imaginary <- grepl("i", imaginary_character)
      
      # Set imaginary numbers to `NA`
      desired_numeric[imaginary] <- NA
      
    }
    
  }else{
    
    # Return all NAs
    desired_numeric <- rep(NA, length(desired_numeric))
    
  }
  
  # Make numeric
  desired_numeric <- as.numeric(desired_numeric)
  
  # Return result
  return(desired_numeric)
  
}

#' @noRd
# Edge count ----
# Counts the number of edges (assumes network is matrix)
# Updated 18.06.2023
edge_count <- function(network, nodes, diagonal = FALSE)
{
  
  # Count edges
  return(
    (sum(network != 0) - ifelse(diagonal, nodes, 0)) / 2
  )
  
}

#' @noRd
# Count table ----
# Provides counts of repeated rows in a data frame
# (clever solution by GPT-4)
# Updated 16.06.2023
count_table <- function(data, proportion = FALSE)
{
  
  # Make data frame
  data <- as.data.frame(data)
  
  # Obtain counts of unique rows (clever solution from GPT-4)
  if(isTRUE(proportion)){
    counts <- table(do.call(paste, data)) / dim(data)[1]
  }else{
    counts <- table(do.call(paste, data))
  }
  
  # Prepare a data frame
  count_df <- as.data.frame(
    do.call(rbind, lapply(strsplit(names(counts), split = " "), as.numeric)),
    stringsAsFactors = FALSE
  )
  
  # Attach counts
  count_df$count <- counts
  
  # Ensure proper naming
  colnames(count_df) <- c(colnames(data), ifelse(proportion, "Proportion", "Count"))
  
  # Return data frame
  return(count_df)
  
}

#' @noRd
# Determine number of categories in data ----
# Updated 02.02.2023
data_categories <- function(data)
{
  
  # Ensure data is matrix
  data <- as.matrix(data)
  
  # Loop over columns
  categories <- apply(
    data, 2, function(x){
      length(na.omit(unique(x)))
    }
  )
  
  # Return categories
  return(categories)
  
}

#' @noRd
# Convert version to number ----
# Updated 02.02.2023
version_conversion <- function(version)
{
  
  # Convert to character
  version <- as.character(version)
  
  # Remove periods
  version <- gsub("\\.", "", version)
  
  # Convert to numeric
  version <- as.numeric(version)
  
  # Return version
  return(version)
  
}

#' @noRd
# All-purpose symmetric checker ----
# Updated 03.02.2023
is_symmetric <- function(data){
  
  # Check for whether rows equal columns
  if(nrow(data) == ncol(data)){
    
    # Convert to matrix
    data_matrix <- as.matrix(data)
    
    # Remove names
    data_matrix <- unname(data_matrix)
    
    # Check that lower triangle equal upper triangle
    lower_triangle <- data_matrix[lower.tri(data_matrix)]
    transpose_matrix <- t(data_matrix) # ensures similar orientation
    upper_triangle <- transpose_matrix[lower.tri(transpose_matrix)]
    
    # Check that all are equal
    all_equal <- all(lower_triangle == upper_triangle, na.rm = TRUE)
    
  }else{
    
    # Not a matrix
    return(FALSE)
    
  }
  
  # Return whether all are equal
  return(all_equal)
  
}

#' @noRd
# Format number with certain decimals ----
# Mainly for naming and printing
# Updated 14.06.2024
format_decimal <- function(numbers, places)
{
  
  return(
    formatC(
      x = numbers,
      digits = places,
      format = "f", flag = "0"
    )
  )

}

#' @noRd
# Format number with certain integer ----
# Mainly for naming and printing
# Updated 14.06.2024
format_integer <- function(numbers, places)
{
  
  return(
    formatC(
      x = numbers,
      digits = places,
      format = "d", flag = "0"
    )
  )
  
}

#' @noRd
# Ensure data has dimension names ----
# Updated 14.06.2023
ensure_dimension_names <- function(data)
{
  
  # Get dimensions
  dimensions <- dim(data)
  
  # Check for column names
  if(is.null(colnames(data))){
    
    # Standardize names
    colnames(data) <- paste0(
      "V", format_integer(
        seq_len(dimensions[2]),
        digits(dimensions[2]) - 1
      )
    )
    
  }
  
  # Check for matrix
  if(dimensions[1] == dimensions[2]){
    
    # Check for row names
    if(is.null(data) | any(row.names(data) != colnames(data))){
      
      # Assign column names to row names
      row.names(data) <- colnames(data)
      
    }
    
  }
  
  # Return named data
  return(data)
  
}

#' @noRd
# Transfer names from data to output ----
# Usually for matrix and data frame
# Updated 18.06.2023
transfer_names <- function(data, output)
{
  
  # Check for variable names
  if(!is.null(colnames(data))){
    
    # Add names to rows and columns
    colnames(output) <- row.names(output) <- colnames(data)
    
  }
  
  # Return named output
  return(output)
  
}

#' @noRd
# No names print ----
# Updated 03.02.2023
no_name_print <- function(object){
  
  # Convert object to data frame
  df <- as.data.frame(object)
  
  # Remove column names
  colnames(df) <- NULL
  
  # Print with no quotes or row names
  print(df, quote = FALSE, row.names = FALSE)
  
}

#' @noRd
# General function to check for packages ----
# Updated 14.06.2023
check_package <- function(packages)
{
  
  # # Check for packages
  # installed <- packages %in% row.names(installed.packages())
  
  # Performs what original `installed.packages()` does
  # but without additional fluff
  installed <- packages %in%
    unlist(
      sapply(.libPaths(), list.files, USE.NAMES = FALSE)
    )
  
  # Determine which packages are not installed
  not_installed <- packages[!installed]
  
  # Print error with missing packages
  if(length(not_installed) != 0){
    
    # Organize packages error output
    if(length(not_installed) > 1){
      
      # Get missing packages
      missing_packages <- paste0("{", packages , "}", collapse = ", ")
      packages <- paste0("\"", packages, "\"", collapse = ", ")
      
      # Stop and tell user to install package
      stop(
        paste0(
          missing_packages, 
          " are not installed but are required for this function. ",
          "Please run \n\n",
          "install.packages(c(", packages, "))",
          "\n\nOnce installed, re-run this function (you may need to restart R/RStudio)."
        )
      )
      
    }else{
      
      # Get missing packages
      missing_packages <- paste0("{", packages, "}")
      packages <- paste0("\"", packages, "\"")
      
      # Stop and tell user to install package
      stop(
        paste0(
          missing_packages, 
          " is not installed but is required for this function. ",
          "Please run \n\n",
          "install.packages(c(", packages, "))",
          "\n\nOnce installed, re-run this function (you may need to restart R/RStudio)."
        )
      )
      
    }
    
  }
  
}

#' @noRd
#'
# General function to silently obtain output ----
# Updated 11.05.2023
silent_call <- function(...){
  
  # Make call
  sink <- capture.output(
    result <- suppressWarnings(
      suppressMessages(
        ...
      )
    )
  )
  
  # Return result
  return(result)
  
}

#' @noRd
#'
# General function to silently load package ----
# Updated 10.06.2023
silent_load <- function(...){
  
  # Return result
  return(
    suppressPackageStartupMessages(...)
  )
  
}

#' @noRd
# Function to obtain arguments ----
# Updated 09.06.2023
obtain_arguments <- function(FUN, FUN.args)
{
  
  # Obtain formal arguments
  FUN.formals <- formals(FUN)
  
  # Check for input arguments
  if(length(FUN.args) != 0){
    
    ## Check for matching arguments
    if(any(names(FUN.args) %in% names(FUN.formals))){
      
      replace.args <- FUN.args[na.omit(match(names(FUN.formals), names(FUN.args)))]
      
      FUN.formals[names(replace.args)] <- replace.args
    }
    
  }
  
  # Remove ellipses
  if("..." %in% names(FUN.formals)){
    FUN.formals[which(names(FUN.formals) == "...")] <- NULL
  }
  
  # Remove call arguments (assume they are supplied elsewhere)
  call_argument <- sapply(FUN.formals, function(x){is(x, "call")})
  
  # Keep non-calls
  FUN.formals <- FUN.formals[!call_argument]
  
  # Return arguments
  return(FUN.formals)
  
}

# Function to check for usable variables ----
#' @noRd
# Updated 13.06.2023
usable_data <- function(data, verbose)
{
  
  # Turn data into a data matrix
  data_matrix <- data.matrix(data)
  
  # All missing after coercions
  remove_columns <- apply(data_matrix, 2, function(x){all(is.na(x))})
  
  # Send warning if there are any removals
  if(any(remove_columns)){
    
    # Send warning
    warning(
      paste0(
        "Some variables were could not to be coerced to numeric values. These variables have been removed from the analysis:\n",
        paste0("'", colnames(data_matrix)[remove_columns], "'", collapse = ", "),
        "\n\nIf these variables were not intended to be removed, then try converting them to numeric values before inputting the data into the function"
      )
    )
    
    # Remove these variables from `data` and `data_matrix`
    data <- data[,-remove_columns]
    data_matrix <- data_matrix[,-remove_columns]
    
  }
  
  # Determine whether each variable was coerced
  coercions <- !apply(data_matrix == data, 2, all, na.rm = TRUE)
  
  # Send warning for any coercions
  if(any(coercions) & isTRUE(verbose)){
    
    # Send warning
    warning(
      paste0(
        "Several variables were coerced to numeric values. These variables were changed to numeric values:\n",
        paste0("'", colnames(data_matrix)[remove_columns], "'", collapse = ", ")
      )
    )
    
    
  }
  
  # Send back usable data
  return(data_matrix)
  
}

#' @noRd
# Obtain data, sample size, correlation matrix ----
# Generic function to get the usual needed inputs
# Updated 23.06.2023
obtain_sample_correlations <- function(data, n, corr, na.data, verbose, ...)
{
  
  # Check if data is a correlation matrix
  if(!is_symmetric(data)){
    
    # Check for appropriate variables
    data <- usable_data(data, verbose)
    
    # Obtain sample size
    n <- dim(data)[1]
    
    # Check for automatic correlations
    if(corr == "auto"){
      
      # Compute correlations
      correlation_matrix <- auto.correlate(
        data = data, corr = "pearson",
        na.data = na.data, verbose = verbose,
        ...
      )
      
    }else{
      
      # Obtain correlations using base R
      correlation_matrix <- cor(data, use = na.data, method = corr)
      
    }
    
  }else{
    
    # Check for sample size
    if(is.null(n)){
      stop("A symmetric matrix was provided in the 'data' argument but the sample size argument, 'n', was not set. Please input the sample size into the 'n' argument.")
    }
    
    # If symmetric and sample size is provided, then
    # data to correlation matrix
    correlation_matrix <- data
    
  }
  
  # Set up results
  results <- list(
    data = data,
    n = n,
    correlation_matrix = correlation_matrix
  )
  
  # Return results
  return(results)
  
}

#' @noRd
# Standard EGA arguments ----
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
# Standard {lavaan} arguments ----
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
# Make unidimensional CFA model ----
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
# Determine estimator arguments ----
# Updated 03.23.2023
estimator_arguments <- function(lavaan_ARGS)
{
  
  # Obtain categories
  categories <- data_categories(
    data = lavaan_ARGS$data
  )
  
  # Check for categories
  if(any(categories <= 7)){
    
    # Set arguments
    lavaan_ARGS$estimator <- "WLSMV"
    lavaan_ARGS$missing <- "pairwise"
    lavaan_ARGS$ordered <- names(categories)[
      categories <= 7
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

#' @noRd
# Set default argument (cleaner missing function) ----
# Updated 15.06.2023
set_default <- function(argument, default, FUN)
{
  
  # Check for function first and foremost
  if(is.function(argument)){
    return(argument)
  }
  
  # Check for value error
  value_error(argument, base::mode(default))
  
  # Check for function
  if(!is.function(FUN)){
    choices <- FUN # assume choices were input to function 
  }else{
    
    # Get formal argument choices from default call
    choices <- formals(FUN)[[substitute(argument)]]
    
    # Check if choices is a call
    if(is.call(choices)){
      
      # Get choices
      choices <- tolower(as.character(choices))
      
      # Remove "c" from call
      choices <- choices[-which(choices == "c")]
      
    }
    
  }
  
  # Check for missing argument
  if(missing(argument)){
    argument <- default
  }else if(length(argument) != 1){
    argument <- default
  }
  
  # Force lowercase argument
  if(is.character(argument)){
    argument <- tolower(argument)
  }
  
  if(!argument %in% choices){
    stop(
      paste0(
        "Invalid argument: ", substitute(argument),
        " = \"", argument, "\"",
        "\n\nPlease choose from: ", 
        paste0(
          "\"", choices, "\"", collapse = ", "
        )
      )
    )
  }
  
  # Return argument
  return(argument)
  
}

#' @noRd
# Legacy Argument Handling ----
# Updated 15.06.2023
legacy_EGA_args <- function(ellipse)
{
  
  # Check for `model.args`
  if("model.args" %in% names(ellipse)){
    
    # Overwrite arguments
    ellipse <- overwrite_arguments(ellipse, ellipse$model.args)
    
    # Remove `model.args`
    ellipse <- ellipse[-which(names(ellipse) == "model.args")]
    
  }
  
  # Check for `algorithm.args`
  if("algorithm.args" %in% names(ellipse)){
    
    # Overwrite arguments
    ellipse <- overwrite_arguments(ellipse, ellipse$algorithm.args)
    
    # Remove `algorithm.args`
    ellipse <- ellipse[-which(names(ellipse) == "algorithm.args")]
    
  }
  
  # Return ellipse
  return(ellipse)
  
}

#' @noRd
# Overwrite Arguments (for Legacy) ----
# Updated 15.06.2023
overwrite_arguments <- function(main, ARGS)
{
  
  # Determine whether any arguments are in both
  if(any(names(ARGS) %in% names(main))){
    
    # Target arguments
    target_args <- intersect(names(ARGS), names(main))
    
    # Replace target arguments
    main[target_args] <- ARGS[target_args]
    
  }
  
  # Return main arguments
  return(main)
  
}

#%%%%%%%%%%%%%%%%%%%%
# PLOT FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Defaults for GGally plotting ----
# For plots and methods
# Updated 15.06.2023
GGally_args <- function(ellipse)
{
  
  # Get default `ggnet2` arguments
  default_args <- formals(
    silent_load(GGally::ggnet2)
    # Should always be first call to {GGally}
  )
  
  # Get default {EGAnet} arguments
  ega_default_args <- list(
    layout.exp = 0.20, label.size = 5,
    edge.label.color = "black", node.alpha = 0.50,
    node.shape = 19, node.size = 12, 
    edge.alpha = "edge.alpha", edge.size = 8
  )
  
  # Replace `ggnet2` arguments with {EGAnet} arguments
  default_args <- overwrite_arguments(default_args, ega_default_args)
  
  # Replace `ggnet2` arguments with arguments input
  default_args <- overwrite_arguments(default_args, ellipse)
  
  # Remove the ellipse
  if("..." %in% names(default_args)){
    default_args <- default_args[-which(names(default_args) == "...")]
  }
  
  # Various possible names for things
  ## Layout
  if("layout" %in% names(ellipse)){
    default_args$mode <- ellipse$layout  
  }
  
  ## Node transparency
  if("alpha" %in% names(ellipse)){
    default_args$node.alpha <- ellipse$alpha
  }
  
  ## Node color
  if("color" %in% names(ellipse)){
    default_args$node.color <- ellipse$color
  }
  
  ## Node shape
  if("shape" %in% names(ellipse)){
    default_args$node.shape <- ellipse$shape
  }
  
  ## Node size
  if("vsize" %in% names(ellipse)){
    default_args$node.size <- ellipse$vsize
  }
  
  ## Edge color
  if(!"edge.color" %in% names(ellipse)){
    default_args$edge.color <- c("darkgreen", "red")
  }
  
  ## Edge line types
  if(!"edge.lty" %in% names(ellipse)){
    default_args$edge.lty <- c("solid", "solid")
  }
  
  ## Color palette
  if(!"color.palette" %in% names(ellipse)){
    default_args$color.palette <- "polychrome"
  }else if(is.character(ellipse$color.palette)){
    
    # Check for gray scale options
    gray_options <- c(
      "greyscale", "grayscale", "colorblind"
    )
    
    # Check for gray scale
    if(tolower(ellipse$color.palette) %in% gray_options){
      default_args$edge.color <- c("#293132", "grey25")
      default_args$edge.lty <- c("solid", "dashed")
    }
      
  }
  # NOTE: gray scale will override node and edge colors
  # as well as line types
  
  # Return arguments
  return(default_args)
  
}

#' @noRd
# Error Checking for GGally plotting ----
# For plots and methods
# Updated 15.06.2023
GGally_errors <- function(
    plot_ARGS, dimensions,
    communities, non_zero_edges
)
{
  
  # Only the most common arguments are checked here
  # Edge case inputs are not considered
  
  # Determine number of nodes and edges
  nodes <- dimensions[2]
  edges <- length(non_zero_edges)
  
  ### Node arguments
  
  # Node Label Alpha
  value_error(plot_ARGS$label.alpha, "numeric")
  length_error(plot_ARGS$label.alpha, c(1, nodes))
  
  # Node Label Color
  value_error(plot_ARGS$label.color, "character")
  length_error(plot_ARGS$label.color, c(1, nodes))
  
  # Node Label Size
  value_error(plot_ARGS$label.size, "numeric")
  length_error(plot_ARGS$label.size, c(1, nodes))
  
  # Node Label
  value_error(plot_ARGS$node.label, "character")
  length_error(plot_ARGS$node.label, c(1, nodes))
  
  # Node Alpha
  value_error(plot_ARGS$node.alpha, "numeric")
  length_error(plot_ARGS$node.shape, c(1, communities, nodes))
  
  # Node Color
  value_error(plot_ARGS$node.color, "character")
  length_error(plot_ARGS$node.color, c(1, communities, nodes))
  
  # Node Shape
  value_error(plot_ARGS$node.shape, "numeric")
  length_error(plot_ARGS$node.shape, c(1, communities, nodes))
  
  # Node Size
  value_error(plot_ARGS$node.size, "numeric")
  length_error(plot_ARGS$node.size, c(1, communities, nodes))
  
  ### Edge arguments
  
  # Edge Alpha
  value_error(plot_ARGS$edge.alpha, "numeric")
  length_error(plot_ARGS$edge.alpha, c(1, edges))
  
  # Edge Color (allow two for positive and negative)
  value_error(plot_ARGS$edge.color, "character")
  length_error(plot_ARGS$edge.color, c(1, 2, edges))
  
  # Edge Size
  value_error(plot_ARGS$edge.size, "numeric")
  length_error(plot_ARGS$edge.size, c(1, edges))
  
  # Edge line type (allow two for positive and negative)
  value_error(plot_ARGS$edge.lty, "character")
  length_error(plot_ARGS$edge.lty, c(1, 2, edges))
  
}

#' @noRd
# Re-scale edges ----
# Updated 14.06.2023
rescale_edges <- function(network, edge_size)
{
  
  # Obtain absolute network (as a vector)
  vector_network <- abs(as.vector(network))
  
  # Set up scaling sequence 
  scale_sequence <- seq(0, 1, 0.0001)
  names(scale_sequence) <- scale_sequence
  
  # Set edge scaling (default `edge.size = 8`)
  edge_scaling <- scale_sequence * edge_size
  
  ## Remove names
  edge_scaling <- unname(edge_scaling[as.character(vector_network)])
  
  # Return scaled edges
  return(edge_scaling)
  
}

#' @noRd
# Readable names ----
# Updated 14.06.2023
readable_names <- function(node_names)
{
  
  # Split names
  name_split <- strsplit(node_names, split = " ")

  # Add return to names
  better_names <- sapply(name_split, function(x){
    
    # Obtain words in name
    words <- length(x)
    
    # Determine if split is necessary
    if(words > 1){
      
      # Determine number of lines
      add_line <- round(words / 2)
      
      # Paste back together name
      name <- paste(
        paste(x[seq_len(add_line)], collapse = " "),
        paste(x[(add_line+1):words], collapse = " "),
        sep = "\n"
      )
      
      # Return name
      return(name)
      
    }else{return(x)}
    
    
  })
  
  # Return names
  return(better_names)
  
}

#' @noRd
# Get network layout ----
# Updated 20.06.2023
get_layout <- function(
    network, dimensions,
    non_zero_index, plot_ARGS, ellipse
)
{
  
  # Determine whether "mode" or "layout" were used
  if(
    !"mode" %in% names(ellipse) &
    !"layout" %in% names(ellipse)
  ){ # Default: {qgraph} Fruchterman-Reingold
    
    # Lower triangle for edge list
    network_lower <- network[lower.tri(network)]
    weights_lower <- network_lower[network_lower != 0]
    
    # Set up edge list
    edge_list <- which(non_zero_index, arr.ind = TRUE)
    edge_list <- edge_list[edge_list[,"row"] < edge_list[,"col"],]
    edge_list <- edge_list[order(edge_list[,"row"]),]
    
    # Set layout (spring)
    network_layout <- qgraph::qgraph.layout.fruchtermanreingold(
      edgelist = edge_list,
      weights = abs(weights_lower / max(abs(weights_lower)))^2,
      vcount = dimensions[2]
    )
    
  }else{ 
    
    # Determine whether "mode" or "layout" is in arguments
    # If both, then override with "mode"
    if(!"mode" %in% names(ellipse) & "layout" %in% names(ellipse)){
      ellipse$mode <- ellipse$layout
    }
    
    # Check for whether mode was provided as character
    # or whether mode was input as 2D distance matrix
    if(is.character(ellipse$mode)){ # Obtain actual "mode" values using {sna}
      
      # Get layout function
      mode_FUN <- switch(
        tolower(plot_ARGS$mode),
        "adj" = sna::gplot.layout.adj,
        "circle" = sna::gplot.layout.circle,
        "circrand" = sna::gplot.layout.circrand,
        "eigen" = sna::gplot.layout.eigen,
        "fruchtermanreingold" = sna::gplot.layout.fruchtermanreingold,
        "geodist" = sna::gplot.layout.geodist,
        "hall" = sna::gplot.layout.hall,
        "kamadakawai" = sna::gplot.layout.kamadakawai,
        "mds" = sna::gplot.layout.mds,
        "princoord" = sna::gplot.layout.princoord,
        "random" = sna::gplot.layout.random,
        "rmds" = sna::gplot.layout.rmds,
        "segeo" = sna::gplot.layout.segeo,
        "seham" = sna::gplot.layout.seham,
        "spring" = sna::gplot.layout.spring,
        "springrepulse" = sna::gplot.layout.springrepulse,
        "target" = sna::gplot.layout.target
      )
      
      # Set network and arguments
      mode_ARGS <- list(
        d = network,
        layout.par = plot_ARGS$layout.par
      )
      
      # Obtain layout
      network_layout <- do.call(
        what = mode_FUN,
        args = mode_ARGS
      )
      
    }else if(is.numeric(ellipse$mode) & is.matrix(ellipse$mode)){
      
      # Assume "mode" is a 2D matrix corresponding to a layout
      network_layout <- ellipse$mode
      
    }
  
  }
  
  # Return layout
  return(network_layout)
  
}

#' @noRd
# Basic set up for plots ----
# Updated 23.06.2023
basic_plot_setup <- function(network, wc = NULL, ...)
{
  
  # Obtain ellipse arguments
  ellipse <- list(...)
  
  # Ensure network is a matrix
  network <- as.matrix(network)
  
  # Make sure network has a zero diagonal
  ## Mainly for `TMFG`
  diag(network) <- 0
  
  # Obtain network dimensions
  dimensions <- dim(network)
  
  # Set insignificant values to zero
  # (prevents `ggnet2` from erroring out)
  network <- round(network, 4)
  # Each digit of accuracy increases time 10x
  
  # Check for empty network
  if(sum(network) == 0){
    
    # Send message
    message("Network is empty. No plot produced.")
    
    # Return NULL
    return(NULL)
    
  }
  
  # Obtain number of communities
  communities <- length(na.omit(unique(wc)))
  
  # Obtain node names
  node_names <- colnames(network)

  # With packages, set up arguments
  plot_ARGS <- GGally_args(ellipse)
  
  # Set up the result of the plot arguments (runs in order of `ggnet2` arguments)
  ## Network
  plot_ARGS$net <- network
  
  # Set up networks for later use
  ## Full network
  non_zero_index <- network != 0
  non_zero_edges <- network[non_zero_index]
  
  ## Mode (layout)
  plot_ARGS$mode <- get_layout(
    network, dimensions, 
    non_zero_index,
    plot_ARGS, ellipse
  )
  
  ### Generic arguments (mostly handled in `GGally_args`)
    
  ## Remove some arguments
  plot_ARGS$alpha <- NULL; plot_ARGS$color <- NULL; plot_ARGS$size <- NULL
  
  ### Node arguments
  
  ## Color palette
  if(all(is.na(wc))){
    palette <- rep("grey", length(wc))
  }else{
    palette <- color_palette_EGA(plot_ARGS$color.palette, wc)
  }
  ## Set missing values to "white"
  palette[is.na(palette)] <- "white"
    
  ## Remove color palette
  color.palette <- plot_ARGS$color.palette
  plot_ARGS$color.palette <- NULL
  
  # Get number of node colors supplied
  node.color_length <- length(plot_ARGS$node.color)
  
  ## Set node color to communities
  if(all(plot_ARGS$node.color == "color")){
    
    # Use predefined palette
    plot_ARGS$node.color <- palette
    
  }else if(node.color_length == communities){
    
    # If number of node colors supplied is
    # for communities, then set them for each node
    plot_ARGS$node.color <- plot_ARGS$node.color[wc]
    
  }
  
  ## Set node label (default)
  if(all(plot_ARGS$node.label == "label")){
    plot_ARGS$node.label <- node_names
  }
  
  ## Set node size to zero (keep original node size)
  node.size <- plot_ARGS$node.size # handled in `GGally_args`
  plot_ARGS$node.size <- 0
  
  ### Edge arguments
  
  ## Set edge alpha (set to "edge.alpha" in `GGally_args`)
  if(all(plot_ARGS$edge.alpha == "edge.alpha")){
    plot_ARGS$edge.alpha <- sqrt(abs(non_zero_edges)) * 0.60
    # Not sure why `* 0.60` is needed to match old behavior
    # but without it the edges appear darker than original plots
  }
  
  ## Set edge color
  if(length(plot_ARGS$edge.color) == 2){
    plot_ARGS$edge.color <- ifelse(non_zero_edges >= 0, plot_ARGS$edge.color[1], plot_ARGS$edge.color[2])
  }
  
  ## Set edge line type
  if(length(plot_ARGS$edge.lty) == 2){
    plot_ARGS$edge.lty <- ifelse(non_zero_edges >= 0, plot_ARGS$edge.lty[1], plot_ARGS$edge.lty[2])
  }
  
  ## Set edge size (scale by `edge.size`)
  if(length(plot_ARGS$edge.size) == 1){
    plot_ARGS$edge.size <- rescale_edges(non_zero_edges, plot_ARGS$edge.size)
  }
  
  ## Edge label size (not used)
  plot_ARGS$edge.label.size <- ifelse(
    plot_ARGS$edge.label.size == "max_size/2",
    node.size / 2,
    plot_ARGS$edge.label.size
  )
  
  # Before call, check all arguments 
  # for any errors
  GGally_errors(
    plot_ARGS = plot_ARGS, dimensions = dimensions,
    communities = communities, non_zero_edges = non_zero_edges
  )
  
  # Get first layer with silent call
  first_layer <- silent_call(
    do.call(GGally::ggnet2, plot_ARGS)
  )
  
  # Return node size to `plot_ARGS` (was removed above)
  plot_ARGS$node.size <- node.size
  
  # Set up node names to be more readable
  node_names <- readable_names(plot_ARGS$node.label)

  # Determine border color
  ## Check for gray scale options
  gray_options <- c(
    "greyscale", "grayscale", "colorblind"
  )
  
  ## Set border color
  if(all(is.na(wc))){ # Plain network (without communities)
    border_color <- rep("black", dimensions[2])
  }else if( 
    length(color.palette) == 1 &
    color.palette %in% gray_options
  ){ # Gray scale network
    border_color <- ifelse(palette == "white", "white", "black")
  }else{ # Same color as nodes
    border_color <- plot_ARGS$node.color
  }
  
  # Custom nodes: transparent insides and dark borders
  second_layer <- first_layer +
    ggplot2::geom_point( # dark borders
      size = node.size, color = border_color,
      shape = 1, stroke = 1.5, alpha = 0.80
    ) +
    ggplot2::geom_point( # transparent insides
      size = node.size + 0.50, shape = 19,
      color = plot_ARGS$node.color,
      alpha = plot_ARGS$node.alpha
    ) +
    ggplot2::geom_text( # put text back on top
      ggplot2::aes(label = node_names), color = "black",
      size = plot_ARGS$label.size
    ) +
    ggplot2::guides( # create legend with these settings
      color = ggplot2::guide_legend(
        override.aes = list(
          color = unique(plot_ARGS$node.color),
          size = median(node.size, na.rm = TRUE),
          alpha = median(plot_ARGS$node.alpha, na.rm = TRUE),
          stroke = 1.5
        ),
        title = ifelse(
          "legend.title" %in% names(ellipse),
          ellipse$legend.title, ""
        )
      )
    )
  
  # Check for title
  if("title" %in% names(ellipse)){
    second_layer <- second_layer +
      ggplot2::labs(title = ellipse$title)
  }
  
  # Check for legend labels
  if("legend.names" %in% names(ellipse)){ # add user assigned names
    second_layer <- second_layer +
      ggplot2::scale_color_manual(
        values = unique(plot_ARGS$node.color),
        labels = ellipse$legend.names
      )
  }else if(all(is.na(wc))){ # no legend (network with no communities plot)
    second_layer <- second_layer +
      ggplot2::theme(legend.position = "none")
  }else{ # add membership names
    second_layer <- silent_call(
      second_layer +
      ggplot2::scale_color_manual(
        values = unique(plot_ARGS$node.color),
        labels = unique(wc)
      )
    )
  }
  
  # Set up return
  ## Hidden argument to return arguments plots
  ## Used most for comparing plots (same node placements)
  if("arguments" %in% names(ellipse) & isTRUE(ellipse$arguments)){
    
    # Set up return list
    results <- list(
      network_plot = second_layer,
      ARGS = plot_ARGS
    )
    
    # Return results
    return(results)
    
  }else{
    
    # Return plot only
    return(second_layer)
    
  }
  
}

#' @noRd
# Basic set up for single plot ----
# Updated 20.06.2023
single_plot <- function(network, wc, ...)
{
  
  # Look for memberships in arguments
  ## If no memberships, then plot network
  # as if all memberships are missing
  if(is.null(wc)){
    wc <- rep(NA, dimensions[2])
  }
  
  # Reorder network and communities
  new_order <- order(wc)
  network <- network[new_order, new_order]
  wc <- wc[new_order]
  
  # Send on and return from `basic_plot_setup`
  return(basic_plot_setup(network, wc, ...))
  
}

#' @noRd
# Dimension comparison for comparison plots ----
# Updated 20.06.2023
dimension_comparison <- function(original, comparison){
  
  # Get dimensions
  original_dimensions <- dim(original)
  comparison_dimensions <- dim(comparison)
  
  # Determine whether network to be compared has same
  # dimensions as the original plotted network
  if(any(original_dimensions != comparison_dimensions)){
    
    # Send error
    stop(
      paste0(
        "The original network's dimensions (",
        paste0(original_dimensions, collapse = " x "),
        ") do not match the comparison network's dimensions (",
        paste0(comparison_dimensions, collapse = " x "),
        ").\n\nDouble check to make sure the network dimensions match."
      )
    )
    
  }
  
  # Get names
  original_names <- colnames(original)
  comparison_names <- colnames(comparison)
  
  # Check for NULL
  if(!is.null(original_names) & is.null(comparison_names)){
    comparison_names <- original_names
  }else if(is.null(original_names) & !is.null(comparison_names)){
    original_names <- comparison_names
  }
  
  # Determine whether network to be compared has same
  # column names as the original plotted network
  not_matched <- !comparison_names %in% original_names 
  
  # Check for names that don't match
  if(any(not_matched)){
    
    # Obtain names that don't match in comparison
    no_match_names <- comparison_names[not_matched]
    
    # Send error
    stop(
      paste0(
        "Some variable names in the comparison network ",
        "did not match the original network: ",
        paste0("\"", no_match_names, "\"", collapse = ", ")
      )
    )
    
  }
  
  
}

#' @noRd
# Basic set up for comparing plots ----
# Updated 22.06.2023
compare_plots <- function(comparison_network, comparison_wc, plot_ARGS)
{
  
  # Comparison network
  comparison_network <- ega_first_level$network
  comparison_wc <- gsub("_.*", "", colnames(comparison_network))
  plot_ARGS <- first_level_plot$ARGS
  
  # Obtain the original network
  original_network <- plot_ARGS$net
  
  # Make sure dimensions are the same before proceeding
  dimension_comparison(original_network, comparison_network)
  
  # Ensure row names to ensure proper ordering
  row.names(comparison_network) <- colnames(comparison_network)
  
  # Put network into same order as original network
  matching_order <- match(
    colnames(original_network), # target to match
    colnames(comparison_network) # adjust to target
  )
  
  # Set comparison network in proper order
  ## Add to plot arguments
  plot_ARGS$network <- comparison_network[
    matching_order, matching_order
  ]
  
  # Also, set comparison memberships in proper order
  ## Add to plot arguments
  plot_ARGS$wc <- comparison_wc[matching_order]
  
  # Remove some arguments from `plot_ARGS`
  ## Essentially, the same call but allows some freedom
  plot_ARGS[c(
    "net", "node.color", "edge.alpha",
    "edge.color", "edge.size"
  )] <- NULL
  
  # Send on and return from `basic_plot_setup`
  return(do.call(basic_plot_setup, plot_ARGS))
  
}

#%%%%%%%%%%%%%%%%%%%%%
# ERROR FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Error for correlation matrix input ----
# Updated 02.02.2023
symmetric_matrix_error <- function(data, n){
  
  # Check for whether rows equal columns
  if(nrow(data) == ncol(data)){
    
    # Check for whether matrix is symmetric
    if(is_symmetric(data)){ # `is_symmetric` is in "helpers-general.R"
      
      # Check for whether "n" argument was provided
      if(is.null(n)){
        stop("A correlation matrix was detected as input into the 'data' argument but no 'n' was provided. Sample size must be provided for a correlation matrix.")
      }
      
    }
    
  }
  
}

#' @noRd
# Error for object type ----
# Updated 30.09.2022
object_error <- function(input, expected_type){
  
  # Check for possible object types
  possible_types <- sapply(
    X = expected_type,
    FUN = is,
    object = input
  )
  
  # Check for object types
  if(all(!possible_types)){
    stop(
      paste(
        "Input into '", deparse(substitute(input)),
        "' argument is not ", paste("'", expected_type, "'", sep = "", collapse = ", "),
        ". Input is ", paste("'", class(input), "'", sep = "", collapse = ", "),
        sep = ""
      )
    )
  }
  
}

#' @noRd
# Error for input value ----
# Updated 08.08.2022
value_error <- function(input, expected_value){
  
  # Check for value
  if(!is(input, expected_value)){
    stop(
      paste(
        "Input into '", deparse(substitute(input)),
        "' argument is not '", expected_type,
        "'. Input is ", paste("'", class(input), "'", sep = "", collapse = ", "),
        sep = ""
      )
    )
  }
  
}

#' @noRd
# Error for input length ----
# Updated 08.08.2022
length_error <- function(input, expected_lengths){
  
  # Check for length of input in expected length
  if(!length(input) %in% expected_lengths){
    stop(
      paste(
        "Length of '", deparse(substitute(input)),
        "' (", length(input),") does not match expected length(s). Length must be: ",
        paste("'", expected_lengths, "'", collapse = " or ", sep = ""),
        sep = ""
      )
    )
  }
  
}

#' @noRd
# Error for input range ----
# Updated 05.09.2022
range_error <- function(input, expected_ranges){
  
  # Obtain expected maximum and minimum values
  expected_maximum <- max(expected_ranges)
  expected_minimum <- min(expected_ranges)
  
  # Obtain maximum and minimum values
  actual_maximum <- round(max(input), 3)
  actual_minimum <- round(min(input), 3)
  
  # Check for maximum of input in expected range
  if(actual_maximum > expected_maximum){
    stop(
      paste(
        "Maximum of '", deparse(substitute(input)),
        "' (", actual_maximum,") does not match expected range(s). Range must be between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        sep = ""
      )
    )
  }
  
  # Check for maximum of input in expected range
  if(actual_minimum < expected_minimum){
    stop(
      paste(
        "Minimum of '", deparse(substitute(input)),
        "' (", actual_minimum,") does not match expected range(s). Range must be between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        sep = ""
      )
    )
  }
  
}

#%%%%%%%%%%%%%%%%%%%%
# MATH FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Trace of matrix ----
# Updated 27.05.2023
trace <- function(object)
{
  
  # Ensure matrix
  if(!is(object, "matrix") | !is(object, "table")){
    object <- as.matrix(object)
  }
  
  # Return trace
  return(sum(diag(object)))
  
}

#%%%%%%%%%%%%%%%%%%%%%%%
# SYSTEM FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Custom Parallelization ----
# Updated 10.05.2023
parallel_process <- function(
    datalist, # list of data
    iter = NULL, # number of iterations
    progress = TRUE, # progress bar
    FUN, # function to use
    FUN_args = list(), # arguments to use in function
    export = NULL, # variables to export (if necessary)
    ncores # number of cores
){
  
  # Obtain arguments
  FUN_args <- obtain.arguments(
    FUN = FUN,
    FUN.args = FUN_args
  )
  
  # Set progress bar up
  if(isTRUE(progress)){
    
    # Calculate total computations
    total_computations <- ifelse(
      is.null(iter), length(datalist), iter
    )
    
    # Count computations
    count_computations <- 0
    
    # Create data splits (necessary for progress bar)
    if(total_computations <= 100){
      
      # Split computations
      split_computations <- ncores
      
      # Set start and end points of data splits
      split_start <- seq(1, total_computations, split_computations)
      split_end <- unique(
        c(
          seq(split_computations, total_computations, split_computations),
          total_computations
        )
      )
      
      # Initialize split list
      data_split <- vector("list", length = length(split_start))
      
      # Populate split list
      for(i in seq_along(data_split)){
        data_split[[i]] <- datalist[split_start[i]:split_end[i]]
      }
      
      # Initialize results list
      results <- vector("list", length = length(data_split))
      
      # Initialize runtime updates
      runtime_update <- seq(
        0, total_computations, ncores
      )
      
      # Obtain runtime updates
      runtime_update <- unique(c(runtime_update, total_computations))
      
    }else{
      
      # Split computations
      split_computations <- ncores
      
      # Set start and end points of data splits
      split_start <- seq(1, total_computations, split_computations)
      split_end <- unique(
        c(
          seq(split_computations, total_computations, split_computations),
          total_computations
        )
      )
      
      # Batch splits
      batch_computations <- round(total_computations / 100)
      
      # Set start and end points of data batches
      batch_start <- seq(1, length(split_end), batch_computations)
      batch_end <- unique(
        c(
          seq(batch_computations, length(split_end), batch_computations),
          length(split_start)
        )
      )
      
      # Initialize batch list
      data_split <- vector("list", length = length(batch_start))
      
      # Populate split list
      for(j in seq_along(data_split)){
        data_split[[j]] <- datalist[
          split_start[batch_start[j]]:split_end[batch_end[j]]
        ]
      }
      
      # Initialize results list
      results <- vector("list", length = length(data_split))
      
      # Initialize runtime updates
      runtime_update <- seq(
        0, total_computations, length(data_split[[1]])
      )
      
      # Obtain runtime updates
      runtime_update <- unique(c(runtime_update, total_computations))
      
    }
    
    # Obtain start time
    if(count_computations == 0){
      start_time <- Sys.time()
    }
    
    # Plan parallelization
    future::plan(
      strategy = "multisession",
      workers = ncores
    )
    
    # Loop through data splits
    for(i in seq_along(data_split)){
      
      # Update progress
      if(count_computations < runtime_update[2]){
        
        # Update progress
        custom_progress(
          i = count_computations,
          max = total_computations,
          start_time = "calculating"
        )
        
      }
      
      # Run parallelization
      results[[i]] <- future.apply::future_lapply(
        X = data_split[[i]],
        FUN = function(x, FUNC, FUN_args){
          
          # Add data
          FUN_args[[names(FUN_args)[1]]] <- x
          
          # Return results from function
          return(do.call(FUNC, as.list(FUN_args)))
          
        },
        FUNC = FUN, FUN_args = FUN_args,
        future.stdout = FALSE,
        future.packages = "EGAnet",
        future.seed = NULL
      )
      
      # Update computation count
      count_computations <- count_computations +
        length(data_split[[i]])
      
      # Update progress
      if(count_computations %in% runtime_update){
        
        custom_progress(
          i = count_computations,
          max = total_computations,
          start_time = start_time
        )
        
      }
      
    }
    
    # Unwrap parallelization
    results <- unlist(results, recursive = FALSE)
    
  }else{ # Run without progress
    
    # Plan parallelization
    future::plan(
      strategy = "multisession",
      workers = ncores
    )
    
    # Run parallelization
    results <- future.apply::future_lapply(
      X = datalist,
      FUN = function(x, FUNC, FUN_args){
        
        # Add data
        FUN_args[[names(FUN_args)[1]]] <- x
        
        # Return results from function
        return(do.call(FUNC, as.list(FUN_args)))
        
      },
      FUNC = FUN, FUN_args = FUN_args,
      future.stdout = FALSE,
      future.packages = "EGAnet",
      future.seed = NULL
    )
    
  }
  
  # Return results
  return(results)
  
}

#' @importFrom utils packageVersion
#' @noRd
# Error Report ----
# Updated 10.05.2023
error.report <- function(result, SUB_FUN, FUN)
{
  # Let user know that an error has occurred
  message(paste("\nAn error has occurred in the '", SUB_FUN, "' function of '", FUN, "':\n", sep =""))
  
  # Give them the error to send to you
  cat(paste(result))
  
  # Tell them where to send it
  message("\nPlease open a new issue on GitHub (bug report): https://github.com/hfgolino/EGAnet/issues/new/choose")
  
  # Give them information to fill out the issue
  OS <- as.character(Sys.info()["sysname"])
  OSversion <- paste(as.character(Sys.info()[c("release", "version")]), collapse = " ")
  Rversion <- paste(R.version$major, R.version$minor, sep = ".")
  EGAversion <- paste(unlist(packageVersion("EGAnet")), collapse = ".")
  
  # Let them know to provide this information
  message(paste("\nBe sure to provide the following information:\n"))
  
  # To reproduce
  message(styletext("To Reproduce:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " Function error occurred in: ", SUB_FUN, " function of ", FUN, sep = ""))
  
  # R and {EGAnet} version
  message(styletext("\nR and {EGAnet} versions:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " R version: ", Rversion, sep = ""))
  message(paste(" ", textsymbol("bullet"), " {EGAnet} version: ", EGAversion, sep = ""))
  
  # Desktop
  message(styletext("\nOperating System:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " OS: ", OS, sep = ""))
  message(paste(" ", textsymbol("bullet"), " Version: ", OSversion, sep = ""))
  
}

#' @noRd
# OS and System Check ----
# Updated 08.09.2020
system.check <- function (...)
{
  OS <- unname(tolower(Sys.info()["sysname"]))
  
  RSTUDIO <- ifelse(Sys.getenv("RSTUDIO") == "1", TRUE, FALSE)
  
  TEXT <- TRUE
  
  if(!RSTUDIO){if(OS != "linux"){TEXT <- FALSE}}
  
  res <- list()
  
  res$OS <- OS
  res$RSTUDIO <- RSTUDIO
  res$TEXT <- TEXT
  
  return(res)
}

#' @noRd
# Colorize text ----
# Updated 08.09.2020
colortext <- function(text, number = NULL, defaults = NULL)
{
  # Check system
  sys.check <- system.check()
  
  if(sys.check$TEXT)
  {
    # Defaults for number (white text)
    if(is.null(number) || number < 0 || number > 231)
    {number <- 15}
    
    # Check for default color
    if(!is.null(defaults))
    {
      # Adjust highlight color based on background color
      if(defaults == "highlight")
      {
        if(sys.check$RSTUDIO)
        {
          
          if(rstudioapi::getThemeInfo()$dark)
          {number <- 226
          }else{number <- 208}
          
        }else{number <- 208}
      }else{
        
        number <- switch(defaults,
                         message = 204,
                         red = 9,
                         orange = 208,
                         yellow = 11,
                         "light green" = 10,
                         green = 34,
                         cyan = 14,
                         blue = 12,
                         magenta = 13,
                         pink = 211,
        )
        
      }
      
    }
    
    return(paste("\033[38;5;", number, "m", text, "\033[0m", sep = ""))
    
  }else{return(text)}
}

#' @noRd
# Style text ----
# Updated 08.09.2020
styletext <- function(text, defaults = c("bold", "italics", "highlight",
                                         "underline", "strikethrough"))
{
  # Check system
  sys.check <- system.check()
  
  if(sys.check$TEXT)
  {
    if(missing(defaults))
    {number <- 0
    }else{
      
      # Get number code
      number <- switch(defaults,
                       bold = 1,
                       italics = 3,
                       underline = 4,
                       highlight = 7,
                       strikethrough = 9
      )
      
    }
    
    return(paste("\033[", number, ";m", text, "\033[0m", sep = ""))
  }else{return(text)}
}

#' @noRd
# Symbols ----
# Updated 24.04.2020
textsymbol <- function(symbol = c("alpha", "beta", "chi", "delta",
                                  "eta", "gamma", "lambda", "omega",
                                  "phi", "pi", "rho", "sigma", "tau",
                                  "theta", "square root", "infinity",
                                  "check mark", "x", "bullet")
)
{
  # Get number code
  sym <- switch(symbol,
                alpha = "\u03B1",
                beta = "\u03B2",
                chi = "\u03C7",
                delta = "\u03B4",
                eta = "\u03B7",
                gamma = "\u03B3",
                lambda = "\u03BB,",
                omega = "\u03C9",
                phi = "\u03C6",
                pi = "\u03C0",
                rho = "\u03C1",
                sigma = "\u03C3",
                tau = "\u03C4",
                theta = "\u03B8",
                "square root" = "\u221A",
                infinity = "\u221E",
                "check mark" = "\u2713",
                x = "\u2717",
                bullet = "\u2022"
  )
  
  return(sym)
}

#' @noRd
# Title Case ----
# Updated 16.06.2023
totitle <- function(string)
{
  
  # Split by spaces
  words <- unlist(strsplit(string, split = " "))
  
  # Set first letters to uppercase
  titleCased <- sapply(words, function(x){
    
    # Stitch together letters
    return(
      paste0(
        toupper(substr(x, 1, 1)),
        tolower(substr(x, 2, nchar(x)))
      )
    )

  })
  
  # Paste words back together
  return(paste(titleCased, collapse = " "))
  
}
