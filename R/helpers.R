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
# Ensure data has dimension names ----
# Updated 03.02.2023
ensure_dimension_names <- function(data)
{
  
  # Check for column names
  if(is.null(colnames(data))){
    
    # Standardize names
    colnames(data) <- paste0(
      "V", formatC(
        x = 1:ncol(data),
        digits = (digits(ncol(data)) - 1),
        format = "d", flag = "0"
      )
    )
    
  }
  
  # Check for matrix
  if(nrow(data) == ncol(data)){
    
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
# Updated 12.04.2023
check_package <- function(packages)
{
  
  # Check for packages
  installed <- packages %in% row.names(installed.packages())
  
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
# Generate multivariate normal data ----
# Removes MASS package dependency (from version 7.3.54)
# Updated 23.12.2021
MASS_mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE)
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  if (EISPACK)
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    drop(X)
  else t(X)
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