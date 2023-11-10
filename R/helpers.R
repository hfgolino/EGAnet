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

#%%%%%%%%%%%%%%%%%%%%%%
# *APPLY FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

# These functions are merely wrappers over the `*apply`
# family to execute specific tasks that are often repeatedly
# done. Their purpose are to provide slightly cleaner code.

# These `*vapply` are specific versions of `vapply` that
# pre-specify the output. These versions are generally
# simpler (in function) and marginally faster than 
# the more commonly used `sapply`. They are more strict
# in that they require knowing the output format (e.g., numeric)

#' @noRd
# Character vector/matrix output apply ----
# Updated 06.07.2023
cvapply <- function(X, FUN, ..., LENGTH = 1, USE.NAMES = TRUE)
{
  return(vapply(X = X, FUN = FUN, FUN.VALUE = character(LENGTH), ..., USE.NAMES = USE.NAMES))
}

#' @noRd
# Logical vector/matrix output apply ----
# Updated 06.07.2023
lvapply <- function(X, FUN, ..., LENGTH = 1, USE.NAMES = TRUE)
{
  return(vapply(X = X, FUN = FUN, FUN.VALUE = logical(LENGTH), ..., USE.NAMES = USE.NAMES))
}

#' @noRd
# Numeric vector/matrix output apply ----
# Updated 06.07.2023
nvapply <- function(X, FUN, ..., LENGTH = 1, USE.NAMES = TRUE)
{
  return(vapply(X = X, FUN = FUN, FUN.VALUE = numeric(LENGTH), ..., USE.NAMES = USE.NAMES))
}

#' @noRd
# Unlist `lapply` ----
# Updated 14.07.2023
ulapply <- function(X, FUN, ..., recursive = TRUE)
{
  return(unlist(lapply(X, FUN, ...), recursive = recursive))
}

# The `row_apply` and `column_apply` functions are
# not necessary faster than `apply` but they are
# more strict in that they return a matrix with
# the same dimensions that was input into 'X'.
# Therefore, there is never a need to transpose
# the output because it will always have the same
# dimensions as what went in

#' @noRd
# Stricter `apply` by row ----
# Slightly faster than `apply` but not appreciably
# Updated 06.07.2023
row_apply <- function(X, FUN, ...)
{
  
  # Ensure matrix
  X <- as.matrix(X)
  
  # Get dimensions of data
  dimensions <- dim(X)
  
  # Get appropriate function
  vapply_FUN <- switch( # in order of most likely
    typeof(X),
    "double" = nvapply,
    "integer" = nvapply,
    "logical" = lvapply,
    "character" = cvapply,
    stop(
      paste0(
        "\"", typeof(X), "\"",
        " is not available for `row_apply`. Only \"character\", \"logical\", and \"numeric\" are available" 
      ), call. = FALSE
    )
  )
  
  # Return call
  return(
    matrix(
      vapply_FUN(
        split(X, seq_len(dimensions[1])),
        FUN = FUN, ..., LENGTH = dimensions[2],
        USE.NAMES = FALSE # names are supplied below
      ), 
      dimensions, dimnames = dimnames(X), byrow = TRUE
    )
  )
  
}

#' @noRd
# Stricter `apply` by column ----
# Slightly faster than `apply` but not appreciably
# Updated 06.07.2023
column_apply <- function(X, FUN, ...)
{
  
  # Ensure matrix
  X <- as.matrix(X)
  
  # Get dimensions of data
  dimensions <- dim(X)
  
  # Get appropriate function
  vapply_FUN <- switch( # in order of most likely
    typeof(X),
    "double" = nvapply,
    "integer" = nvapply,
    "logical" = lvapply,
    "character" = cvapply,
    stop(
      paste0(
        "\"", typeof(X), "\"",
        " is not available for `row_apply`. Only \"character\", \"logical\", and \"numeric\" are available" 
      ), call. = FALSE
    )
  )
  
  # Return call
  return(
    matrix(
      vapply_FUN(
        split(t(X), seq_len(dimensions[2])),
        FUN = FUN, ..., LENGTH = dimensions[1],
        USE.NAMES = FALSE # names are supplied below
      ), 
      dimensions, dimnames = dimnames(X), byrow = FALSE
    )
  )
  
}

#' @noRd
# 3D Array apply ----
# Apply's a function across a 3D array to obtain a 2D array
# About 2x faster than `apply(simplify2array(X), 1:2, FUN)`
# Updated 06.07.2023
symmetric_matrix_lapply <- function(X, FUN, ...){
  
  # Get appropriate function
  vapply_FUN <- switch( # in order of most likely
    typeof(X[[1]]),
    "double" = nvapply,
    "integer" = nvapply,
    "logical" = lvapply,
    "character" = cvapply,
    stop(
      paste0(
        "\"", typeof(X[[1]]), "\"",
        " is not available for `symmetric_matrix_lapply`. Only \"character\", \"logical\", and \"numeric\" are available" 
      ), call. = FALSE
    )
  )
  
  # Get dimensions of single matrix in list
  dimensions <- dim(X[[1]])
  
  # Pre-obtain lower triangle indices
  # Equivalent to: (`lower.tri(matrix, diag = TRUE)`)
  lower_triangle_index <- .row(dimensions) >= .col(dimensions)
  
  # Get lower triangles
  lower_matrix <- do.call(rbind, lapply(X, function(x){x[lower_triangle_index]}))
  
  # Apply function
  values <- vapply_FUN(as.data.frame(lower_matrix), FUN, ...)
  
  # Initialize return matrix
  new_matrix <- matrix(
    nrow = dimensions[1], ncol = dimensions[2],
    dimnames = dimnames(X[[1]])
  )
  
  # Add values to lower triangle
  new_matrix[lower_triangle_index] <- values
  
  # Transpose 
  new_matrix <- t(new_matrix)
  
  # Add values again to lower triagle
  new_matrix[lower_triangle_index] <- values
  
  # Return new matrix
  return(new_matrix)
  
}

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
# Assistance in writing C and interfacing with R
# was provided by GPT-4

#' @noRd
# Generate integer seeds ----
# (non-randomly based on `seed`)
# `NULL` is reserved for actual pseudo-randomness
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 27.07.2023
reproducible_seeds <- function(n, seed = NULL)
{
  
  # Return call from C
  return(
    .Call(
      "r_xoshiro_seeds",
      as.integer(n), swiftelse(is.null(seed), 0, seed),
      PACKAGE = "EGAnet"
    )
  )
  
}

#' @noRd
# Generate uniform data ----
# Allows adjustment of range
# Updated 26.10.2023
runif_xoshiro <- function(n, min = 0, max = 1, seed = NULL)
{
  
  # Get values
  values <- .Call(
    "r_xoshiro_uniform",
    as.integer(n),
    swiftelse(is.null(seed), 0, seed),
    PACKAGE = "EGAnet"
  )
  
  # Check for changes to minimum and maximum
  if(min != 0 || max != 1){ # transform
    values <- min + (max - min) * values
  }
  
  
  # Return call from C
  return(values)
  
}

#' @noRd
# Shuffle (without replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 30.07.2023
shuffle <- function(x, size = length(x), seed = NULL)
{
  
  # Return call from C
  return(
    x[.Call(
        "r_xoshiro_shuffle",
        as.integer(seq_along(x)), 
        swiftelse(is.null(seed), 0, seed),
        PACKAGE = "EGAnet"
    )][seq_len(size)]
  )
  
}

#' @noRd
# Shuffle (with replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 30.07.2023
shuffle_replace <- function(x, size = length(x), seed = NULL)
{
  
  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle_replace",
      x, swiftelse(is.null(seed), 0, seed),
      PACKAGE = "EGAnet"
    )][seq_len(size)]
  )
  
}

#' @noRd
# Random normal generation with Ziggurat ----
# https://people.sc.fsu.edu/~jburkardt/cpp_src/ziggurat/ziggurat.html
# Updated 27.07.2023
rnorm_ziggurat <- function(n, seed = NULL)
{
  
  # Return call from C
  return(
    .Call(
      "r_ziggurat",
      as.integer(n), 
      swiftelse(is.null(seed), 0, seed),
      PACKAGE = "EGAnet"
    )
  )
  
}

#' @noRd
# Generate multivariate normal data ----
# Removes MASS package dependency (from version 7.3.60)
# Original code here for legacy and bug checking
# Updated 03.07.2023
MASS_mvrnorm <- function(
    n = 1, mu, Sigma,
    tol = 1e-06, empirical = FALSE, EISPACK = FALSE
)
{
  
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  if (EISPACK) 
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1]))) 
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
    nm <- dn[[1]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
  
}

#' @noRd
# Generate multivariate normal data (quicker) ----
# A very aggressive optimization
# Removes 'tol', 'empirical' and 'EISPACK' arguments
# Removes safeguards like positive definite (should be from `auto.correlate`)
# 'empirical' is ALWAYS = FALSE
# Pre-computes `np` (n * p)
# Avoids repeated calls to `eigen` in `MASS_mvrnorm`
# Pre-computes `eS$vectors %*% diag(sqrt(pmax(ev, 0)), p)`
# Updated 27.07.2023
mvrnorm_precompute <- function(cases, Sigma)
{
 
  # Get p (number of variables)
  p <- dim(Sigma)[2]
  
  # Compute eigenvalues
  eS <- eigen(Sigma, symmetric = TRUE)
  
  # Get eigenvalues
  eigenvalues <- eS$values
  
  # Get non-zero eigenvalues
  non_zero <- eigenvalues > 0
  
  # Set negative eigenvalues to zero
  eigenvalues[!non_zero] <- 0
  
  # Get square root of non-zero eigenvalues
  eigenvalues[non_zero] <- sqrt(eigenvalues[non_zero])

  # Return pre-computed values
  return(
    list(
      p = p,
      np = cases * p,
      coV = eS$vectors %*% diag(eigenvalues, p)
    )
  )
  
}

#' @noRd
# Generate multivariate normal data (quick) ----
# Updated 27.07.2023
MASS_mvrnorm_quick <- function(seed = NULL, p, np, coV)
{
  return(t(tcrossprod(coV, matrix(rnorm_ziggurat(np, seed), ncol = p))))
}

#' @noRd
# Generate reproducible bootstrap data ----
# Wrapper for `reproducible_parametric` and `reproducible_resampling`
# Updated 30.07.2023
reproducible_bootstrap <- function(
    seed = NULL, data = NULL, case_sequence = NULL,
    mvrnorm_parameters = NULL,
    type = c("parametric", "resampling")
)
{
  
  # Based on bootstrap type, generate data
  if(type == "parametric"){
    
    # Return parametric samples
    return(
      do.call(
        what = MASS_mvrnorm_quick,
        args = c(
          list(seed = seed),
          mvrnorm_parameters
        )
      )
    )
    
  }else if(type == "resampling"){
    
    # Return resampling samples
    return(data[shuffle_replace(x = case_sequence, seed = seed),])
    
  }
  
}

#%%%%%%%%%%%%%%%%%%%%%
# PARALLELIZATION ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Clear memory ----
# Updated 07.11.2023
clear_memory <- function()
{
  sink <- capture.output(
    gc(verbose = FALSE, reset = TRUE, full = TRUE)
  )
}

#' @noRd
# Get available memory ----
# Updated 18.10.2023
available_memory <- function()
{

  # Get operating system
  OS <- tolower(Sys.info()["sysname"])

  # Branch based on OS
  if(OS == "windows"){ # Windows
    
    # Alternative (outputs memory in kB)
    bytes <- as.numeric(
      trimws(system("wmic OS get FreePhysicalMemory", intern = TRUE))[2]
    ) * 1e+03

  }else if(OS == "linux"){ # Linux
    
    # Split system information
    info_split <- strsplit(system("free", intern = TRUE), split = " ")
    
    # Remove "Mem:" and "Swap:"
    info_split <- lapply(info_split, function(x){gsub("Mem:", "", x)})
    info_split <- lapply(info_split, function(x){gsub("Swap:", "", x)})
    
    # Get actual values
    info_split <- lapply(info_split, function(x){x[x != ""]})
    
    # Bind values
    info_split <- do.call(rbind, info_split[1:2])
    
    # Get free values (Linux reports in *kilo*bytes -- thanks, Aleksandar Tomasevic)
    bytes <- as.numeric(info_split[2, info_split[1,] == "available"]) * 1e+03
    
  }else{ # Mac
    
    # System information
    system_info <- system("top -l 1 -s 0 | grep PhysMem", intern = TRUE)

    # Get everything after comma
    unused <- gsub(" .*,", "", system_info)
    
    # Get values only
    value <- gsub(" unused.", "", gsub("PhysMem: ", "", unused))
    
    # Check for bytes
    if(grepl("M", value)){
      bytes <- as.numeric(gsub("M", "", value)) * 1e+06
    }else if(grepl("G", value)){
      bytes <- as.numeric(gsub("G", "", value)) * 1e+09
    }else if(grepl("K", value)){
      bytes <- as.numeric(gsub("K", "", value)) * 1e+03
    }else if(grepl("B", value)){ # edge case
      bytes <- as.numeric(gsub("B", "", value)) * 1
    }else if(grepl("T", value)){ # edge case
      bytes <- as.numeric(gsub("T", "", value)) * 1e+12
    }
      
  }
  
  # Return bytes
  return(bytes)

}

#' @noRd
# Store random state ----
# Updated 04.08.2023
store_state <- function()
{
  if(exists(".Random.seed", envir = parent.frame(2))){
    assign("saved_state", .Random.seed, envir = parent.frame(2))
  }
}

#' @noRd
# Restore and remove random state ----
# Updated 04.08.2023
restore_state <- function()
{
  if(exists("saved_state", envir = parent.frame(2))){
    saved_state <- get("saved_state", envir = parent.frame(2))
    assign(".Random.seed", saved_state, envir = parent.frame(2))
    rm("saved_state", envir = parent.frame(2))
  }
}

#' @noRd
# Wrapper for parallelization ----
# Updated 24.10.2023
parallel_process <- function(
    iterations, # number of iterations
    datalist = NULL, # list of data
    FUN, # function to use
    ..., # ellipse arguments to pass on
    export = TRUE, # variables to export (if necessary)
    packages = "EGAnet", # always uses {EGAnet}
    ncores, # number of cores
    progress = TRUE, # progress bar
    clear = FALSE # whether progress bar should be cleared
){
  
  # Get available memory
  memory_available <- try(
    available_memory(),
    silent = TRUE
  )

  # In case the memory check fails
  if(!is(memory_available, "try-error")){
    
    # Check for global environment size
    if(export){ # needs `isTRUE` in case of character vector
      global_memory_usage <- sum(nvapply(ls(),function(x){object.size(get(x))}))
    }else{
      global_memory_usage <- sum(nvapply(ls()[ls() %in% export],function(x){object.size(get(x))}))
    }
    
    # Check for memory overload
    if(memory_available < global_memory_usage * ncores){
      stop(
        paste0(
          "Available memory (", byte_digits(memory_available), ") is less than ",
          "the amount of memory needed to perform parallelization: ",
          byte_digits(global_memory_usage * ncores), ".\n\n",
          "Lower the number of cores (`ncores`) or perform ",
          "batches of your operation."
        ), call. = FALSE
      )
    }
    
    # Set max size
    options(future.globals.maxSize = memory_available)
    
  }
  
  # Set up plan
  future::plan(
    strategy = "multisession",
    workers = ncores
  )
  
  # Check for progress
  if(isTRUE(progress)){
    
    # Set up handler
    progressr::handlers(
      progressr::handler_progress(
        format = ":spin [:bar] :percent elapsed: :elapsed ~remaining: :eta",
        clear = clear
      )
    )
    
    # Run progress locally
    progressr::with_progress({
      
      # Initialize progress bar
      progressbar <- progressr::progressor(iterations / ncores)
      
      # Perform parallel processing
      results <- future.apply::future_lapply(
        X = seq_len(iterations),
        function(iteration){
          
          # Update progress with full cores completion
          # (rather than every completion)
          if(iteration %% ncores == 0){
            progressbar()
          }
          
          # Return results
          if(!is.null(datalist)){
            return(
              silent_call( # Ensures quiet run with all arguments passed on
                FUN(datalist[[iteration]], ...)
              )
            )
          }else{
            return(
              silent_call( # Ensures quiet run with all arguments passed on
                FUN(...)
              )
            )
          }
          
        },
        future.globals = export,
        future.packages = packages,
        future.seed = NULL
      )
      
    }, enable = TRUE)
    
  }else{
    
    # Perform parallel processing
    results <- future.apply::future_lapply(
      X = seq_len(iterations),
      function(iteration){
        
        # Return results
        if(!is.null(datalist)){
          return(
            silent_call( # Ensures quiet run with all arguments passed on
              FUN(datalist[[iteration]], ...)
            )
          )
        }else{
          return(
            silent_call( # Ensures quiet run with all arguments passed on
              FUN(...)
            )
          )
        }
        
      },
      future.globals = export,
      future.packages = packages,
      future.seed = NULL
    )
    
  }
  
  # Force close of connections
  future::plan(strategy = "sequential")
  
  # Return results
  return(results)
  
}

#%%%%%%%%%%%%%%%%%%%%%%
# FASTER SEQUENCES ----
#%%%%%%%%%%%%%%%%%%%%%%

# These functions aren't often used because
# dimensions of the data are usually used
# elsewhere in {EGAnet}; however, these 
# functions use `seq_len` and `dim` which
# have minimal speed advantages over their
# respective `1:nrow` and `1:ncol`
#
# The primary benefit of using this functions
# is `seq_len` which throws an error with
# `seq_len(0)` where `1:0` will not
# (at least not at the top of the call)

#' @noRd
# Faster row sequences ----
# Updated 06.07.2023
nrow_sequence <- function(data)
{
  return(seq_len(dim(data)[1]))
}


#' @noRd
# Faster column sequences ----
# Updated 06.07.2023
ncol_sequence <- function(data)
{
  return(seq_len(dim(data)[2]))
}

#' @noRd
# Faster `ifelse` ----
# For single value replacements:
# 1.5x faster with 1 value
# 2.5x faster with 10 values
# >= 18x faster with >= 100 values
# Updated 24.07.2023
swiftelse <- function(condition, true, false)
{
  
  # Get condition length
  condition_length <- length(condition)
  
  # Check for single value
  if(condition_length == 1){
  
    # If TRUE
    if(condition){
      return(true)
    }else{ # Otherwise, FALSE
      return(false)
    }
    
  }
  
  # Initialize result
  result <- vector(mode(true), condition_length)
  
  # Set TRUE condition
  if(length(true) == 1){
    result[condition] <- true
  }else{
    result[condition] <- true[condition]
  }
  
  # Set FALSE condition
  if(length(false) == 1){
    result[!condition] <- false
  }else{
    
    # Get opposite condition (slightly faster than repeated calls)
    opposite_condition <- !condition
    result[opposite_condition] <- false[opposite_condition]
  }
  
  
  # Return result
  return(result)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# FORMATTING FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Determine decimal places before non-zero ----
# From StackOverflow: https://stackoverflow.com/questions/35553244/count-leading-zeros-between-the-decimal-point-and-first-nonzero-digit
# Updated 10.11.2023
leading_zero <- function(number)
{
  floor(-log10(.Machine$double.eps + abs(number) - floor(abs(number))))
}

#' @noRd
# Determine number of digits in a number ----
# Updated 24.07.2023
digits <- function(number)
{
  # Obtain the lowest value of log base 10 and add 1
  return(floor(log10(abs(number))) + 1)
  
}

#' @noRd
# Determine format of bytes ----
# Updated 22.07.2023
byte_digits <- function(bytes)
{
  
  # Get number of digits
  ndigits <- digits(bytes)
  
  # Loop over values
  if(ndigits > 11){
    value <- paste0(round(bytes / 1e+12, 2), " TB")
  }else if(ndigits > 8){
    value <- paste0(round(bytes / 1e+09, 2), " GB")
  }else if(ndigits > 5){
    value <- paste0(round(bytes / 1e+06, 2), " MB")
  }else if(ndigits > 2){
    value <- paste0(round(bytes / 1e+03, 2), " KB")
  }else{
    value <- paste0(round(bytes, 2), " B")
  }
  
  # Return configured value
  return(value)
  
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

#%%%%%%%%%%%%%%%%%%%%%%
# OBJECT FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Gets base `EGA` object that are necessary for most
# functions in {EGAnet}
# Updated 21.07.2023
get_EGA_object <- function(object)
{
  
  # `dynEGA` is a special case
  if(is(object, "dynEGA")){
    
    # Regardless of level, return all
    return(
      list(
        population = object$dynEGA$population,
        group = object$dynEGA$group,
        individual = object$dynEGA$individual
      )
    )
    
  }
  
  # `hierEGA` is also a special case
  if(is(object, "hierEGA")){
    
    # Return levels
    return(
      list(
        lower_order = object$lower_order,
        higher_order = object$higher_order
      )
    )
    
  }
  
  # `bootEGA` has some precedence with
  # empirical `EGA` object nested
  if(is(object, "bootEGA")){
    object <- object$EGA
  }
  
  # Based on class, get EGA object
  return(
    switch(
      tolower(class(object)),
      "ega" = object,
      "ega.fit" = object$EGA,
      "riega" = object$EGA,
      "hierega" = object,
      stop(
        "Could not find `EGA` object nested in object",
        call. = FALSE
      )
    )
  )
  
}

#' @noRd
# Get appropriate `dynEGA` objects ----
# Updated 09.07.2023
get_dynEGA_object <- function(dynEGA.object)
{
  
  # Check for legacy
  if(is(dynEGA.object, "dynEGA.ind.pop")){
    
    # Set up objects
    dynEGA.pop <- dynEGA.object$dynEGA.pop
    dynEGA.ind <- dynEGA.object$dynEGA.ind
    
    # Remove Methods
    if("Methods" %in% names(dynEGA.ind$dynEGA)){
      dynEGA.ind$dynEGA <- dynEGA.ind$dynEGA[names(dynEGA.ind$dynEGA) != "Methods"]
    }
    
    # Return proper objects
    return(
      list(
        population = dynEGA.pop$dynEGA,
        individual = dynEGA.ind$dynEGA
      )
    )
    
  }
  
  # Get EGA objects
  ega_objects <- get_EGA_object(dynEGA.object)
  
  # Determine which names are available
  available_objects <- c("population", "individual") %in% names(ega_objects)
  names(available_objects) <- c("population", "individual")
  
  # Get proper objects
  if(all(available_objects)){
    
    # Return population and individual
    return(
      list(
        population = ega_objects$population,
        individual = ega_objects$individual
      )
    )
    
  }else{ # Otherwise, send an error
    
    # Set up proper messaging
    missing_objects <- swiftelse(
      sum(available_objects) == 0,
      paste0("\"population\" and \"individual\""),
      paste0("\"", names(available_objects)[!available_objects], "\"")
    )
    
    # Send error
    stop(
      paste0(
        "The input for 'dynEGA.object' was class \"dynEGA\" but did not contain object(s): ",
        missing_objects, "\n\n",
        "To avoid this error, perform `dynEGA` using the argument: level = c(\"population\", \"individual\")"
      ),
      call. = FALSE
    )
    
  }
  
}

#' @noRd
# Faster data frame initialization ----
# Initializes a matrix and converts to a data frame
# Data frames are initialized as individual vectors
# carrying extra overhead 
# Updated 06.07.2023
fast.data.frame <- function(
    data = NA, nrow = 1, ncol = 1,
    rownames = NULL, colnames = NULL,
    ...
)
{
  
  return(
    as.data.frame(
      matrix(
        data = data,
        nrow = nrow, ncol = ncol,
        dimnames = list(rownames, colnames)
      ), make.names = FALSE,
      ...
    )
  )
  
}

#' @noRd
# Determine object type ----
# Updated 05.07.2023
get_object_type <- function(object)
{
  
  # Remove attributes and class from object
  unclassed_object <- remove_attributes(object)
  
  # Add object type to object's attributes
  ## In order of precedence
  if(is(object, "tbl")){ # needs precedence over 'data.frame'
    return("tibble")
  }else if(is.data.frame(object)){ # needs precedence over 'list'
    return("data.frame")
  }else if(is.list(unclassed_object)){ # needs precedence over 'vector'
    return("list")
  }else if(is.factor(object)){
    return("factor")
  }else if(is.vector(unclassed_object)){
    return("vector")
  }else if(is.matrix(unclassed_object)){
    return("matrix")
  }else if(is.array(unclassed_object)){
    return("array")
  }
  
}

#' @noRd
# Force vector ----
# For 1D matrix or data frames only
# Updated 05.07.2023
force_vector <- function(object)
{
  
  # Get object type
  object_type <- get_object_type(object)
  
  # Branch across types
  if(object_type == "vector"){
    return(object) # already a vector
  }else if(object_type %in% c("matrix", "data.frame", "tibble")){
    
    # Ensure matrix
    new_matrix <- as.matrix(object)
    
    # Get dimensions
    dimensions <- dim(new_matrix)
    
    # Check for more than one dimension
    if(all(dimensions) > 1){
      stop("Cannot sufficiently force object into 'vector'.", call. = FALSE)
    }
    
    # Set as vector
    new_vector <- as.vector(new_matrix)
    
    # Set names based on dimension
    names(new_vector) <- unlist(dimnames(object)[dimensions != 1])
    
    # Return new vector
    return(new_vector)
    
  }
  
}

#' @noRd
# Force numeric ----
# For factors
# Updated 05.07.2023
force_numeric <- function(object)
{
  
  # Check for whether object is already numeric
  if(is.numeric(object)){
    return(object) # just return
  }else if(is.factor(object)){ # factors are trickier...
    
    # Switch based on type in the levels
    if(is.numeric(levels(object))){
      
      # Get original values rather than factor values
      return(as.numeric(as.character(object)))
      
    }else{
      
      # Just used factor assigned values
      return(as.numeric(object))
      
    }
    
  }
  
}

#%%%%%%%%%%%%%%%%%%%%%
# COUNT FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Unique Length ----
# Counts length of unique variables (without `NA`)
# Updated 03.07.2023
unique_length <- function(x)
{
  return(length(unique(x)) - anyNA(x))
}

#' @noRd
# Edge count ----
# Counts the number of edges
# (assumes network is _symmetric_ matrix)
# Updated 11.10.2023
edge_count <- function(network, nodes, diagonal = FALSE)
{
  return((sum(network != 0, na.rm = TRUE) - swiftelse(diagonal, nodes, 0)) * 0.50)
}

#' @noRd
# Faster table ----
# Skips many checks in `table`
# Updated 22.07.2023
fast_table <- function(values)
{
  
  # Get factor values
  factor_values <- swiftelse(
    is.factor(values), values, factor(values)
  )
  
  # Get factor levels
  factor_levels <- levels(factor_values)
  
  # Return with names
  return(
    setNames(
      tabulate(factor_values, length(factor_levels)),
      factor_levels
    )
  )

}

#' @noRd
# Converts entire vectors to a single value ----
# Same logic as {plyr}'s `id_var` and `ninteraction`
# Avoids {plyr} dependency for single function
# Updated 06.07.2023
vector2factor <- function(data)
{
  
  # Use "paste" method
  stringed_values <- do.call( 
    # Get total unique values each variable takes
    paste, rev(lapply( # `ninteraction`
      
      # Ensure data frame
      as.data.frame(data), function(x){
      
      # Get matched ID (`id_var`)
      matched_id <- match(x, sort(unique(x), na.last = TRUE))
      
      # Attached number of unique as attribute
      attr(matched_id, "n") <- max(matched_id)
      
      # Return result
      return(matched_id)
      
    }))
  )

  # Return unique ID
  return(factor(match(stringed_values, unique(stringed_values))))
  
}

#' @noRd
# Count table ----
# Provides counts of repeated rows in a data frame
# Same logic as {plyr}'s `count`
# Avoids {plyr} dependency for single function
# Updated 22.07.2023
count_table <- function(data, proportion = FALSE)
{

  # Ensure data frame
  data <- as.data.frame(data)
  
  # Get unique ID
  unique_id <- vector2factor(data)
  
  # Get unique solutions
  unique_solutions <- data[!duplicated(unique_id),,drop = FALSE]
  
  # Tabulate data
  Value <- fast_table(unique_id) # tabulate(unique_id, length(levels(unique_id)))
  
  # Check for proportion
  if(proportion){
    Value <- Value / dim(data)[1]
  }
  
  # Return data frame
  return(
    as.data.frame(
      cbind(unique_solutions[order(unique(unique_id)),], Value)
    )
  )
  
}

#%%%%%%%%%%%%%%%%%%%%
# DATA FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Determine number of categories in data ----
# Updated 06.07.2023
data_categories <- function(data)
{
  return(nvapply(as.data.frame(data), unique_length))
}

#' @noRd
# All-purpose symmetric checker ----
# Faster but less robust than `isSymmetric`
# Defaults to floating point tolerance
# Updated 01.07.2023
is_symmetric <- function(data, tolerance = 1e-06){

  # Get dimensions
  dimensions <- dim(data)
  
  # Check for whether rows equal columns
  if(dimensions[1] != dimensions[2]){ # Not a (square) matrix
    return(FALSE)
  }else{
    
    # Convert to matrix
    data_matrix <- as.matrix(data)
    
    # Check that all are equal
    return(all(abs(data_matrix - t(data_matrix)) < tolerance))
    
  }
  
}

#' @noRd
# Ensure data has dimension names ----
# Updated 24.07.2023
ensure_dimension_names <- function(data)
{
  
  # Determine object type
  object_type <- get_object_type(data)
  
  # Branch for vector vs. matrix and data frame
  if(object_type == "vector"){
    
    # Get length of vector
    vector_length <- length(data)
    
    # Check for names
    if(is.null(names(data))){
      names(data) <- paste0(
        "V", format_integer(
          seq_len(vector_length),
          digits(vector_length) - 1
        )
      )
    }
    
  }else if(object_type %in% c("matrix", "data.frame")){ 
    
    # Get dimensions
    dimensions <- dim(data)
    
    # Get dimension names
    dimension_names <- dimnames(data)
    
    # Check for column names
    if(is.null(dimension_names[[2]])){
      
      # Standardize names
      dimnames(data)[[2]] <- paste0(
        "V", format_integer(
          seq_len(dimensions[2]),
          digits(dimensions[2]) - 1
        )
      )
      
    }
    
    # Check for symmetric matrix
    if(is_symmetric(data)){
      
      # Check for row names
      if(is.null(dimnames(data)[[1]]) | any(dimension_names[[1]] != dimension_names[[2]])){
        
        # Assign column names to row names
        dimnames(data)[[1]] <- dimnames(data)[[2]]
        
      }
      
    }
    
  }
  
  # Return named data
  return(data)
  
}

#' @noRd
# Transfer names from data to output ----
# Usually for matrix and data frame
# Updated 25.06.2023
transfer_names <- function(data, output)
{
  
  # Obtain column names
  column_names <- dimnames(data)[[2]]
  
  # Check for variable names
  if(!is.null(column_names)){
    
    # Add names to rows and columns
    dimnames(output) <- list(column_names, column_names)
    
  }
  
  # Return named output
  return(output)
  
}

#' @noRd
# Function to check for usable variables ----
# Updated 24.07.2023
usable_data <- function(data, verbose)
{
  
  # Turn data into a data matrix
  data_matrix <- data.matrix(data)
  
  # All missing after coercions
  remove_columns <- lvapply(as.data.frame(is.na(data_matrix)), all)
  
  # Send warning if there are any removals
  if(any(remove_columns)){
    
    # Send warning
    warning(
      paste0(
        "Some variables could not to be coerced to numeric values. These variables have been removed from the analysis:\n",
        paste0("'", dimnames(data_matrix)[[2]][remove_columns], "'", collapse = ", "),
        "\n\nIf these variables were not intended to be removed, then try converting them to numeric values before inputting the data into the function"
      ), call. = FALSE
    )
    
    # Remove these variables from `data` and `data_matrix`
    keep_columns <- !remove_columns
    data <- data[,keep_columns]
    data_matrix <- data_matrix[,keep_columns]
    
  }
  
  # Determine whether each variable was coerced
  coercions <- lvapply(as.data.frame(data_matrix != data), any, na.rm = TRUE)

  # Send warning for any coercions
  if(any(coercions) & verbose){
    
    # Send warning
    warning(
      paste0(
        "Several variables were coerced to numeric values. These variables were changed to numeric values:\n",
        paste0("'", dimnames(data_matrix)[[2]][coercions], "'", collapse = ", ")
      ), call. = FALSE
    )
    
  }
  
  # Send back usable data
  return(data_matrix)
  
}

#' @noRd
# Obtain {igraph} and matrix networks ----
# Updated 24.07.2023
obtain_networks <- function(network, signed)
{
  
  # Get network matrix first
  if(is(network, "igraph")){
    network_matrix <- igraph2matrix(network)
  }else{ # If not {igraph}, then ensure matrix
    network_matrix <- as.matrix(network)
  }
  
  # Then, check for sign
  if(!signed){
    network_matrix <- abs(network_matrix)
  }

  # Return networks
  return(
    list(
      # (Re-)convert to {igraph} network
      igraph_network = convert2igraph(network_matrix),
      network_matrix = network_matrix
    )
  )
  
}

#' @noRd
# Whether usable data is needed ----
# Updated 07.09.2023
needs_usable <- function(ellipse)
{
  return(!"needs_usable" %in% names(ellipse) || ellipse$needs_usable)
}

#' @noRd
# Obtain data, sample size, correlation matrix ----
# Generic function to get the usual needed inputs
# Updated 04.09.2023
obtain_sample_correlations <- function(data, n, corr, na.data, verbose, ...)
{
  
  # Check if data is a correlation matrix
  if(is_symmetric(data)){
    
    # Check for sample size
    if(is.null(n)){
      stop(
        "A symmetric matrix was provided in the 'data' argument but the sample size argument, 'n', was not set. Please input the sample size into the 'n' argument.",
        call. = FALSE
      )
    }
    
    # If symmetric and sample size is provided, then
    # data to correlation matrix
    correlation_matrix <- data
    
  }else{
    
    # Check for usable data
    if(needs_usable(list(...))){
      data <- usable_data(data, verbose)
    }
    
    # Obtain sample size
    n <- dim(data)[1]
    
    # Check for automatic correlations
    if(corr == "auto"){
      correlation_matrix <- auto.correlate(
        data = data, corr = "pearson", na.data = na.data, 
        verbose = verbose, ...
      )
    }else if(corr == "cor_auto"){
      
      # Get arguments for `cor_auto`
      cor_auto_ARGS <- obtain_arguments(
        FUN = qgraph::cor_auto,
        FUN.args = list(...)
      )
      
      # Set 'data' and 'verbose' arguments
      cor_auto_ARGS[c("data", "verbose")] <- list(data, verbose)
      
      # Obtain correlations
      correlation_matrix <- do.call(
        what = qgraph::cor_auto,
        args = cor_auto_ARGS
      )
      
    }else{
      correlation_matrix <- cor(data, use = na.data, method = corr)
    }
    
  }
  
  # Return results
  return(
    list(
      data = data, n = n,
      correlation_matrix = correlation_matrix
    )
  )
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%
# QUIET CALLS & LOADS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%

# Lots of R packages and functions print
# messages that add to the noise of {EGAnet}
# outputs; these functions make other
# packages' messages go away
#
# NOTE: In some cases, some package messages
# are important and should not be ignored
# (don't use these functions in those cases)

#' @noRd
# General function to silently obtain output ----
# Updated 01.07.2023
silent_call <- function(...)
{
  
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
# General function to silently load package ----
# Updated 01.07.2023
silent_load <- function(...)
{
  return(suppressPackageStartupMessages(...))
}

#' @noRd
# Specific function to silently plot ----
# Updated 06.07.2023
silent_plot <- function(...)
{
  return(silent_call(silent_load(plot(...))))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FUNCTIONS & ARGUMENTS FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Strip classes and attributes
# Updated 09.07.2023
remove_attributes <- function(object)
{

  # Remove all except special attributes (see `?structure`)
  attributes(object) <- attributes(object)[
    names(attributes(object)) %in% c(
      "dim", "dimnames", "names", "tsp", "levels"
    )
  ]
  
  # Return unclassed object
  return(unclass(object))
  
}

#' @noRd
# Set default argument (cleaner missing function) ----
# Updated 28.07.2023
set_default <- function(argument, default, FUN, several.ok = FALSE)
{

  # Check for type error
  typeof_error(argument, c(typeof(default), "closure"), "set_default")
  
  # Return argument if it's a function
  if(is.function(argument)){
    return(argument)
  }
  
  # Return default if NULL
  if(is.null(argument)){
    return(default)
  }
  
  # Return default if length > 1 and several are NOT okay
  if(length(argument) > 1 & !several.ok){
    return(default)
  }
  
  # Get argument name
  argument_name <- as.character(substitute(argument))
  
  # Get choices for the function
  if(!is.function(FUN)){ # choices are provided in 'FUN'
    choices <- FUN
  }else{ # get formal arguments from the function itself
    
    # Get formal argument choices from default call
    choices <- formals(FUN)[[argument_name]]
    
    # Check if choices is a call
    if(is.call(choices)){
      
      # Get choices
      choices <- as.character(choices)
      
      # Remove "c" from call
      choices <- choices[-which(choices == "c")]
      
    }
    
  }
  
  # Make argument lowercase (if necessary)
  argument <- tolower(argument); choices <- tolower(choices)
  
  # Get arguments not in choices
  if(any(!argument %in% choices)){ # Repeated calls OK since returning error
    
    # Get arguments not in choices
    not_in_choices <- !argument %in% choices
    
    # Get number not in choices
    number_not <- sum(not_in_choices)
    
    # Send error for bad argument
    if(number_not > 0){
      
      # Only one not in choices
      if(number_not == 1){
        argument_text <- paste0(" = \"", argument[not_in_choices], "\"")
      }else{ # Multiple not in choices
        argument_text <- paste0(
          " = c(",
          paste0("\"", argument[not_in_choices], "\"", collapse = ", "),
          ")"
        )
      }
      
      # Throw error
      stop(
        paste0(
          paste0("Invalid argument: ", argument_name, argument_text),
          "\n\nPlease choose from: ", 
          paste0("\"", choices, "\"", collapse = ", ")
        ), call. = FALSE
      )
      
    }
    
  }
  
  # Return argument
  return(argument)
  
}

#' @noRd
# Overwrite Arguments ----
# Updated 24.06.2023
overwrite_arguments <- function(main, ARGS)
{
  
  # Obtain replacement arguments
  target_args <- ARGS[names(ARGS) %in% names(main)]
  
  # Replace the arguments with the same name
  main[names(target_args)] <- target_args
  
  # Return main arguments
  return(main)
  
}

#' @noRd
# Function to obtain arguments ----
# Updated 06.07.2023
obtain_arguments <- function(FUN, FUN.args)
{

  # Overwrite arguments
  FUN.formals <- overwrite_arguments(formals(FUN), FUN.args)
  
  # Remove ellipses
  FUN.formals <- FUN.formals[names(FUN.formals) != "..."]
  
  # Remove call arguments (assume they are supplied elsewhere)
  return(FUN.formals[!lvapply(FUN.formals, is, "call")])
  
}

#' @noRd
# Legacy Argument Overwrite ----
# Preference given to legacy
# Updated 25.07.2023
legacy_overwrite <- function(ellipse, legacy_ARG)
{

  # Check for legacy argument in ellipse
  if(legacy_ARG %in% names(ellipse)){
    
    # Get legacy arguments
    legacy_ARGS <- ellipse[[legacy_ARG]]
    
    # Overwrite arguments
    ellipse[names(legacy_ARGS)] <- legacy_ARGS
    
    # Remove legacy argument
    ellipse <- ellipse[names(ellipse) != legacy_ARG]
    
  }
  
  # Return ellipse
  return(ellipse)
  
}

#' @noRd
# Legacy Argument Handling ----
# Updated 24.06.2023
legacy_EGA_args <- function(ellipse)
{
  
  # `model.args`
  ellipse <- legacy_overwrite(ellipse, "model.args")
  
  # `algorithm.args`
  ellipse <- legacy_overwrite(ellipse, "algorithm.args")
  
  # `plot.args`
  ellipse <- legacy_overwrite(ellipse, "plot.args")
  
  # Return ellipse
  return(ellipse)
  
}

#' @noRd
# Make unidimensional CFA model ----
# Updated 08.11.2023
make_unidimensional_cfa <- function(variable_names)
{
  return(
    paste(
      "LF =~",
      swiftelse(
        length(variable_names) == 2,
        paste0("a*", variable_names, collapse = " + "),
        paste(variable_names, collapse = " + ")
      )
    )
  )
}

#' @noRd
# Determine estimator arguments ----
# Updated 24.07.2023
estimator_arguments <- function(lavaan_ARGS, ellipse)
{
  
  # Check for `ordinal.categories` in ellipse arguments
  ordinal.categories <- swiftelse(
    !"ordinal.categories" %in% names(ellipse), 7, # default
    ellipse$ordinal.categories
  )
  
  # Obtain categories
  categories <- data_categories(lavaan_ARGS$data)
  
  # Get categorical variables
  categorical_variables <- categories <= ordinal.categories
  
  # Check for categories
  if(any(categorical_variables)){
    
    # Set arguments
    lavaan_ARGS$estimator <- "WLSMV"
    lavaan_ARGS$missing <- "pairwise"
    lavaan_ARGS$ordered <- names(categories)[categorical_variables]
    
  }else{
    
    # Set arguments
    lavaan_ARGS$estimator <- "MLR"
    lavaan_ARGS$missing <- "fiml"
    lavaan_ARGS$ordered <- FALSE
    
  }
  
  # Return arguments
  return(lavaan_ARGS)
  
  
}

#%%%%%%%%%%%%%%%%%%%%
# PLOT FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Defaults for GGally plotting ----
# For plots and methods
# Updated 04.08.2023
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
    label.color = "black", mode = "qgraph",
    edge.label.color = "black", node.alpha = 0.50,
    node.shape = 19, node.size = 12, 
    edge.alpha = "edge.alpha", edge.size = 8
  )
  
  # Legacy arguments
  ellipse <- legacy_EGA_args(ellipse)
  
  # Replace `ggnet2` arguments with {EGAnet} arguments
  default_args <- overwrite_arguments(default_args, ega_default_args)
  
  # Replace `ggnet2` arguments with arguments input
  default_args <- overwrite_arguments(default_args, ellipse)
  
  # Remove the ellipse
  default_args <- default_args[names(default_args) != "..."]
  
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
    if(any(tolower(ellipse$color.palette) %in% gray_options)){
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
# Updated 13.08.2023
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
  typeof_error(plot_ARGS$label.alpha, "numeric", "plot.EGAnet")
  length_error(plot_ARGS$label.alpha, c(1, nodes), "plot.EGAnet")
  
  # Node Label Color
  typeof_error(plot_ARGS$label.color, "character", "plot.EGAnet")
  length_error(plot_ARGS$label.color, c(1, nodes), "plot.EGAnet")
  
  # Node Label Size
  typeof_error(plot_ARGS$label.size, "numeric", "plot.EGAnet")
  length_error(plot_ARGS$label.size, c(1, nodes), "plot.EGAnet")
  
  # Node Label
  typeof_error(plot_ARGS$node.label, "character", "plot.EGAnet")
  length_error(plot_ARGS$node.label, c(1, nodes), "plot.EGAnet")
  
  # Node Alpha
  typeof_error(plot_ARGS$node.alpha, "numeric", "plot.EGAnet")
  length_error(plot_ARGS$node.shape, c(1, communities, nodes), "plot.EGAnet")
  
  # Node Color
  typeof_error(plot_ARGS$node.color, "character", "plot.EGAnet")
  length_error(plot_ARGS$node.color, c(1, communities, nodes), "plot.EGAnet")
  
  # Node Shape
  typeof_error(plot_ARGS$node.shape, "numeric", "plot.EGAnet")
  length_error(plot_ARGS$node.shape, c(1, communities, nodes), "plot.EGAnet")
  
  # Node Size
  typeof_error(plot_ARGS$node.size, "numeric", "plot.EGAnet")
  length_error(plot_ARGS$node.size, c(1, communities, nodes), "plot.EGAnet")
  
  ### Edge arguments
  
  # Edge Alpha
  typeof_error(plot_ARGS$edge.alpha, "numeric", "plot.EGAnet")
  length_error(plot_ARGS$edge.alpha, c(1, edges), "plot.EGAnet")
  
  # Edge Color (allow two for positive and negative)
  typeof_error(plot_ARGS$edge.color, "character", "plot.EGAnet")
  length_error(plot_ARGS$edge.color, c(1, 2, edges), "plot.EGAnet")
  
  # Edge Size
  typeof_error(plot_ARGS$edge.size, "numeric", "plot.EGAnet")
  length_error(plot_ARGS$edge.size, c(1, edges), "plot.EGAnet")
  
  # Edge line type (allow two for positive and negative)
  typeof_error(plot_ARGS$edge.lty, "character", "plot.EGAnet")
  length_error(plot_ARGS$edge.lty, c(1, 2, edges), "plot.EGAnet")
  
}

#' @noRd
# Re-scale edges ----
# Updated 01.07.2023
rescale_edges <- function(network, edge_size)
{
  
  # Set up scaling sequence 
  scale_sequence <- seq.int(0, 1, 0.0001)
  names(scale_sequence) <- scale_sequence
  
  # Set edge scaling (default `edge.size = 8`)
  edge_scaling <- scale_sequence * edge_size
  
  # Return scaled edges
  return(unname(edge_scaling[as.character(abs(network))]))
  
}

#' @noRd
# Readable names ----
# Updated 01.07.2023
readable_names <- function(node_names)
{
  
  # Add return to names
  return(
    cvapply(
      strsplit(node_names, split = " "), function(x){
        
        # Obtain words in name
        words <- length(x)
        
        # Determine if split is necessary
        if(words > 1){
          
          # Determine number of lines
          add_line <- round(words / 2)
          
          # Paste back together name
          name <- paste(
            paste(x[seq_len(add_line)], collapse = " "),
            paste(x[(add_line + 1):words], collapse = " "),
            sep = "\n"
          )
          
          # Return name
          return(name)
          
        }else{return(x)}
        
      }
    )
  )
  
}

#' @noRd
# Get network layout ----
# Updated 04.08.2023
get_layout <- function(network, dimensions, non_zero_index, plot_ARGS)
{
  
  # Determine whether "mode" was used
  if(is.character(plot_ARGS$mode)){
    
    # Check for {qgraph}
    if(plot_ARGS$mode == "qgraph"){
      
      # Lower triangle for edge list
      network_lower <- network[lower.tri(network)]
      weights_lower <- abs(network_lower[network_lower != 0])
      
      # Set up edge list
      edge_list <- which(non_zero_index, arr.ind = TRUE)
      edge_list <- edge_list[edge_list[,"row"] < edge_list[,"col"],, drop = FALSE]
      
      # Set layout (spring)
      network_layout <- qgraph::qgraph.layout.fruchtermanreingold(
        edgelist = edge_list[order(edge_list[,"row"]),, drop = FALSE],
        weights = (weights_lower / max(weights_lower))^2,
        vcount = dimensions[2]
      )
      
    }else{
      
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
      
    }
    
  }else{ # Assume "mode" is a 2D matrix corresponding to a layout
    network_layout <- plot_ARGS$mode
  }

  # Return layout
  return(network_layout)
  
}

#' @noRd
# Basic set up for plots ----
# Updated 20.08.2023
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
  communities <- unique_length(wc)
  
  # Obtain node names
  node_names <- dimnames(network)[[2]]

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
    non_zero_index, plot_ARGS
  )
  
  ### Generic arguments (mostly handled in `GGally_args`)
    
  ## Remove some arguments
  plot_ARGS[c("alpha", "color", "size")] <- NULL
  
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
    plot_ARGS$edge.color <- swiftelse(non_zero_edges >= 0, plot_ARGS$edge.color[1], plot_ARGS$edge.color[2])
  }
  
  ## Set edge line type
  if(length(plot_ARGS$edge.lty) == 2){
    plot_ARGS$edge.lty <- swiftelse(non_zero_edges >= 0, plot_ARGS$edge.lty[1], plot_ARGS$edge.lty[2])
  }
  
  ## Set edge size (scale by `edge.size`)
  if(length(plot_ARGS$edge.size) == 1){
    plot_ARGS$edge.size <- rescale_edges(non_zero_edges, plot_ARGS$edge.size)
  }
  
  ## Edge label size (not used)
  plot_ARGS$edge.label.size <- swiftelse(
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
  
  # Set up node names to be more readable
  node_names <- readable_names(plot_ARGS$node.label)
  
  # Remove node labels
  plot_ARGS$node.label <- NULL
  
  # Get first layer with silent call
  first_layer <- silent_call(
    do.call(GGally::ggnet2, plot_ARGS)
  )
  
  # Return node size to `plot_ARGS` (was removed above)
  plot_ARGS$node.size <- node.size

  # Determine border color
  ## Check for gray scale options
  gray_options <- c(
    "greyscale", "grayscale", "colorblind"
  )
  
  ## Set border color
  if(all(is.na(wc))){ # Plain network (without communities)
    border_color <- rep("black", dimensions[2])
  }else if( 
    length(color.palette) == 1 &&
    color.palette %in% gray_options
  ){ # Gray scale network
    border_color <- swiftelse(palette == "white", "white", "black")
  }else{ # Same color as nodes
    border_color <- plot_ARGS$node.color
  }
  
  # Custom nodes: transparent insides and dark borders
  second_layer <- first_layer +
    ggplot2::geom_point( # transparent insides
      size = node.size + 0.50, shape = 19,
      color = plot_ARGS$node.color,
      alpha = plot_ARGS$node.alpha,
      show.legend = FALSE
    ) +
    ggplot2::geom_point( # dark borders
      size = node.size, color = border_color,
      shape = 1, stroke = 1.5, alpha = 0.80
    ) +
    ggplot2::geom_text( # put text back on top
      ggplot2::aes(label = node_names), 
      color = plot_ARGS$label.color,
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
        title = swiftelse(
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
    return(
      list(
        network_plot = second_layer,
        ARGS = plot_ARGS
      )
    )
    
  }else{
    
    # Return plot only
    return(second_layer)
    
  }
  
}

#' @noRd
# Basic set up for single plot ----
# Updated 25.07.2023
single_plot <- function(network, wc = NULL, ...)
{
  
  # Look for memberships in arguments
  ## If no memberships, then plot network
  # as if all memberships are missing
  if(is.null(wc)){
    wc <- rep(NA, dim(network)[2])
  }
  
  # Reorder network and communities
  new_order <- order(wc)
  
  # Send on and return from `basic_plot_setup`
  return(
    basic_plot_setup(
      network[new_order, new_order], 
      wc[new_order], 
      ...
    )
  )
  
  # I'm not sure why the memberships were ordered
  # in the original code (before refactoring)
  # `single_plot` doesn't affect anything in
  # terms of reordering but it's not good
  # for more than one plot
  #
  # `basic_plot_setup` is preferred for multiple plots
  
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
      ),
      call. = FALSE
    )
    
  }
  
  # Get names
  original_names <- dimnames(original)[[2]]
  comparison_names <- dimnames(comparison)[[2]]
  
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
      ),
      call. = FALSE
    )
    
  }
  
}

#' @noRd
# Basic set up for comparing plots ----
# Updated 26.10.2023
compare_plots <- function(comparison_network, comparison_wc, plot_ARGS)
{
  
  # original network = plot_ARGS$net

  # Make sure dimensions are the same before proceeding
  dimension_comparison(plot_ARGS$net, comparison_network)
  
  # Get comparison network names
  comparison_names <- dimnames(comparison_network)[[2]]
  
  # Ensure row names to ensure proper ordering
  dimnames(comparison_network)[[1]] <- comparison_names
  
  # Put network into same order as original network
  matching_order <- match(
    dimnames(plot_ARGS$net)[[2]], # target to match
    comparison_names # adjust to target
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
  plot_ARGS[c("net", "node.color")] <- NULL
  
  # Check for edges
  ## Assume more than one edge alpha is default
  if(length(plot_ARGS[["edge.alpha"]]) > 1){
    plot_ARGS[["edge.alpha"]] <- NULL
  }
  
  ## Assume more than two edge color is default
  if(length(plot_ARGS[["edge.color"]]) > 2){
    plot_ARGS[["edge.color"]] <- NULL
  }
  
  ## Assume more than two edge line type is default
  if(length(plot_ARGS[["edge.lty"]]) > 2){
    plot_ARGS[["edge.lty"]] <- NULL
  }
  
  ## Assume more than one edge size is default
  if(length(plot_ARGS[["edge.size"]]) > 1){
    plot_ARGS[["edge.size"]] <- NULL
  }
  
  # Send on and return from `basic_plot_setup`
  return(do.call(basic_plot_setup, plot_ARGS))
  
}

#%%%%%%%%%%%%%%%%%%%%%
# ERROR FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Experimental warning ----
# Updated 03.08.2023
experimental <- function(function_name)
{
  warning(
    paste0(
      "This implementation of `", function_name,
      "` is ", styletext("experimental", "italics"), ". \n\nThe underlying function ",
      "and/or output may change until the results have been appropriately vetted and validated."
    ), call. = FALSE
  )
}

#' @noRd
# Not refactored warning ----
# Updated 04.08.2023
not_refactored <- function(function_name)
{
  warning(
    paste0(
      "This implementation of `", function_name, "` was not refactored in the {EGAnet} 2.0.0 update. ",
      "There are no guarantees that this function will work.\n\n",
      "Please do not create a GitHub issue or bug report for this function"
    ), call. = FALSE
  )
}

#' @noRd
# Error for class ----
# Updated 04.09.2023
class_error <- function(input, expected_class, function_name){
  
  # Check for object types
  if(!is(input, expected_class)){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)),
        "' is an object with ", paste0("'", class(input), "'", collapse = ", "), " class(es)",
        ". Input is expected to be ", paste0("'", expected_class, "'", collapse = ", "), " class(es)",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#class-error"
      ), 
      call = function_name
    )
  }
  
}

#' @noRd
# Error for object type ----
# Updated 04.09.2023
object_error <- function(input, expected_type, function_name){
  
  # Get input type
  input_type <- get_object_type(input)
  
  # Check for object types
  if(!input_type %in% expected_type){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)),
        "' is a ", paste0("'", input_type, "'", collapse = ", "), " object",
        ". Input is expected to be ", paste0("'", expected_type, "'", collapse = ", "), " object",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#object-error"
      ), 
      call = function_name
    )
  }
  
}

#' @noRd
# Error for `typeof` ----
# Updated 04.09.2023
typeof_error <- function(input, expected_value, function_name){
  
  # Switch out "closure" with "function"
  if("closure" %in% expected_value){
    expected_value[expected_value == "closure"] <- "function"
  }
  
  # Get type of input
  typeof_input <- typeof(input)
  
  # Convert "closure" to "function"
  if("closure" %in% typeof_input){
    typeof_input[typeof_input == "closure"] <- "function"
  }
  
  # Convert "integer" and "double" to "numeric"
  ## Input
  typeof_input <- swiftelse(
    typeof_input %in% c("integer", "double"),
    "numeric", typeof_input
  )
  ## Expected value
  expected_value <- swiftelse(
    any(expected_value %in% c("integer", "double")),
    "numeric", expected_value
  )
  
  # Check for value
  if(!typeof_input %in% expected_value){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Input into '", deparse(substitute(input)),
        "' is ", paste("'", typeof_input, "'", sep = "", collapse = ", "), " type",
        ". Input is expected to be ",
        paste0("'", expected_value, "'", collapse = " or "), " type",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#typeof-error"
        # can use "or" because `typeof` won't ever be more than two
      ), 
      call = function_name
    )
  }
  
}

#' @noRd
# Error for `length` ----
# Updated 04.09.2023
length_error <- function(input, expected_lengths, function_name){
  
  # Check for length of input in expected length
  if(!length(input) %in% expected_lengths){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Length of '", deparse(substitute(input)),
        "' (", length(input),") does not match expected length(s). Length must be: ",
        paste0("'", expected_lengths, "'", collapse = " or "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#length-error"
      ), 
      call = function_name
    )
  }
  
}

#' @noRd
# Error for `range` ----
# Updated 04.09.2023
range_error <- function(input, expected_ranges, function_name){
  
  # Obtain expected maximum and minimum values
  expected_maximum <- max(expected_ranges, na.rm = TRUE)
  expected_minimum <- min(expected_ranges, na.rm = TRUE)
  
  # Obtain maximum and minimum values
  actual_maximum <- round(max(input, na.rm = TRUE), 3)
  actual_minimum <- round(min(input, na.rm = TRUE), 3)
  
  # Check for maximum of input in expected range
  if(actual_maximum > expected_maximum){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Maximum value of '", deparse(substitute(input)),
        "' (", actual_maximum,") does not match expected range(s). Values must range between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#range-error"
      ), 
      call = function_name
    )
  }
  
  # Check for maximum of input in expected range
  if(actual_minimum < expected_minimum){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Minimum value of '", deparse(substitute(input)),
        "' (", actual_minimum,") does not match expected range(s). Values must range between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#range-error"
      ), 
      call = function_name
    )
  }
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COVARIANCE CONVERSION FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Zero-order correlations to partial correlations ----
# Updated 29.06.2023
cor2pcor <- function(correlation_matrix)
{
  
  # Convert to inverse correlations to partial correlations
  partial_correlations <- -cov2cor(solve(correlation_matrix))
  
  # Set diagonal to zero
  diag(partial_correlations) <- 0
  
  # Return
  return(partial_correlations)
  
}

#' @noRd
# Partial correlations to zero-order correlations ----
# Updated 29.06.2023
pcor2cor <- function(partial_correlations)
{
  
  # Set diagonal to negative 1
  diag(partial_correlations) <- -1
  
  # Return partial correlations as zero-order correlations
  return(cov2cor(solve(-partial_correlations)))
  
}

#' @noRd
# Partial correlations to inverse covariances ----
# Updated 29.06.2023
pcor2inv <- function(partial_correlations)
{
  
  # Set diagonal to negative 1
  diag(partial_correlations) <- -1
  
  # Return inverse covariance matrix
  return(solve(-partial_correlations))
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%
## NETWORK FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Create sparse network ----
# Updated 08.11.2023
sparse_network <- function(network)
{
  
  # Get number of nodes
  nodes <- dim(network)[2]
  
  # Get node sequence
  node_sequence <- seq_len(nodes)
  
  # Create data frame
  sparse_df <- data.frame(
    row = rep(node_sequence, each = nodes),
    col = rep(node_sequence, times = nodes),
    weight = as.vector(network)
  )
  
  # Return lower triangle
  return(sparse_df[sparse_df$row < sparse_df$col,])
  
}

#' @noRd
# Scramble networks ----
# Updated 10.11.2023
network_scramble <- function(base, comparison)
{
  
  # Get number of nodes
  nodes <- dim(base)[2]
  
  # Get sparse networks
  base_sparse <- sparse_network(base)
  comparison_sparse <- sparse_network(comparison)
  
  # Get edges
  base_edges <- base_sparse$weight != 0
  comparison_edges <- comparison_sparse$weight != 0
  
  # Get shared edges
  shared_total <- sum(base_edges & comparison_edges)
  
  # Decide on how many to add
  if(sum(base_edges) >= sum(comparison_edges)){
    
    # Assign edges
    base_sparse$weight[-shuffle(which(base_edges), shared_total)] <- 0
    base_sparse$weight[unique_index] <- comparison_sparse$weight[
      which(!base_edges & comparison_edges)
    ]
    
    # Set equivalent edges
    equivalent_sparse <- base_sparse
    
  }else{
    
    # Assign edges
    comparison_sparse$weight[-shuffle(which(comparison_edges), shared_total)] <- 0
    comparison_sparse$weight[unique_index] <- base_sparse$weight[
      which(base_edges & !comparison_edges)
    ]
    
    # Set equivalent edges
    equivalent_sparse <- comparison_sparse
    
  }
  
  # Remove zero edges from equivalent
  equivalent_sparse <- equivalent_sparse[
    equivalent_sparse$weight != 0,
  ]
  
  # Initialize network to return
  return_network <- matrix(0, nrow = nodes, ncol = nodes)
  
  # Loop over sparse equivalent
  for(i in nrow_sequence(equivalent_sparse)){
    
    # Set indices
    index <- equivalent_sparse[i,]
    
    # Populate return network
    return_network[index[,1], index[,2]] <-
      return_network[index[,2], index[,1]] <-
      index[,3]
    
  }
  
  # Return the network
  return(return_network)
  
}

#' @noRd
# Rewiring based on {igraph} ----
# Updated 30.10.2023
igraph_rewire <- function(network, prob, noise = 0)
{
  
  # Get nodes
  nodes <- dim(network)[2]
  
  # Assume NAs are zero
  network[is.na(network)] <- 0
  
  # Get rewired network
  rewired_network <- igraph2matrix(
    igraph::rewire(
      graph = convert2igraph(network),
      with = igraph::each_edge(prob = prob)
    )
  )
  
  # Add noise (if any)
  if(noise != 0){
    
    # Get absolute noise
    abs_noise <- abs(noise)
    
    # Get lower triangle
    lower_triangle <- lower.tri(rewired_network)
    
    # Get lower triangle
    rewired_lower <- rewired_network[lower_triangle]
    
    # Get non-zero
    non_zero <- rewired_lower != 0
    
    # Set noise
    rewired_lower[non_zero] <- rewired_lower[non_zero] +
      runif_xoshiro(sum(non_zero), min = -abs_noise, max = abs_noise)
    
    # Create new matrix
    rewired_network <- matrix(0, nrow = nodes, ncol = nodes)
    
    # Add to lower triangle
    rewired_network[lower_triangle] <- rewired_lower
    
    # Make symmetric
    rewired_network <- rewired_network + t(rewired_network)
    
  }
  
  # Return rewired network
  return(rewired_network)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MATH & STATS FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# There are some redundancies in the math functions because
# it's often faster to call some functions directly
# rather than nesting non-base R functions in functions

#' @noRd
# Matrix entropy calculation ----
# The default `exp(1)` base is pre-calculated to reduce function calls
# Returns entropy value
# Updated 29.06.2023
matrix_entropy <- function(density_matrix, base = 2.718282)
{
  return(
    -sum(
      diag(
        density_matrix %*% log(density_matrix, base = base)
      ), na.rm = TRUE
    )
  )

}

#' @noRd
# Basic entropy calculation ----
# Returns entropy value
# Updated 29.06.2023
entropy <- function(values, base = 2.718282)
{
  return(-sum(values * log(values, base = base), na.rm = TRUE))
}

#' @noRd
# Positive definite matrix ----
# Logical for whether a matrix is positive definite
# Updated 29.06.2023
is_positive_definite <- function(data)
{
  return(all(eigen(data, symmetric = TRUE, only.values = TRUE)$values > 0))
}

#' @noRd
# Wrapper for `eigen` ----
# Extracts eigenvalues only
# Updated 29.06.2023
matrix_eigenvalues <- function(data)
{
  return(eigen(data, symmetric = TRUE, only.values = TRUE)$values)
}

#' @noRd
# Trace of matrix ----
# Updated 25.06.2023
trace <- function(object)
{
  return(sum(diag(object), na.rm = TRUE))
}

#' @noRd
# Cohen's d ----
# Updated 13.07.2023
d <- function(sample1, sample2, paired = FALSE)
{
  
  # Check for paired
  if(isTRUE(paired)){
    
    # Get differences
    differences <- sample1 - sample2
  
    # Return paired Cohen's d
    return(
      mean(differences, na.rm = TRUE) /
      sd(differences, na.rm = TRUE)
    )
    
  }
  
  # Get usable indices
  usable1 <- sample1[!is.na(sample1)]
  usable2 <- sample2[!is.na(sample2)]
  
  # Numerator
  num <- mean(usable1) - mean(usable2)
  
  # Degrees of freedom
  df1 <- length(usable1) - 1
  df2 <- length(usable2) - 1
  
  # Denominator
  denom <- sqrt(
    (
      (df1 * var(usable1)) + (df2 * var(usable2))
    ) / (df1 + df2)
  )
  
  # Return Cohen's d
  return(abs(num / denom))
  
}

#' @noRd
# Adaptive Alpha ----
# Needs desparate updating
# Updated 01.08.2022
adapt.a <- function (test = c("anova","chisq","cor","one.sample","two.sample","paired"),
                     ref.n = NULL, n = NULL, alpha = .05, power = .80,
                     efxize = c("small","medium","large"), groups = NULL, df = NULL)
{
  
  # Need a test
  if(missing(test)){
    stop("test must be selected")
  }else{test <- match.arg(test)}
  
  # Assign medium effect size
  if(missing(efxize)){
    efxize <- "medium"
    message("No effect size selected. Medium effect size computed.")
  }else{efxize <- efxize}
  
  # ANOVA
  if(test == "anova"){
    
    # Check for groups
    if(is.null(groups)){
      stop("ANOVA is selected. Number of groups must be set")
    }
    
    # Set effect size
    efxize <- switch(
      efxize,
      "small" = 0.10,
      "medium" = 0.25,
      "large" = 0.40
    )
    
    # Determine reference sample size
    if(is.null(ref.n)){
      ref.n <- pwr::pwr.anova.test(f=efxize,power=power,sig.level=alpha,k=groups)$n
      message("ref.n is observations per group")
    }
    
    # Numerator
    num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    
  }else if(test == "chisq"){ # Chi-square
    
    # Needs degrees of freedom
    if(is.null(df)){
      stop("Chi-square is selected. Degrees of freedom must be set")
    }
    
    # Set effect size
    efxize <- switch(
      efxize,
      "small" = 0.10,
      "medium" = 0.30,
      "large" = 0.50
    )
    
    # Determine reference sample size
    if(is.null(ref.n)){
      ref.n <- pwr::pwr.chisq.test(w=efxize,df=df,power=power,sig.level=alpha)$N
    }
    # Numerator
    num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    
  }else if(test == "cor"){ # Correlation
    
    # Set effect size
    efxize <- switch(
      efxize,
      "small" = 0.10,
      "medium" = 0.30,
      "large" = 0.50
    )
    
    # Determine reference sample size
    if(is.null(ref.n)){
      ref.n <- pwr::pwr.r.test(r=efxize,power=power,sig.level=alpha)$n
    }
    
    # Numerator
    num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    
  }else if(any(c("one.sample", "two.sample", "paired") %in% test)){# t-test
    
    # Set effect size
    efxize <- switch(
      efxize,
      "small" = 0.20,
      "medium" = 0.50,
      "large" = 0.80
    )
    
    # Determine reference sample size
    if(is.null(ref.n)){
      ref.n <- pwr::pwr.t.test(d=efxize,power=power,sig.level=alpha,type=test)$n
    }
    
    # Numerator
    num <- sqrt(ref.n*(log(ref.n)+qchisq((1-alpha),1)))
    
  }else{stop("test does not exist")}
  
  # Denominator
  denom <- (sqrt(n*(log(n)+qchisq((1-alpha),1))))
  
  # Adjusted alpha calculation
  adj.a <- alpha*num/denom
  
  # Critical values
  if(test == "anova"){
    
    critical.f <- function (groups, n, a)
    {
      df1 <- groups - 1
      df2 <- n - groups
      cvf <- qf(a, df1, df2, lower.tail = FALSE)
      return(cvf)
    }
    
    cv <- critical.f(groups, n, adj.a)
    
  }else if(test == "chisq"){
    
    critical.chi <- function (df, a)
    {
      cvchi <- qchisq(a, df, lower.tail = FALSE)
      return(cvchi)
    }
    
    cv <- critical.chi(df, adj.a)
    
  }else if(test == "cor"){
    
    critical.r <- function (n, a)
    {
      df <- n - 2
      critical.t <- qt( a/2, df, lower.tail = FALSE )
      cvr <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
      return(cvr)
    }
    
    cv <- critical.r(n, adj.a)
    
  }else if(any(c("one.sample", "two.sample", "paired") %in% test)){
    
    critical.t <- function (n, a)
    {
      df <- n - 2
      cvt <- qt( a/2, df, lower.tail = FALSE )
      return(cvt)
    }
    
    cv <- critical.t(n, adj.a)
    
  }
  
  # Output
  output <- list(
    adapt.a = adj.a, crit.value = cv,
    orig.a = alpha, ref.n = ref.n,
    exp.n = n, power = power,
    efxize = efxize
  )
  # Check for ANOVA or Chi-square
  if(test == "anova"){
    output$groups <- groups
    output$df <- c((groups - 1), (n - groups))
    
  }else if(test=="chisq"){
    output$df <- df
  }
  # Add test
  output$test <- test
  
  return(output)
}

#%%%%%%%%%%%%%%%%%%%%%%%
# SYSTEM FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# OS and System Check ----
# Updated 08.09.2020
system.check <- function (...)
{
  OS <- unname(tolower(Sys.info()["sysname"]))
  
  RSTUDIO <- swiftelse(Sys.getenv("RSTUDIO") == "1", TRUE, FALSE)
  
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
        number <- 208
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
# Updated 26.07.2023
styletext <- function(
    text, defaults = c(
      "bold", "italics", "highlight",
      "underline", "strikethrough"
    )
)
{

  if(system.check()$TEXT)
  {
    
    if(missing(defaults)){
      number <- 0
    }else{
      
      number <- switch(
        defaults,
        bold = 1, italics = 3,
        underline = 4, highlight = 7,
        strikethrough = 9
      )
      
    }
    
    return(paste("\033[", number, ";m", text, "\033[0m", sep = ""))
    
  }else{
    return(text)
  }
  
}

#' @noRd
# Symbols ----
# Updated 26.07.2024
textsymbol <- function(
    symbol = c(
      "alpha", "beta", "chi", "delta",
      "eta", "gamma", "lambda", "omega",
      "phi", "pi", "rho", "sigma", "tau",
      "theta", "square root", "infinity",
      "check mark", "x", "bullet"
    )
)
{
  # Return code
  return(
    switch(
      symbol,
      alpha = "\u03B1", beta = "\u03B2", chi = "\u03C7",
      delta = "\u03B4", eta = "\u03B7", gamma = "\u03B3",
      lambda = "\u03BB,", omega = "\u03C9", phi = "\u03C6",
      pi = "\u03C0", rho = "\u03C1", sigma = "\u03C3", tau = "\u03C4", 
      theta = "\u03B8", "square root" = "\u221A", infinity = "\u221E", 
      "check mark" = "\u2713", x = "\u2717", bullet = "\u2022"
    )
  )
  
}

#' @noRd
# Title Case ----
# Not truly title case -- just capitializes
# the first letter of each word
# Uses `tools::toTitleCase` for actual title case
# Updated 25.06.2023
totitle <- function(string)
{
  
  # Split by spaces
  words <- unlist(strsplit(string, split = " "))
  
  # Set first letters to uppercase
  titleCased <- cvapply(words, function(x){
    
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

#' @noRd
# General function to check for packages ----
# Updated 13.07.2023
check_package <- function(packages)
{

  # Performs what original `installed.packages()` does
  # but without additional fluff
  installed <- packages %in% ulapply(.libPaths(), list.files)
  
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
          "`install.packages(c(", packages, "))`",
          "\n\nOnce installed, re-run this function (you may need to restart R/RStudio)."
        ), call. = FALSE
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
          "`install.packages(", packages, ")`",
          "\n\nOnce installed, re-run this function (you may need to restart R/RStudio)."
        ), call. = FALSE
      )
      
    }
    
  }
  
}
