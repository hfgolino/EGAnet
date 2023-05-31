#%%%%%%%%%%%%%%%%%%%%%%%
# SYSTEM FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
#'
# General function to perform
# system-specific parallelization on lists
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

#' Error report
#'
#' @description Gives necessary information for user reporting error
#'
#' @param result Character.
#' The error from the result
#'
#' @param SUB_FUN Character.
#' Sub-routine the error occurred in
#'
#' @param FUN Character.
#' Main function the error occurred in
#'
#' @return Error and message to send to GitHub
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
#' @importFrom utils packageVersion
#'
# Error Report
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

#' System check for OS and RSTUDIO
#'
#' @description Checks for whether text options are available
#'
#' @param ... Additional arguments
#'
#' @return \code{TRUE} if text options are available and \code{FALSE} if not
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# System Check
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

#' Colorfies Text
#'
#' Makes text a wide range of colors (8-bit color codes)
#'
#' @param text Character.
#' Text to color
#'
#' @return Colorfied text
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Color text
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

#' Stylizes Text
#'
#' Makes text bold, italics, underlined, and strikethrough
#'
#' @param text Character.
#' Text to stylized
#'
#' @return Sytlized text
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# Style text
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

#' Text Symbols
#'
#' Makes text symbols (star, checkmark, square root)
#'
#' @param symbol Character.
#'
#' @return Outputs symbol
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# Symbols
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
