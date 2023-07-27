#' @title Visually Compare Two or More \code{\link{EGAnet}} plots
#' 
#' @description Organizes EGA plots for comparison. Ensures that 
#' nodes are placed in the same layout to maximize comparison
#'
#' @param ... Handles multiple arguments:
#' 
#' \itemize{
#' 
#' \item{\code{*EGA} objects}
#' {can be dropped in without any argument
#' designation. The function will search across input to find
#' necessary \code{\link{EGAnet}} objects}
#' 
#' \item{\code{\link[GGally]{ggnet2}} arguments}
#' {can be passed along to \code{\link[GGally]{ggnet2}}}
#' 
#' \item{\code{\link[sna]{gplot.layout}}}
#' {can be specified using \code{mode = } or
#' \code{layout = } using the name of the layout
#' (e.g., \code{mode = "circle"} will produce the 
#' circle layout from \link[sna]{gplot.layout}).
#' By default, the layout is the same as \code{\link{qgraph}}}
#'
#' }
#' 
#' @param input.list List.
#' Bypasses \code{...} argument in favor of using a list
#' as an input
#' 
#' @param base.plot Numeric (length = 1).
#' Plot to be used as the base for the configuration of the networks.
#' Uses the number of the order in which the plots are input.
#' Defaults to \code{1} or the first plot
#' 
#' @param labels Character (same length as input).
#' Labels for each \code{\link{EGAnet}} object
#' 
#' @param rows Numeric (length = 1).
#' Number of rows to spread plots across
#' 
#' @param columns Numeric (length = 1).
#' Number of columns to spread plots down
#'
#' @return Visual comparison of \code{\link{EGAnet}} objects
#' 
#' @examples
#' # Obtain WMT-2 data
#' wmt <- wmt2[,7:24]
#' 
#' # Draw random samples of 300 cases
#' sample1 <- wmt[sample(1:nrow(wmt), 300),]
#' sample2 <- wmt[sample(1:nrow(wmt), 300),]
#' 
#' # Estimate EGAs
#' ega1 <- EGA(sample1)
#' ega2 <- EGA(sample2)
#' 
#' # Compare EGAs via plot
#' compare.EGA.plots(
#'   ega1, ega2,
#'   base.plot = 1, # use "ega1" as base for comparison
#'   labels = c("Sample 1", "Sample 2"),
#'   rows = 1, columns = 2
#' )
#' 
#' # Change layout to circle plots
#' compare.EGA.plots(
#'   ega1, ega2,
#'   labels = c("Sample 1", "Sample 2"),
#'   mode = "circle"
#' )
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#
# Compare EGA plots ----
# Updated 27.07.2023
compare.EGA.plots <- function(
  ..., input.list = NULL, base.plot = 1,
  labels = NULL, rows = NULL, columns = NULL
)
{
  
  # Start with ellipse objects/arguments
  ellipse <- list(...)
  
  # Determine whether it's necessary to search for
  # plots in the ellipse
  if(is.null(input.list)){
    
    # Get classes for all input
    classes <- lapply(ellipse, class)
    
    # Determine which inputs are `EGA`
    ega_object <- lvapply(classes, function(x){grepl("EGA", x)})
    
    # Check for no `EGA` objects
    if(!any(ega_object)){
      stop(
        "The 'input.list' was `NULL` and there were no `EGA` class objects detected.",
        call. = FALSE
      )
    }
    
    # Extract `EGA` objects from list
    input.list <- ellipse[ega_object]
    ellipse <- ellipse[!ega_object]
    
  }
  
  # With the input list, there could be different types of
  # `EGA` objects... let's figure that out
  ega_classes <- cvapply(input.list, class)
  
  # Separate `dynEGA` from rest of `EGA` functions
  dynega_classes <- grepl("dynEGA", ega_classes)
  
  # Separate input list
  ega_input <- input.list[!dynega_classes]
  dynega_input <- input.list[dynega_classes]
  
  # For `EGA` input, extract `EGA` objects
  if(length(ega_input) != 0){
    ega_list <- lapply(ega_input, get_EGA_object)
  }else{ # Set up empty list
    ega_list <- list()
  }
  
  # For `dynEGA` input, extract `EGA` objects
  if(length(dynega_input) != 0){
    
    # Extract `dynEGA` objects
    dynega_objects <- ulapply(
      dynega_input, get_EGA_object, recursive = FALSE
    )
    
    # Check for `population` objects
    population_objects <- names(dynega_objects) == "population"
    
    # Separate `population` objects from the rest
    dynega_population <- dynega_objects[population_objects]
    dynega_other <- unlist( # Extract `EGA` objects from rest
      dynega_objects[!population_objects], recursive = FALSE
    )
      
    # Put all objects together
    dynega_list <- c(
      dynega_population, dynega_other
    )
    
  }else{ # Set up empty list
    dynega_list <- list()
  }
  
  # Compile all `EGA` objects
  input.list <- c(ega_list, dynega_list)
  
  # Organize input list with base plot
  if(is.character(base.plot)){
    base.plot <- which(names(input.list) == base.plot)
  }
  
  # Re-assemble input list
  input.list <- c(
    input.list[base.plot], input.list[-base.plot]
  )
  
  # Assign "labels"
  if(!is.null(labels)){
    names(input.list) <- labels
  }
  
  # Set up as `dynEGA.Individual` object to
  # leverage the code already implemented there
  class(input.list) <- "dynEGA.Individual"
  
  # For `ellipse`, get legacy arguments
  ellipse <- legacy_EGA_args(ellipse)
  
  # Get length of input list
  input_length <- length(input.list)
  
  # Handle rows and columns
  if(is.null(rows) & is.null(columns)){
    
    # Set rows first and then columns
    rows <- floor(sqrt(input_length))
    columns <- ceiling(input_length / rows)
    # use `ceiling` since any non-zero 
    # remainder means an extra column is necessary
    
  }else if(is.null(rows)){ # Set rows based on columns
    rows <- ceiling(input_length / columns)
  }else if(is.null(columns)){ # Set columns based on rows
    columns <- ceiling(input_length / rows)
  }
  
  # Add labels, rows, and columns
  ellipse[c("nrow", "ncol")] <- list(rows, columns)
  
  # Plot using all "IDs"
  return(
    do.call(
      what = plot,
      args = c(
        list(
          x = input.list,
          id = seq_len(input_length)
        ),
        ellipse
      )
    )
  )

}

