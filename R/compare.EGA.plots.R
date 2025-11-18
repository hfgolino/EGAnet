#' @title Visually Compare Two or More \code{EGAnet} plots
#'
#' @description Organizes EGA plots for comparison. Ensures that
#' nodes are placed in the same layout to maximize comparison
#'
#' @param ... Handles multiple arguments:
#'
#' \itemize{
#'
#' \item \code{*EGA} objects --- can be dropped in without any argument
#' designation. The function will search across input to find
#' necessary \code{EGAnet} objects
#'
#' \item \code{\link[GGally]{ggnet2}} arguments --- can be passed along to \code{\link[GGally]{ggnet2}}
#'
#' \item \code{\link[sna]{gplot.layout}} --- can be specified using \code{mode = } or
#' \code{layout = } using the name of the layout
#' (e.g., \code{mode = "circle"} will produce the
#' circle layout from \link[sna]{gplot.layout}).
#' By default, the layout is the same as \code{qgraph}
#'
#' }
#'
#' @param input.list List.
#' Bypasses \code{...} argument in favor of using a list
#' as an input
#'
#' @param base Numeric (length = 1).
#' Plot to be used as the base for the configuration of the networks.
#' Uses the number of the order in which the plots are input.
#' Defaults to \code{1} or the first plot
#'
#' @param same.layout Boolean (length = 1).
#' Whether the nodes should be in the same position for all networks.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} for their individual layouts
#'
#' @param labels Character (same length as input).
#' Labels for each \code{EGAnet} object
#'
#' @param rows Numeric (length = 1).
#' Number of rows to spread plots across
#'
#' @param columns Numeric (length = 1).
#' Number of columns to spread plots down
#'
#' @param plot.all Boolean (length = 1).
#' Whether plot should be produced or just output.
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to avoid plotting (but still obtain
#' plot objects)
#'
#' @return Visual comparison of \code{EGAnet} objects
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
#' \donttest{
#' # Compare EGAs via plot
#' compare.EGA.plots(
#'   ega1, ega2,
#'   base = 1, # use "ega1" as base for comparison
#'   labels = c("Sample 1", "Sample 2"),
#'   rows = 1, columns = 2
#' )
#'
#' # Change layout to circle plots
#' compare.EGA.plots(
#'   ega1, ega2,
#'   labels = c("Sample 1", "Sample 2"),
#'   mode = "circle"
#' )}
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{EGAnet}
#'
#' @export
#
# Compare EGA plots ----
# Updated 18.11.2025
compare.EGA.plots <- function(
  ..., input.list = NULL, base = 1, same.layout = TRUE,
  labels = NULL, rows = NULL, columns = NULL,
  plot.all = TRUE
)
{

  # Error check 'same.layout'
  length_error(same.layout, 1, "compare.EGA.plots")
  typeof_error(same.layout, "logical", "compare.EGA.plots")

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
      .handleSimpleError(
        h = stop,
        msg = "The 'input.list' was `NULL` and there were no `EGA` class objects detected.",
        call = "compare.EGA.plots"
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

    # Remove ".Ord*" labels
    dynega_list <- lapply(dynega_list, function(x){

      # New names
      new_names <- gsub(".Ord*.", "", dimnames(x$network)[[2]])

      # Change the network names
      dimnames(x$network) <- list(new_names, new_names)

      # Change the membership names
      names(x$wc) <- new_names

      # Return full result
      return(x)

    })

  }else{ # Set up empty list
    dynega_list <- list()
  }

  # Compile all `EGA` objects
  input.list <- c(ega_list, dynega_list)

  # Organize input list with base plot
  if(is.character(base)){
    base <- which(names(input.list) == base)
  }

  # Re-assemble input list
  input.list <- c(
    input.list[base], input.list[-base]
  )

  # Assign "labels"
  if(!is.null(labels)){
    names(input.list) <- labels
  }

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

  # Set up for individual plots

  # Remove base from individuals
  base_object <- input.list[[base]]

  # Extract other groups
  other_objects <- input.list[-base]

  # Get sequence length of other objects
  sequence_length <- seq_len(length(other_objects))

  # Get base plot
  base_plot <- silent_load(
    do.call(
      what = basic_plot_setup,
      args = c(
        list(
          network = base_object$network,
          wc = base_object$wc,
          arguments = TRUE
        ), ellipse
      )
    )
  )

  # Check for same layout
  if(same.layout){

    # Check if any arguments in `ellipse`
    # match with `base_plot`
    base_plot$ARGS <- overwrite_arguments(base_plot$ARGS, ellipse)

    # Set removal arguments
    removal_ARGS <- c(
      "node.color", "edge.alpha",
      "edge.color", "edge.lty", "edge.size"
    )

    # Check for if any arguments in still need
    # to be removed from `base_plot`
    removal_ARGS <- removal_ARGS[!removal_ARGS %in% names(ellipse)]

    # Check for any remaining arguments to remove
    if(length(removal_ARGS) != 0){

      # Remove some arguments from `base_plot`
      base_plot$ARGS <- base_plot$ARGS[
        !names(base_plot$ARGS) %in% removal_ARGS
      ]

    }

    # Set up comparison plots
    comparison_plots <- lapply(
      sequence_length, function(i){
        compare_plots(
          comparison_network = other_objects[[i]]$network,
          comparison_wc = other_objects[[i]]$wc,
          plot_ARGS = base_plot$ARGS
        )
      }
    )

  }else{

    # Create comparison plots without same layout
    comparison_plots <- lapply(
      sequence_length, function(i){
        silent_load(
          do.call(
            what = basic_plot_setup,
            args = c(
              list(
                network = other_objects[[i]]$network,
                wc = other_objects[[i]]$wc
              ), ellipse
            )
          )
        )
      }
    )


  }

  # Set up plot list
  plotlist <- c(list(base_plot$network_plot), comparison_plots)

  # `ggarrange` does not like non-plot arguments for its ellipse
  ellipse <- ellipse[
    names(ellipse) %in% names(formals(ggpubr::ggarrange))
  ]

  # Store plots all-in-one
  all_in_one <- do.call(
    what = ggpubr::ggarrange,
    args = c(
      list(
        plotlist = plotlist,
        labels = labels,
        legend = "bottom"
      ), ellipse
    )
  )

  # Should the plot be produced?
  if(isTRUE(plot.all)){
    silent_plot(all_in_one)
  }

  # Check for labels
  if(is.null(labels)){

    # Plot numbers
    plot_numbers <- seq_along(plotlist)

    # Get non-base plot labels
    not_base_labels <- plot_numbers[!plot_numbers %in% base]

    # Set labels
    names(plotlist) <- c(base, not_base_labels)

  }

  # Return plot lists
  return(
    list(
      all = all_in_one,
      individual = plotlist
    )
  )

}

