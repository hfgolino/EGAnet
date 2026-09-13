#%%%%%%%%%%%%%%%%%%%%%%%%%#
#### {EGAnet} Plotting ####
#%%%%%%%%%%%%%%%%%%%%%%%%%#

# Internal functions that set up, validate, and render {EGAnet}'s
# network plots (built on top of {GGally}'s `ggnet2` and, for
# composite/multi-panel plots, {ggpubr}'s `ggarrange`).
#
# Moved out of `helpers.R` into their own file since they form a
# large, self-contained section of the package.

#%%%%%%%%%%%%%%%%%%%%%%
# GENERIC HELPERS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Create a memoized, zero-argument accessor for `compute()` ----
# Used to cache values (e.g., a package function's `formals()`, or a
# constant used in a hot path) that never change within a session but
# would otherwise be recomputed on every single plot call
# Updated 13.09.2026
memoize_once <- function(compute)
{

  # Cache (populated on first call)
  cached <- NULL

  # Return the memoized accessor
  return(
    function()
    {

      # Check for cache
      if(is.null(cached)){
        cached <<- compute()
      }

      # Return cache
      return(cached)

    }
  )

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

  # Get default `ggnet2` arguments (cached -- see `ggnet2_formals`)
  default_args <- ggnet2_formals()

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
      default_args$edge.color <- c("grey75", "grey25")
      default_args$edge.lty <- c("solid", "longdash")
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
# Cached base [0, 1] scaling sequence (with names) used by `rescale_edges` ----
# Independent of `edge_size`, so it's built once instead of on every plot
# Updated 13.09.2026
rescale_edges_sequence <- memoize_once(function(){
  scale_sequence <- seq.int(0, 1, 0.0001)
  names(scale_sequence) <- scale_sequence
  return(scale_sequence)
})

#' @noRd
# Re-scale edges ----
# Updated 13.09.2026
rescale_edges <- function(network, edge_size)
{

  # Set edge scaling (default `edge.size = 8`)
  edge_scaling <- rescale_edges_sequence() * edge_size

  # Return scaled edges
  return(unname(edge_scaling[as.character(abs(network))]))

}

#' @noRd
# Readable names ----
# Updated 13.02.2026
readable_names <- function(node_names)
{

  # Check for nodes without names
  if(all(node_names == "")){
    return(node_names)
  }

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
# Updated 13.02.2026
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

    # Cheat the network
    network[] <- 0.0001
    diag(network) <- 0

  }

  # Obtain number of communities
  communities <- unique_length(wc)

  # Whether the network has no communities (used a few times below)
  no_communities <- all(is.na(wc))

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
  if(no_communities){
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

  ## Node label size: coerce a bare `NA` (logical) to `NA_real_`
  ## so it passes as "numeric" in `GGally_errors` below
  if(length(plot_ARGS$label.size) == 1 && is.na(plot_ARGS$label.size)){
    plot_ARGS$label.size <- NA_real_
  }

  ## `label.size` of `0` or `NA` means "no node labels"
  show_labels <- !all(is.na(plot_ARGS$label.size) | plot_ARGS$label.size == 0)

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
  if(no_communities){ # Plain network (without communities)
    border_color <- rep("grey", dimensions[2])
  }else if(
    length(color.palette) == 1 &&
    color.palette %in% gray_options
  ){ # Gray scale network
    border_color <- swiftelse(palette == "white", "white", "grey")
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
    )

  # Only add node labels back on top if `label.size` isn't `0` or `NA`
  if(show_labels){
    second_layer <- second_layer +
      ggplot2::geom_text(
        ggplot2::aes(label = node_names),
        color = plot_ARGS$label.color,
        size = plot_ARGS$label.size
      )
  }

  second_layer <- second_layer +
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
  }else if(no_communities){ # no legend (network with no communities plot)
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
# Updated 21.11.2025
compare_plots <- function(comparison_network, comparison_wc, plot_ARGS, ...)
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

  # Check if any arguments in `ellipse` match with `plot_ARGS`
  plot_ARGS <- overwrite_arguments(plot_ARGS, list(...))

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ARGUMENT VALIDATION AND PASSING FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Cached `formals()` for external plotting functions ----
# Shared by the "allowed names" functions below and by
# `GGally_args`/`ggarrange_args`, so each is only ever computed once
# Updated 13.09.2026
ggnet2_formals <- memoize_once(function() formals(silent_load(GGally::ggnet2)))
ggarrange_formals <- memoize_once(function() formals(silent_load(ggpubr::ggarrange)))
theme_formals <- memoize_once(function() formals(silent_load(ggplot2::theme)))

#' @noRd
# Allowed `ggnet2` argument names ----
# Real `ggnet2` formals plus {EGAnet}'s own friendly
# aliases, extras, legacy shims, and internal control flag
# Updated 13.09.2026
ggnet2_allowed_names <- memoize_once(function(){
  c(
    # Real `ggnet2` formal names
    setdiff(names(ggnet2_formals()), "..."),
    # {EGAnet} friendly aliases (see `GGally_args`)
    "layout", "alpha", "color", "shape", "vsize",
    # {EGAnet} extras consumed directly from `ellipse`
    # in `basic_plot_setup` (never forwarded to `ggnet2`)
    "title", "legend.title", "legend.names",
    # Legacy list-shims flattened by `legacy_EGA_args`
    "model.args", "algorithm.args", "plot.args",
    # Documented flag (see `EGAnet-plot`'s "Argument Passing" section) that
    # changes the return value to `list(network_plot=, ARGS=)` instead of
    # just the plot (never forwarded to `ggnet2` itself)
    "arguments"
  )
})

#' @noRd
# Allowed `ggpubr::ggarrange` argument names ----
# Updated 13.09.2026
ggarrange_allowed_names <- memoize_once(function(){
  setdiff(names(ggarrange_formals()), "...")
})

#' @noRd
# Allowed `ggplot2::theme` argument names ----
# `theme()`'s element names are real formals (not just `...`), so
# these are introspected the same way as `ggnet2`/`ggarrange` --
# no hand-maintained list, stays in sync with the installed
# {ggplot2} version automatically. Used only where a single
# `...` legitimately feeds `ggnet2`/`ggarrange` AND `theme()`
# (e.g., `plot.bootEGA`)
# Updated 13.09.2026
theme_allowed_names <- memoize_once(function(){
  setdiff(names(theme_formals()), "...")
})

#' @noRd
# Keep only the `ggnet2`-relevant names in an `ellipse` ----
# Used where a single `...` legitimately feeds both an internally
# dispatched `ggnet2`-only method (e.g., `plot.EGA`) and a separate
# `ggarrange`/`theme` destination -- only the `ggnet2` subset should
# be forwarded to the former, since it validates strictly on its own
# Updated 13.09.2026
filter_ggnet2_ellipse <- function(ellipse)
{
  ellipse[names(ellipse) %in% ggnet2_allowed_names()]
}

#' @noRd
# Error for unrecognized plotting argument names ----
# Catches typos that would otherwise be silently dropped
# by `overwrite_arguments`/`obtain_arguments`
# Updated 13.09.2026
argument_name_error <- function(ellipse, allowed_names, function_name)
{

  # Ignore unnamed/positional arguments (e.g., `*EGA` objects
  # passed into `compare.EGA.plots`)
  supplied_names <- names(ellipse)
  supplied_names <- supplied_names[nzchar(supplied_names)]

  # Determine which supplied names are not recognized
  unrecognized <- supplied_names[!supplied_names %in% allowed_names]

  # Throw error listing every unrecognized name at once
  if(length(unrecognized) != 0){
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Unrecognized argument", swiftelse(length(unrecognized) > 1, "s", ""),
        ": ", paste0("'", unrecognized, "'", collapse = ", "),
        ". Check for typos -- see `?GGally::ggnet2`",
        swiftelse(
          any(allowed_names %in% ggarrange_allowed_names()),
          ", `?ggpubr::ggarrange`,", ","
        ), " and `?EGAnet::EGAnet-plot` for valid argument names.",
        "\n\n For more details on how to fix this error, see:\n",
        "https://r-ega.net/articles/errors.html#argument-name-error"
      ),
      call = function_name
    )
  }

}

#' @noRd
# Defaults for `ggpubr::ggarrange` plotting ----
# Mirrors `GGally_args`, for composite (multi-panel) plot methods
# Updated 13.09.2026
ggarrange_args <- function(ellipse, site_defaults = list(), plotlist = NULL)
{

  # Get default `ggarrange` arguments
  default_args <- ggarrange_formals()

  # Replace `ggarrange` arguments with this call site's own defaults
  # (e.g., `ncol = 2, nrow = 1, labels = , legend = "bottom"`)
  default_args <- overwrite_arguments(default_args, site_defaults)

  # Replace `ggarrange` arguments with the user's arguments
  # (only real `ggarrange` names survive -- already name-validated upstream)
  default_args <- overwrite_arguments(default_args, ellipse)

  # Remove the ellipsis
  default_args <- default_args[names(default_args) != "..."]

  # Remove any remaining un-evaluated defaults (e.g., `hjust`, `font.label`,
  # `align`) so `do.call` doesn't choke on a raw language object --
  # `ggarrange` will apply its own default for anything omitted here
  # (same convention as `obtain_arguments`)
  default_args <- default_args[!lvapply(default_args, is, "call")]

  # Set plot list last (never user-overridable), if supplied
  if(!is.null(plotlist)){
    default_args$plotlist <- plotlist
  }

  # Return arguments (ready for `do.call(ggpubr::ggarrange, ...)`)
  return(default_args)

}
