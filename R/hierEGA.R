#' Hierarchical \code{\link[EGAnet]{EGA}}
#'
#' Estimates EGA using the lower-order solution of \code{\link[igraph]{cluster_louvain}}
#' to identify the lower-order dimensions and then uses factor or network loadings to
#' estimate factor or network scores, which are used to estimate the higher-order dimensions
#'
#' @param data Matrix or data frame.
#' Variables (down columns) only.
#' Does not accept correlation matrices
#'
#' @param scores Character.
#' How should scores for the higher-order structure be estimated?
#' Defaults to \code{"network"} for network scores computed using
#' the \code{\link[EGAnet]{net.scores}} function.
#' Set to \code{"factor"} for factor scores computed using
#' \code{\link[psych]{fa}}. Factors are assumed to be correlated
#' using the \code{"oblimin"} rotation. \emph{NOTE}: Factor scores
#' use the number of communities from \code{\link[EGAnet]{EGA}}.
#' Estimated factor may not align with these communities. The plots
#' using factor scores with have higher order factors that may not
#' completely map onto the lower order communities. Look at the
#' \code{$hierarchical$higher_order$lower_loadings} to determine
#' the composition of the lower order factors.
#'
#' By default, both factor and network scores are computed and stored
#' in the output. The selected option only appears in the main output (\code{$hierarchical})
#'
#' @param rotation Character.
#' A rotation to use, like factor loadings, to obtain
#' a simple structure. For a list of rotations,
#' see \link{GPArotation}
#'
#' @param consensus.iter Numeric.
#' Number of iterations to perform in consensus clustering
#' (see Lancichinetti & Fortunato, 2012).
#' Defaults to \code{1000}
#'
#' @param consensus.method Character.
#' What consensus clustering method should be used?
#' Defaults to \code{"highest_modularity"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{highest_modularity}}}
#' {Uses the community solution that achieves the highest modularity
#' across iterations}
#'
#' \item{\strong{\code{most_common}}}
#' {Uses the community solution that is found the most
#' across iterations}
#'
#' \item{\strong{\code{iterative}}}
#' {Identifies the most common community solutions across iterations
#' and determines how often nodes appear in the same community together.
#' A threshold of 0.30 is used to set low proportions to zero.
#' This process repeats iteratively until all nodes have a proportion of
#' 1 in the community solution.
#' }
#'
#' \item{\code{lowest_tefi}}
#' {Uses the community solution that achieves the lowest \code{\link[EGAnet]{tefi}}
#' across iterations}
#' 
#' \item{\code{most_common_tefi}}
#' {Uses the most common number of communities detected across the number
#' of iterations. After, if there is more than one solution for that number
#' of communities, then the solution with the lowest \code{\link[EGAnet]{tefi}
#' is used}}
#' 
#' }
#'
#' By default, all \code{consensus.method} options are computed and
#' stored in the output. The selected method will be used to
#' plot and appear in the main output (\code{$hierarchical})
#'
#' @param uni.method Character.
#' What unidimensionality method should be used? 
#' Defaults to \code{"LE"}.
#' Current options are:
#' 
#' \itemize{
#'
#' \item{\strong{\code{expand}}}
#' {Expands the correlation matrix with four variables correlated .50.
#' If number of dimension returns 2 or less in check, then the data 
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This is the method used in the Golino et al. (2020)
#' \emph{Psychological Methods} simulation.}
#'
#' \item{\strong{\code{LE}}}
#' {Applies the Leading Eigenvalue algorithm (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvalue solution is used; otherwise, regular EGA
#' is used. This is the final method used in the Christensen, Garrido,
#' and Golino (2021) simulation.}
#' 
#' \item{\strong{\code{louvain}}}
#' {Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix using a resolution parameter = 0.95.
#' If the number of dimensions is 1, then the Louvain solution is used; otherwise,
#' regular EGA is used. This method was validated in the Christensen (2022) simulation.}
#' 
#' }
#'
#' @param corr Type of correlation matrix to compute. The default uses \code{\link[qgraph]{cor_auto}}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{cor_auto}}}
#' {Computes the correlation matrix using the \code{\link[qgraph]{cor_auto}} function from
#' \code{\link[qgraph]{qgraph}}}.
#'
#' \item{\strong{\code{pearson}}}
#' {Computes Pearson's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#'
#' \item{\strong{\code{spearman}}}
#' {Computes Spearman's correlation coefficient using the pairwise complete observations via
#' the \code{\link[stats]{cor}}} function.
#' }
#'
#' @param model Character.
#' A string indicating the method to use.
#' Defaults to \code{"glasso"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{glasso}}}
#' {Estimates the Gaussian graphical model using graphical LASSO with
#' extended Bayesian information criterion to select optimal regularization parameter}
#'
#' \item{\strong{\code{TMFG}}}
#' {Estimates a Triangulated Maximally Filtered Graph}
#'
#' }
#'
#' @param model.args List.
#' A list of additional arguments for \code{\link[EGAnet]{EBICglasso.qgraph}}
#' or \code{TMFG}
#'
#' @param algorithm A string indicating the algorithm to use or a function from \code{\link{igraph}}
#' Defaults to \code{"walktrap"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
#' 
#' \item{\strong{\code{leiden}}}
#' {Computes the Leiden algorithm using \code{\link[igraph]{cluster_leiden}}.
#' Defaults to \code{objective_function = "modularity"}}
#'
#' \item{\strong{\code{louvain}}}
#' {Computes the Louvain algorithm using \code{\link[igraph]{cluster_louvain}}}
#'
#' }
#'
#' @param algorithm.args List.
#' A list of additional arguments for \code{\link[igraph]{cluster_walktrap}}, \code{\link[igraph]{cluster_louvain}},
#' or some other community detection algorithm function (see examples)
#'
#' @param lower.louvain Boolean.
#' Should lower Louvain solution be used at the higher order level?
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to use the lower Louvain solution 
#'
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#' @param plot.args List.
#' A list of additional arguments for the network plot. See \code{\link[GGally]{ggnet2}} for
#' full list of arguments:
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#' \item{\strong{\code{label.size}}}
#' {Size of the labels. Defaults to 5.}
#'
#' \item{\strong{\code{alpha}}}
#' {The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.7.}
#'
#' \item{\strong{\code{edge.alpha}}}
#' {The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.4.}
#'
#'  \item{\strong{\code{legend.names}}}
#' {A vector with names for each dimension}
#'
#' \item{\strong{\code{color.palette}}}
#' {The color palette for the nodes. For custom colors,
#' enter HEX codes for each dimension in a vector.
#' See \code{\link[EGAnet]{color_palette_EGA}} for
#' more details and examples}
#'
#' }
#'
#' @return Returns a list of lists containing: \cr
#' \cr
#' \strong{Main Results} \cr
#'
#' \item{hierarhical}{The main results list containing:
#'
#' \itemize{
#'
#' \item{\code{lower_order}}
#' {Lower order \code{\link[EGAnet]{EGA}} results for the selected methods}
#'
#' \item{\code{higher_order}}
#' {Higher order \code{\link[EGAnet]{EGA} results for the selected methods}}
#'
#' If \code{plot.EGA = TRUE}, then:
#'
#' \item{\code{lower_plot}}
#' {Plot of the lower order results}
#'
#' \item{\code{higher_plot}}
#' {Plot of the higher order results}
#'
#' \item{\code{hier_plot}}
#' {Plot of the lower and higher order results together, side-by-side}
#'
#' }
#'
#' }
#'
#' \strong{Secondary Results} \cr
#'
#' \item{lower_ega}{A list containing the lower order \code{\link[EGAnet]{EGA}}
#' results. The \code{$wc} does not contain valid results. Do not use its output.
#' }
#'
#' \item{lower_wc}{A list containing consensus clustering results:
#'
#' \itemize{
#'
#' \item{\code{highest_modularity}}
#' {Community memberships based on the highest modularity across the
#' \code{\link[igraph]{cluster_louvain}} applications}
#'
#' \item{\code{most_common}}
#' {Community memberships based on the most commonly found memberships across the
#' \code{\link[igraph]{cluster_louvain}} applications}
#'
#' \item{\code{iterative}}
#' {Community memberships based on consensus clustering described by
#' Lancichinetti & Fortunato (2012)}
#'
#' \item{\code{lowest_tefi}}
#' {Community memberships based on the lowest \code{\link[EGAnet]{tefi}}
#' across the \code{\link[igraph]{cluster_louvain}} applications}
#'
#' \item{\code{summary_table}}
#' {A data frame summarizing the unique community solutions across the iterations. Down the
#' columns indicate: number of dimensions (\code{N_Dimensions}),
#' proportion of times each community solution was identified (\code{Proportion}),
#' modularity of each community solution (\code{Modularity}),
#' total entropy fit index of each community solution (\code{\link[EGAnet]{tefi}}),
#' and the memberships for each item. Across the rows indicate each
#' unique community solution
#' }
#'
#' }
#'
#' }
#'
#' \item{factor_results}{A list containing higher order results based on factor scores.
#' A list for each \code{consensus.method} is provided with their \code{\link[EGAnet]{EGA}} results}
#'
#' \item{network_results}{A list containing higher order results based on network scores.
#' A list for each \code{consensus.method} is provided with their \code{\link[EGAnet]{EGA}} results}
#'
#' @references
#' Lancichinetti, A., & Fortunato, S. (2012).
#' Consensus clustering in complex networks.
#' \emph{Scientific Reports}, \emph{2}(1), 1-7.
#'
#' @author
#' Marcos Jimenez <marcosjnezhquez@gmailcom>,
#' Francisco J. Abad <fjose.abad@uam.es>,
#' Eduardo Garcia-Garzon <egarcia@ucjc.edu>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>, and
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu.do>
#'
#' @examples
#' # Obtain example data
#' data <- optimism
#' 
#' \dontrun{
#' # hierEGA example
#' opt.hier<- hierEGA(
#'   data = optimism,
#'   algorithm = "louvain"
#' )}
#'
#' @export
#'
# Hierarchical EGA ----
# Updated 31.07.2023
hierEGA <- function(
    data, 
    # `net.scores` arguments
    loading.method = c("BRM", "experimental"),
    rotation = NULL, scores = c("factor", "network"),
    loading.structure = c("simple", "full"),
    impute = c("mean", "median", "none"),
    # `EGA` arguments
    corr = c("auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),
    lower.algorithm = "louvain",
    higher.algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    plot.EGA = TRUE, verbose = FALSE,
    ...
)
{

  # Get ellipse arguments
  ellipse <- list(...)
  
  # Check for missing arguments (argument, default, function)
  ## `net.scores`
  loading.method <- set_default(loading.method, "brm", net.loads)
  scores <- set_default(scores, "network", hierEGA)
  loading.structure <- set_default(loading.structure, "simple", hierEGA)
  impute <- set_default(impute, "none", net.scores)
  ## `EGA`
  corr <- set_default(corr, "auto", c("auto", "cor_auto", "pearson", "spearman"))
  corr <- swiftelse(corr == "cor_auto", "auto", corr) # deprecate `cor_auto`
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  
  ## Determine if single algorithm is entered
  ## Overwrite 'lower.algorithm' and 'higher.algorithm'
  if("algorithm" %in% names(ellipse)){
    
    # Set lower and higher order algorithm to single algorithm
    lower.algorithm <- ellipse$algorithm
    higher.algorithm <- ellipse$algorithm
    
    # Remove 'algorithm' from ellipse to avoid conflicts
    ellipse <- ellipse[names(ellipse) != "algorithm"]
    
  }
  
  ## Handle 'lower.algorithm' and 'higher.algorithm' arguments (workaround)
  algorithm <- lower.algorithm
  lower.algorithm <- set_default(algorithm, "louvain", community.detection)
  algorithm <- higher.algorithm
  higher.algorithm <- set_default(algorithm, "walktrap", community.detection)
  
  # Get EGA
  lower_order_result <- do.call(
    what = EGA,
    args = c(
      # Standard `EGA` arguments
      list(
        data = data, corr = corr, na.data = na.data,
        model = model, algorithm = lower.algorithm, uni.method = uni.method,
        plot.EGA = FALSE, verbose = verbose, order = "lower"
      ),
      # Send ellipse arguments
      ellipse
    )
  )

  # Ensure data has names (needed in `net.scores`)
  if(is.null(dimnames(data)[[2]])){
    dimnames(data)[[2]] <- dimnames(lower_order_result$network)[[2]]
  }

  # Check for scores
  if(scores == "factor"){
    
    # Send warning
    warning(
      paste0(
        "Factor dimensions may not align with `EGA` detected communities. ",
        "Interpret results with caution. ",
        "\n\nUse `$higher_order$lower_loadings` to guide interpretations"
      ),
      call. = FALSE
    )
    
    # Get arguments for EFA
    efa_ARGS <- obtain_arguments(efa_scores, list(...))
    
    # Set data, correlation matrix, and number of factors
    efa_ARGS[
      c("data", "correlation_matrix", "nfactors", "impute")
    ] <- list(
      data, lower_order_result$correlation, 
      lower_order_result$n.dim, impute
    )
    
    # Estimate EFA
    efa_output <- do.call(efa_scores, efa_ARGS)
    
    # Get scores
    score_estimates <- efa_output$scores
  
    # Get loading names
    loading_names <- dimnames(efa_output$loadings)
    
    # Obtain highest loadings for each variable
    lower_wc <- max.col(abs(efa_output$loadings), "first")
    
    # Get max `wc`
    max_wc <- max(lower_wc, na.rm = TRUE)
    
    # Format `wc` to be characters
    lower_wc <- format_integer(lower_wc, digits(max_wc) - 1)
    names(lower_wc) <- loading_names[[1]]
    
    # Set up dimension names to mirror network output
    dimnames(efa_output$loadings)[[2]] <- format_integer(
      as.numeric(gsub("MR", "", loading_names[[2]])),
      digits(max_wc) - 1
    )
    dimnames(score_estimates)[[2]] <- dimnames(efa_output$loadings)[[2]]
    
    # Put loadings in descending order
    lower_loadings <- descending_order(
      standardized = efa_output$loadings, wc = lower_wc,
      unique_communities = sort(dimnames(efa_output$loadings)[[2]]),
      node_names = loading_names[[1]]
    )
    
  }else if(scores == "network"){
    
    # Compute network scores
    network_output <- net.scores(
      data = data, A = lower_order_result,
      rotation = rotation, loading.method = loading.method,
      scoring.method = "network", 
      loading.structure = loading.structure,
      impute = impute,
      ...
    )
    
    # Score estimates
    if(is.null(rotation)){
      score_estimates <- network_output$scores$std.scores
      lower_loadings <- network_output$loadings$std
    }else{
      score_estimates <- network_output$scores$rot.scores
      lower_loadings <- network_output$loadings$rotated
    }
    
  }
  
  # Check for unidimensional lower order result
  # (rare but can happen)
  if(lower_order_result$n.dim == 1){
    
    # Set up results
    results <- list(
      lower_order = lower_order_result,
      higher_order = lower_order_result,
      parameters = list(
        lower_loadings = lower_loadings,
        lower_scores = score_estimates
      )
    )

  }else{
    
    # Set up results
    results <- list(
      lower_order = lower_order_result,
      higher_order = do.call(
        what = EGA,
        args = c(
          # Standard `EGA` arguments
          list(
            data = score_estimates, corr = corr, na.data = na.data,
            model = model, algorithm = higher.algorithm, uni.method = uni.method,
            plot.EGA = FALSE, verbose = verbose
          ),
          # Send ellipse arguments
          ellipse
        )
      ),
      parameters = list(
        lower_loadings = lower_loadings,
        lower_scores = score_estimates
      )
    )
    
  }
  
  # Perform parallel PCA to check for no general factors
  # sink <- capture.output(
  #   pca <-
  #     psych::fa.parallel(
  #       x = results$higher_order$correlation,
  #       fa = "pc",
  #       n.obs = dim(data)[1],
  #       plot = FALSE
  #     )
  # )
  # Look into alternative solutions
  
  # Add dimension variables like `EGA`
  results$dim.variables <- fast.data.frame(
    c(
      dimnames(data)[[2]], 
      as.vector(results$lower_order$wc),
      as.vector(
        single_revalue_memberships( # function in `bootEGA` internals
          results$lower_order$wc, results$higher_order$wc
        )
      )
    ),
    nrow = length(results$lower_order$wc), ncol = 3,
    colnames = c("items", "lower", "higher")
  )
  
  # Set "methods" attributes
  attr(results, "methods") <- list(
    loading.method = loading.method, 
    rotation = rotation,
    scores = scores
  )
  
  # Set class
  class(results) <- "hierEGA"
  
  # Add TEFI to the result
  results$TEFI <- tefi(results)$VN.Entropy.Fit
  
  # Check for plot
  if(lower_order_result$n.dim != 1 && isTRUE(plot.EGA)){
    
    # Get plot
    results$plot.hierEGA <- plot(results, ...)
    
    # Actually send plot
    silent_plot(results$plot.hierEGA)
  
  }

  # Return results
  return(results)
  
}

# Bug checking ----
# data = NetworkToolbox::neoOpen; loading.method = "BRM"
# rotation = "geominQ"; scores = "network"
# scoring.method = "network"; impute = "none"
# corr = "auto"; na.data = "pairwise"; model = "glasso"
# algorithm = "walktrap"; uni.method = "louvain"
# plot.EGA = FALSE; verbose = FALSE

#' @exportS3Method 
# S3 Print Method ----
# Updated 31.07.2023
print.hierEGA <- function(x, ...)
{
  
  # Print lower order
  cat(
    styletext(
      text = styletext(
        text =  "Lower Order\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Print lower network
  print(x$lower_order)
  
  # Add breakspace
  cat("\n\n------------\n\n")
  
  # Print higher order
  cat(
    styletext(
      text = styletext(
        text =  "Higher Order\n\n", 
        defaults = "underline"
      ),
      defaults = "bold"
    )
  )
  
  # Print higher network
  print(x$higher_order)
  
  # Add break space
  cat("\n\n----\n\n")
  
  # Print TEFI
  cat(paste0("Generalized TEFI: ", round(x$TEFI, 3)))
  
}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 17.07.2023
summary.hierEGA <- function(object, ...)
{
  print(object, ...) # same as print
}


#' @exportS3Method 
# S3 Plot Method ----
# Updated 02.08.2023
plot.hierEGA <- function(
    x, plot.type = c("multilevel", "separate"),
    color.match = TRUE, ...
)
{
  
  # Get ellipse arguments for defaults
  ellipse <- list(...)
  
  # Get attributes
  methods_attributes <- attr(x, "methods")
  
  # Check for factor scores
  plot.type <- swiftelse(
    methods_attributes$scores == "factor",
    "separate", # must be separate
    set_default(plot.type, "multilevel", plot.hierEGA)
  )
  
  # Multilevel plot
  if(plot.type == "multilevel"){
    
    # Set edge size
    if(!"edge.size" %in% ellipse){
      ellipse$edge.size <- 8 # default in `basic_plot_setup`
    }
    
    # Get names for levels
    lower_names <- names(x$lower_order$wc)
    higher_names <- names(x$higher_order$wc)
    all_names <- c(lower_names, higher_names)
    
    # Number of total nodes
    total_nodes <- length(all_names)
    
    # Initialize network
    hierarchical_network <- matrix(
      0, nrow = total_nodes, ncol = total_nodes,
      dimnames = list(all_names, all_names)
    )
    
    # Population hybrid network with lower and higher order networks
    hierarchical_network[lower_names, lower_names] <- x$lower_order$network
    hierarchical_network[higher_names, higher_names] <- x$higher_order$network
    
    # Add assignment loading to hierarchical network
    for(community in dimnames(x$parameters$lower_loadings)[[2]]){
      
      # Get assignment loadings
      assignment_loadings <- x$parameters$lower_loadings[
        names(x$lower_order$wc)[x$lower_order$wc == as.numeric(community)],
        community
      ]
      
      # Add to hierarchical network
      hierarchical_network[names(assignment_loadings), community] <- 
        hierarchical_network[community, names(assignment_loadings)] <- 
        assignment_loadings
      
    }
    
    # Get digits based on lower order communities (includes `NA`)
    places <- digits(length(unique(x$lower_order$wc))) - 1
    
    # Set up names for higher order memberships
    higher_order_names <- paste0(
      "Higher_", format_integer(x$higher_order$wc, places)
    )
    
    # Replace "NA" names in higher order
    # Forces singleton representation
    higher_order_names[grep("NA", higher_order_names)] <-
      paste0(
        "Higher_", 
        format_integer(
          names(x$higher_order$wc)[is.na(x$higher_order$wc)], places
        )
      )
    
    # Set up plot list as standard `EGA`
    plot_list <- list(
      # Round hierarchical network to 4
      # (same as what's used in `basic_plot_setup`)
      network = round(hierarchical_network, 4),
      wc = c(
        paste0(
          "Lower_", 
          format_integer(x$lower_order$wc, places)
        ),
        higher_order_names
      )
    )
    
    # Set class
    class(plot_list) <- "EGA"
    
    # Create the initial plot
    initial_plot <- plot(plot_list, ..., arguments = TRUE)
    
    # Make a copy of the hierarchical network
    hierarchical_copy <- plot_list$network
    
    # Add assignment loading to hierarchical network
    for(community in dimnames(x$parameters$lower_loadings)[[2]]){
      
      # Get assignment loadings
      assignment_loadings <- x$parameters$lower_loadings[
        names(x$lower_order$wc)[x$lower_order$wc == as.numeric(community)],
        community
      ]
      
      # Add to hierarchical network
      hierarchical_copy[names(assignment_loadings), community] <- 
        hierarchical_copy[community, names(assignment_loadings)] <- 
        1 # Set all values to 1
      
    }
    
    # Update plot list
    plot_list$network <- hierarchical_copy
    
    # Create the second plot
    second_plot <- plot(plot_list, arguments = TRUE)
    
    # Update multilevel edge appearances
    
    # Get edge indices
    edge_index <- second_plot$ARGS$edge.size == ellipse$edge.size
    
    # Set edge color
    edge_color <- initial_plot$ARGS$edge.color
    edge_color[edge_index] <- "grey"
      
    ## Edge line type
    line_type <- initial_plot$ARGS$edge.lty
    line_type[edge_index] <- "dashed"
    
    ## Edge size
    edge_size <- initial_plot$ARGS$edge.size
    edge_size[edge_index] <- edge_size[edge_index] * 0.50
    
    ## Edge alpha
    edge_alpha <- initial_plot$ARGS$edge.alpha
    edge_alpha[edge_index] <- 0.50
    
    # Get mode (layout)
    mode <- scale(initial_plot$ARGS$mode, scale = FALSE)
    
    # Get maximum y-position
    y_position <- max(mode[,2]) + 15 # for some space
    
    # Get higher order indices
    higher_index <- seq_along(higher_names)
    
    # Shift higher order nodes up
    mode[higher_index, 2] <- mode[higher_index, 2] + y_position +
      abs(min(mode[higher_index, 2]))
    
    # Re-scale x-position
    for(index in higher_index){
      
      # Get lower mode position
      lower_mode_position <- match(
        lower_names[x$lower_order$wc == index],
        initial_plot$ARGS$node.label
      )
      
      # Get higher mode position
      higher_mode_position <- match(
        higher_names[index],
        initial_plot$ARGS$node.label
      )
      
      # Update x-position for higher order node
      mode[higher_mode_position, 1] <-
        median(mode[lower_mode_position, 1])
      
    }
    
    # Check for community matching between layers
    if(color.match){
      
      # Get first plot list
      first_plot_list <- silent_plot(
        plot_list,
        mode = mode,
        edge.size = edge_size,
        edge.color = edge_color,
        edge.alpha = edge_alpha,
        edge.lty = line_type,
        node.size = NA,
        ...,
        arguments = TRUE
      )
      
      # Separate plot and arguments
      first_layer <- first_plot_list$network_plot
      plot_ARGS <- first_plot_list$ARGS
      
      # Create copy of node colors for border colors
      border_color <- plot_ARGS$node.color
      
      # Get targets
      higher_order_targets <- match(names(x$higher_order$wc), plot_ARGS$node.label)
      
      # Target higher order colors
      border_color[higher_order_targets] <- unique(
        border_color[match(names(x$lower_order$wc), plot_ARGS$node.label)]
      )
      
      # Set stroke size to be thicker for higher order only
      ## Get stroke size
      stroke_size <- rep(
        first_layer$guides$colour$override.aes$stroke, total_nodes
      )
      ## Adjust higher order only
      stroke_size[higher_order_targets] <- 2.5

      # Custom nodes: transparent insides and dark borders
      second_layer <- first_layer +
        ggplot2::geom_point( # transparent insides
          size = second_plot$ARGS$node.size + 0.50, shape = 19,
          color = plot_ARGS$node.color,
          alpha = plot_ARGS$node.alpha,
          show.legend = FALSE
        ) +
        ggplot2::geom_point( # dark borders
          size = second_plot$ARGS$node.size,
          color = border_color,
          shape = 1, stroke = stroke_size, alpha = 1
        ) +
        ggplot2::geom_text( # put text back on top
          ggplot2::aes(label = first_layer$data$label), 
          color = "black",
          size = plot_ARGS$label.size
        ) +
        ggplot2::guides( # create legend with these settings
          color = ggplot2::guide_legend(
            override.aes = list(
              shape = 21,
              fill = unique(plot_ARGS$node.color),
              size = median(second_plot$ARGS$node.size, na.rm = TRUE),
              alpha = median(plot_ARGS$node.alpha, na.rm = TRUE),
              stroke = 1.5
            ),
            title = swiftelse(
              "legend.title" %in% names(ellipse),
              ellipse$legend.title, ""
            )
          )
        )
      
      # Get current warning options
      current_warning <- getOption("warn")
      
      # Set them off
      options(warn = -1)
      
      # Send it
      silent_plot(second_layer)
      
      # Set them back to the user's setting
      options(warn = current_warning)
      
    }else{
      silent_plot(
        plot_list,
        mode = mode,
        edge.size = edge_size,
        edge.color = edge_color,
        edge.alpha = edge_alpha,
        edge.lty = line_type,
        ...
      )
    }

  }else if(plot.type == "separate"){ # Separate plot
    
    # Set labels
    if(!"labels" %in% ellipse){
      ellipse$labels <- c("Lower", "Higher")
    }
    
    # Plot lower and higher order side-by-side
    ggpubr::ggarrange(
      silent_plot(x$lower_order, ...),
      silent_plot(x$higher_order, ...),
      labels = ellipse$labels,
      ...
    )
    
  }
  
}

#' @noRd
# Oblimin fit ----
# Updated 15.07.2023
oblimin <- function(L, gamma = 0)
{
  
  # Get dimensions of loadings
  dimensions <- dim(L)
  
  # Check for special case of gamma = 0
  if(gamma == 0){ # Set IgC matrix as identity
    IgC <- diag(dimensions[1])
  }else{
    
    # Get gC matrix
    gC <- matrix(
      gamma / dimensions[1],
      nrow = dimensions[1], ncol = dimensions[1]
    )
    
    # Subtract gC from identity matrix
    IgC <- diag(dimensions[1]) - gC
    
  }
  
  # Initialize N matrix
  N <- matrix(1, nrow = dimensions[2], ncol = dimensions[2])
  
  # Set diagonal to zero
  diag(N) <- 0
  
  # Square loadings matrix
  L2 <- L^2 # nanoseconds faster than L * L
  
  # Return f
  return(sum(diag(crossprod(L2, IgC %*% L2 %*% N))) * 0.25)
  
}

#' @noRd
# Oblimin rotation ----
# Updated 15.07.2023
oblimin_rotate <- function(
    loadings, nfactors, 
    n.rotations = 10, maxit = 10000,
    eps = 1e-5, gamma = 0,
    rotate = "oblimin"
)
{
  
  # Check for unidimensional structure (28.11.2022)
  if(nfactors == 1){
    return(list(loadings = loadings, Phi = 1))
  }
  
  # Get factor sequence
  factor_sequence <- seq_len(nfactors)
  
  # Initialize minimum fit value
  fit_value <- Inf
  
  # Perform loop
  for(i in seq_len(n.rotations)){
    
    # Random values
    X <- nvapply(factor_sequence, function(x){rnorm(nfactors)}, LENGTH = nfactors)
    
    # Loading results
    result <- GPArotation::GPFoblq(
      A = loadings, method = rotate, Tmat = qr.Q(qr(X)),
      maxit = maxit, eps = eps
    )
    
    # Check for improvement in fit
    if(oblimin(result$loadings, gamma = gamma) < fit_value){
      return_result <- result[c("loadings", "Phi")]
    }
    
  }
  
  # Return results
  return(return_result)
  
}

#' @noRd
# Exploratory Factor Analysis Scores ----
# Updated 15.07.2023
efa_scores <- function(
    data, correlation_matrix, nfactors, 
    fm = "minres", rotate = "oblimin",
    n.rotations = 10, maxit = 1e4, 
    factor.scores = "Thurstone", gamma = 0,
    eps = 1e-5, impute
)
{

  # Start with {psych}'s factor analysis
  fa <- psych::fa(
    correlation_matrix, 
    nfactors = nfactors, 
    fm = fm, rotate = "none"
  )
  
  # Perform oblimin rotation
  rotation <- oblimin_rotate(
    loadings = fa$loadings, nfactors = nfactors,
    n.rotations = n.rotations, maxit = maxit,
    eps = eps, gamma = gamma, rotate = rotate
  )

  # Return result
  return(
    list(
      loadings = rotation$loadings,
      scores = psych::factor.scores(
        x = data, f = rotation$loadings,
        Phi = rotation$Phi, method = factor.scores,
        rho = correlation_matrix, missing = TRUE,
        impute = impute
      )$scores
    )
  )

  # At this point, only scores are used and reported
  # return(
  #   list(
  #     loadings = rotation$loadings, Phi = rotation$Phi,
  #     Shat = fa$model, uniquenesses = fa$uniquenesses,
  #     factor.scores = psych::factor.scores(
  #       x = data, f = rotation$loadings,
  #       Phi = rotation$Phi, method = factor.scores,
  #       rho = correlation_matrix, missing = TRUE,
  #       impute = impute
  #     )
  #   )
  # )

}



