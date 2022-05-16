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
#' {Applies the leading eigenvalue algorithm (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the leading eigenvalue solution is used; otherwise, regular EGA
#' is used. This is the final method used in the Christensen, Garrido,
#' and Golino (2021) simulation.}
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
#' Defaults to \code{"louvain"}.
#' Current options are:
#'
#' \itemize{
#'
#' \item{\strong{\code{walktrap}}}
#' {Computes the Walktrap algorithm using \code{\link[igraph]{cluster_walktrap}}}
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
#' @param plot.EGA Boolean.
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#' @param plot.type Character.
#' Plot system to use.
#' Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
#' Defaults to \code{"GGally"}
#'
#' @param plot.args List.
#' A list of additional arguments for the network plot.
#' For \code{plot.type = "qgraph"}:
#'
#' \itemize{
#'
#' \item{\strong{\code{vsize}}}
#' {Size of the nodes. Defaults to 6.}
#'
#'}
#' For \code{plot.type = "GGally"} (see \code{\link[GGally]{ggnet2}} for
#' full list of arguments):
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
#' # hierEGA example (no plots)
#' opt.hier<- hierEGA(
#'   data = optimism,
#'   algorithm = "louvain",
#'   plot.EGA = FALSE # no plots for CRAN check
#' )
#'
#' \dontrun{
#'
#' # hierEGA example (plots)
#' opt.hier <- hierEGA(
#'   data = optimism,
#'   algorithm = "louvain"
#' )
#'
#' # Save plots
#' ggplot2::ggsave(
#'   filename = "hierEGA_plot.png", # name of plot
#'   plot = opt.res$hierarchical$hier_plot, # plot to save
#'   height = 8, # figure height
#'   width = 10, # figure width
#'   dpi = 600 # dots per inch
#' )
#'
#' }
#'
#'
#' @export
#'
# Hierarchical EGA
# Updated 16.05.2022
hierEGA <- function(
    data, scores = c("factor", "network"),
    consensus.iter = 1000,
    consensus.method = c(
      "highest_modularity",
      "most_common",
      "iterative",
      "lowest_tefi"
    ),
    uni.method = c("expand", "LE"),
    corr = c("cor_auto", "pearson", "spearman"),
    model = c("glasso", "TMFG"), model.args = list(),
    algorithm = c("walktrap", "louvain"), algorithm.args = list(),
    plot.EGA = TRUE, plot.type = c("GGally", "qgraph"),
    plot.args = list()
)
{

  #### ARGUMENTS HANDLING ####

  if(missing(scores)){
    scores <- "network"
  }else{scores <- match.arg(scores)}

  if(missing(consensus.method)){
    consensus.method <- "highest_modularity"
  }else{consensus.method <- match.arg(consensus.method)}

  if(missing(uni.method)){
    uni.method <- "LE"
  }else{uni.method <- match.arg(uni.method)}

  if(missing(corr)){
    corr <- "cor_auto"
  }else{corr <- match.arg(corr)}

  if(missing(model)){
    model <- "glasso"
  }else{model <- match.arg(model)}

  if(missing(algorithm)){
    algorithm <- "louvain"
  }else if(!is.function(algorithm)){algorithm <- tolower(match.arg(algorithm))}

  if(missing(plot.type)){
    plot.type <- "GGally"
  }else{plot.type <- match.arg(plot.type)}

  #### ARGUMENTS HANDLING ####

  # Ensure data is a matrix
  data <- as.matrix(data)

  # Check for symmetric matrix
  if(isSymmetric(data)){
    stop(
      paste(
        "Symmetric matrix detected. Data with variables down the columns and cases across the rows are required to compute",
        scores,
        "scores for this analysis."
      )
    )
  }

  # Obtain EGA defaults
  ega_defaults <- formals(EGA)

  # Remove "..."
  ega_defaults <- ega_defaults[-length(ega_defaults)]

  # Start lower-order ----

  # Set lower-order defaults
  lower_order_defaults <- ega_defaults
  lower_order_defaults$data <- data
  lower_order_defaults$uni.method <- uni.method
  lower_order_defaults$corr <- corr
  lower_order_defaults$model <- model
  lower_order_defaults$model.args <- model.args
  lower_order_defaults$algorithm <- "louvain" # for lower order communities
  lower_order_defaults$algorithm.args <- algorithm.args
  lower_order_defaults$plot.lower_order <- FALSE # do not plot
  lower_order_defaults$plot.type <- plot.type
  lower_order_defaults$plot.args <- plot.args
  lower_order_defaults$lower.louvain <- TRUE # provides lower order Louvain

  # Send message
  message(
    "Obtaining lower-order dimensions...",
    appendLF = FALSE
  )

  # Get EGA
  lower_order_result <- suppressMessages(
    do.call(
      EGA.estimate, lower_order_defaults
    )
  )

  # Perform consensus clustering
  lower_order_wc <- consensus_clustering(
    network = lower_order_result$network,
    corr = lower_order_result$cor.data,
    order = "lower",
    consensus.iter = consensus.iter
  )

  # End message
  message("done.")

  # Estimate scores
  ## Factor ----
  # Set up for all three results
  factor_results <- vector("list", length = length(lower_order_wc) - 1)
  names(factor_results) <- names(lower_order_wc)[
    -length(lower_order_wc)
  ]

  # Loop through
  for(i in 1:length(factor_results)){

    # Obtain memberships
    memberships <- lower_order_wc[[names(factor_results)[i]]]
    unique_memberships <- unique(memberships)

    # Estimate factor model
    fm <- suppressWarnings(
      psych::fa(
        r = lower_order_result$cor.data, # correlation matrix
        n.obs = nrow(data),
        nfactors = length(na.omit(unique_memberships)) # number of factors
      )
    )

    # Score estimates
    score_est <- psych::factor.scores(
      x = data,
      f = fm,
      method = "tenBerge"
    )$scores

    # Lower-order loadings
    lower_loads <- fm$loadings[,1:length(na.omit(unique_memberships))]

    # Reorder lower loadings and scores
    lower_loads <- lower_loads[,order(colnames(lower_loads))]
    score_est <- score_est[,order(colnames(score_est))]
    
    # Obtain highest loading to reorder rows
    highest_loads <- apply(abs(lower_loads), 1, which.max)
    highest_loads[1:length(highest_loads)] <- paste("MR", highest_loads, sep = "")
    lower_loads <- descend.ord(lower_loads, highest_loads)

    # Start higher-order with factor ----

    # Set the rest of the arguments
    ega_defaults$data <- score_est # set data as score estimates
    ega_defaults$uni.method <- uni.method
    ega_defaults$corr <- corr
    ega_defaults$model <- model
    ega_defaults$model.args <- model.args
    ega_defaults$algorithm <- "walktrap"
    ega_defaults$algorithm.args <- algorithm.args
    ega_defaults$plot.EGA <- FALSE
    ega_defaults$plot.type <- plot.type
    ega_defaults$plot.args <- plot.args

    # Get EGA
    higher_order_result <- suppressWarnings(
      suppressMessages(
        do.call(
          EGA, ega_defaults
        )
      )
    )

    # Make consensus
    higher_order_wc <- consensus_clustering(
      network = higher_order_result$network,
      corr = higher_order_result$correlation,
      order = "higher",
      consensus.iter = consensus.iter
    )[[names(factor_results)[i]]]

    # Set up result
    factor_results[[names(factor_results)[i]]]$lower_scores <- round(score_est, 3)
    factor_results[[names(factor_results)[i]]]$lower_loadings <- round(lower_loads, 3)

    ## Walktrap
    factor_results[[names(factor_results)[i]]]$walktrap <- higher_order_result

    ## Louvain
    ## Adjust memberships in higher order
    higher_order_result$wc <- higher_order_wc
    higher_order_result$n.dim <- length(na.omit(unique(higher_order_wc)))
    higher_order_result$dim.variables[,"dimension"] <- higher_order_wc[
      higher_order_result$dim.variables[,"items"]
    ]
    factor_results[[names(factor_results)[i]]]$louvain <- higher_order_result

  }

  # Skip network scores?
  ## Identify memberships with only one node
  one_node <- unlist(lapply(lower_order_wc[-length(lower_order_wc)], function(x){
    
    ### Obtain frequencies
    wc_freq <- table(x)
    
    ### Check for one node memberships
    if(any(wc_freq <= 1)){
      one <- TRUE
    }else{
      one <- FALSE
    }
    
    ### Return
    return(one)
    
  }))
  
  # Perform network scores
  perform_network <- TRUE
  
  ## Check for all
  if(all(one_node)){
    
    ### Warn that network scores cannot be used
    warning("All consensus methods detected a single node community. Network scores cannot be computed.")
    
    ### Check if scores was what user wanted
    if(scores == "network"){
      stop("Network scores could not be computed. Please use factor scores.")
    }else{
      perform_network <- FALSE
    }
    
  }else if(any(one_node)){
    
    # One node consensus methods
    one_node_names <- names(lower_order_wc[-length(lower_order_wc)])[one_node]
    
    warning(
      paste(
        "Some consensus methods detected a single node community. Network scores could not be computed for:\n",
        paste(one_node_names, collapse = "\n"),
        sep = ""
      )
    )
    
    # Remove from list
    lower_order_wc <- lower_order_wc[-match(one_node_names, names(lower_order_wc))]
    
    # Check if method was what user wanted
    if(one_node_names %in% consensus.method){
      stop(
        paste(
          "Cannot use consensus method due to at least one single node community. Please use one of the following methods:\n",
          paste(
            setdiff(c(
              "highest_modularity", "most_common",
              "iterative", "lowest_tefi"
            ), one_node_names),
            collapse = "\n"
          ),
          sep = ""
        )
      )
    }
    
  }
  
  if(isTRUE(perform_network)){
    
    ## Network ----
    # Set up for all three results
    network_results <- vector("list", length = length(lower_order_wc) - 1)
    names(network_results) <- names(lower_order_wc)[
      -length(lower_order_wc)
    ]
    
    # Loop through
    for(i in 1:length(network_results)){
      
      # Obtain memberships
      memberships <- lower_order_wc[[names(network_results)[i]]]
      
      # Compute network scores
      nt <- suppressWarnings(
        net.scores(
          data = data,
          A = lower_order_result$network,
          wc = memberships
        )
      )
      
      # Score estimates
      score_est <- nt$std.scores
      
      # Lower-order loadings
      lower_loads <- nt$loads
      
      # Start higher-order with network ----
      
      # Set the rest of the arguments
      ega_defaults$data <- score_est # set data as score estimates
      ega_defaults$uni.method <- uni.method
      ega_defaults$corr <- corr
      ega_defaults$model <- model
      ega_defaults$model.args <- model.args
      ega_defaults$algorithm <- "walktrap"
      ega_defaults$algorithm.args <- algorithm.args
      ega_defaults$plot.EGA <- FALSE
      ega_defaults$plot.type <- plot.type
      ega_defaults$plot.args <- plot.args
      
      # Get EGA
      higher_order_result <- suppressWarnings(
        suppressMessages(
          do.call(
            EGA, ega_defaults
          )
        )
      )
      
      # Make consensus
      higher_order_wc <- consensus_clustering(
        network = higher_order_result$network,
        corr = higher_order_result$correlation,
        order = "higher",
        consensus.iter = consensus.iter
      )[[names(network_results)[i]]]
      
      # Set up result
      network_results[[names(network_results)[i]]]$lower_scores <- round(score_est, 3)
      network_results[[names(network_results)[i]]]$lower_loadings <- round(lower_loads, 3)
      
      ## Walktrap
      network_results[[names(network_results)[i]]]$walktrap <- higher_order_result
      
      ## Louvain
      ## Adjust memberships in higher order
      higher_order_result$wc <- higher_order_wc
      higher_order_result$n.dim <- length(na.omit(unique(higher_order_wc)))
      higher_order_result$dim.variables[,"dimension"] <- higher_order_wc[
        higher_order_result$dim.variables[,"items"]
      ]
      network_results[[names(network_results)[i]]]$louvain <- higher_order_result
      
    }
    
  }else{
    network_results <- "Did not perform network scores due to at least one single node community in all consensus methods"
  }

  # Return results
  results <- list()
  results$lower_ega <- lower_order_result
  results$lower_wc <- lower_order_wc
  results$factor_results <- factor_results
  results$network_results <- network_results

  # Set up return results
  hierarchical <- list()

  # Obtain lower order results
  lower_order_result$wc <- lower_order_wc[[consensus.method]]
  lower_order_result$n.dim <- length(na.omit(unique(lower_order_result$wc)))

  # Get S3 print information for lower order
  if(is.null(colnames(data))){
    dim.variables <- data.frame(
      items = paste("V", 1:ncol(data), sep = ""),
      dimension = lower_order_result$wc
    )
  }else{
    dim.variables <- data.frame(
      items = colnames(data),
      dimension = lower_order_result$wc
    )
  }
  dim.variables <- dim.variables[order(dim.variables[, 2]),]
  lower_order_result$dim.variables <- dim.variables

  # Set class for lower-order
  class(lower_order_result) <- "EGA"

  # Insert into hierarchical result
  hierarchical$lower_order <- lower_order_result

  # Obtain higher order result
  if(scores == "factor"){

    hierarchical$higher_order <- factor_results[[consensus.method]][
      c("lower_scores", "lower_loadings", algorithm)
    ]

    names(hierarchical$higher_order)[
      which(names(hierarchical$higher_order) == algorithm)
    ] <- "EGA"

  }else if(scores == "network"){

    hierarchical$higher_order <- network_results[[consensus.method]][
      c("lower_scores", "lower_loadings", algorithm)
    ]

    names(hierarchical$higher_order)[
      which(names(hierarchical$higher_order) == algorithm)
    ] <- "EGA"

  }

  # Insert hierarchical result into results
  results$hierarchical <- hierarchical

  # Obtain higher order result for plot
  higher_order_result <- hierarchical$higher_order$EGA

  # Set up plots
  if(isTRUE(plot.EGA)){

    # Set up plots
    lower_plot <- suppressPackageStartupMessages(
      plot(lower_order_result, produce = FALSE)
    )
    higher_plot <- suppressPackageStartupMessages(
      plot(higher_order_result, produce = FALSE)
    )

    # Set up output
    hier_plot <- suppressPackageStartupMessages(
      ggpubr::ggarrange(
        lower_plot, # plot lower-order
        higher_plot, # plot higher-order
        labels = c("Lower-order", "Higher-order")
      )
    )

    # Output plots
    suppressPackageStartupMessages(
      plot(hier_plot)
    )

    # Add to main results
    results$hierarchical$lower_plot <- lower_plot
    results$hierarchical$higher_plot <- higher_plot
    results$hierarchical$hier_plot <- hier_plot

    # Send factor warning
    if(scores == "factor"){

      warning("Lower order factor loadings may not map to lower order network dimensions.\nPlease see `$lower_loadings` to map lower order dimensions to higher order dimensions in the plot")

    }

  }

  # Reorder results so main results are first
  results <- results[
    c("hierarchical", "lower_ega", "lower_wc",
      "factor_results", "network_results")
  ]

  # Set up Methods
  methods <- list(
    corr = corr, model = model, model.args = model.args,
    algorithm = algorithm, algorithm.args = algorithm.args,
    uni.method = uni.method, scores = scores,
    consensus.method = consensus.method, consensus.iter = consensus.iter
  )

  # Insert methods into hierarhical
  results$hierarchical$Methods <- methods

  # Make class "hierEGA"
  class(results) <- "hierEGA"

  return(results)
}
