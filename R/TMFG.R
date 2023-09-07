#' @title Triangulated Maximally Filtered Graph
#'
#' @description Applies the Triangulated Maximally Filtered Graph (TMFG) filtering method
#' (see Massara et al., 2016). The TMFG method uses a structural
#' constraint that limits the number of zero-order correlations included in the network
#' (3\emph{n} - 6; where \emph{n} is the number of variables). The TMFG algorithm begins by
#' identifying four variables which have the largest sum of correlations to all other
#' variables. Then, it iteratively adds each variable with the largest sum of three
#' correlations to nodes already in the network until all variables have been added to
#' the network. This structure can be associated with the inverse correlation matrix
#' (i.e., precision matrix) to be turned into a GGM (i.e., partial correlation network)
#' by using Local-Global Inversion Method (LoGo; see Barfuss et al., 2016 for more details).
#' See \strong{Details} for more information
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or correlation matrix
#'
#' @param n Numeric (length = 1).
#' Sample size for when a correlation matrix is input into \code{data}.
#' Defaults to \code{NULL}.
#' \code{n} is not necessary and is provided for better functionality in
#' \code{\link{EGAnet}}
#' 
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"auto"} --- }
#' {Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)}
#' 
#' \item{\code{"cor_auto"} --- }
#' {Uses \code{\link[qgraph]{cor_auto}} to compute correlations. Arguments
#' can be passed along to the function}
#' 
#' \item{\code{"pearson"} --- }
#' {Pearson's correlation is computed for all variables regardless of
#' categories}
#' 
#' \item{\code{"spearman"} --- }
#' {Spearman's rank-order correlation is computed for all variables
#' regardless of categories}
#' 
#' }
#' 
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#' 
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"} --- }
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"} --- }
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param partial Boolean (length = 1).
#' Whether partial correlations should be output.
#' Defaults to \code{FALSE}.
#' The TMFG method is based on the zero-order correlations;
#' the Local-Global Inversion Method (LoGo; see Barfuss et al., 2016 for more details)
#' uses the decomposability of the TMFG network to obtain the inverse covariance
#' structure of the network (which is then converted to partial correlations).
#' Set to \code{TRUE} to obtain the partial correlations from the LoGo method
#' 
#' @param returnAllResults Boolean (length = 1).
#' Whether all results should be returned.
#' Defaults to \code{FALSE} (network only).
#' Set to \code{TRUE} to access separators and cliques
#' 
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}}
#'
#' @return Returns a network or list containing:
#'
#' \item{network}{The filtered adjacency matrix}
#'
#' \item{separators}{The separators (3-cliques) in the network}
#'
#' \item{cliques}{The cliques (4-cliques) in the network}
#'
#' @details The TMFG method applies a structural constraint on the network,
#' which restrains the network to retain a certain number of edges (3\emph{n}-6, where \emph{n}
#' is the number of nodes; Massara et al., 2016). The network is also composed of 3- and 4-node
#' cliques (i.e., sets of connected nodes; a triangle and tetrahedron, respectively). The
#' TMFG method constructs a network using zero-order correlations and the resulting network
#' can be associated with the inverse covariance matrix
#' (yielding a GGM; Barfuss, Massara, Di Matteo, & Aste, 2016).
#' Notably, the TMFG can use any association measure and thus does not assume the data is multivariate normal.
#'
#' Construction begins by forming a tetrahedron of the four nodes that have
#' the highest sum of correlations that are greater than the average correlation in the
#' correlation matrix. Next, the algorithm iteratively identifies the node that maximizes
#' its sum of correlations to a connected set of three nodes (triangles) already included
#' in the network and then adds that node to the network. The process is completed once
#' every node is connected in the network. In this process, the network automatically
#' generates what's called a planar network. A planar network is a network that could be
#' drawn on a sphere with no edges crossing (often, however, the networks are depicted
#' with edges crossing; Tumminello, Aste, Di Matteo, & Mantegna, 2005).
#'
#' @examples
#' # TMFG filtered network
#' TMFG(wmt2[,7:24])
#' 
#' # Partial correlations using the LoGo method
#' TMFG(wmt2[,7:24], partial = TRUE)
#'
#' @references
#' \strong{Local-Global Inversion Method} \cr
#' Barfuss, W., Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Parsimonious modeling with information filtering networks.
#' \emph{Physical Review E}, \emph{94}, 062306.
#'
#' \strong{Psychometric network introduction to TMFG} \cr
#' Christensen, A. P., Kenett, Y. N., Aste, T., Silvia, P. J., & Kwapil, T. R. (2018).
#' Network structure of the Wisconsin Schizotypy Scales-Short Forms: Examining psychometric network filtering approaches.
#' \emph{Behavior Research Methods}, \emph{50}, 2531-2550.
#'
#' \strong{Triangulated Maximally Filtered Graph} \cr
#' Massara, G. P., Di Matteo, T., & Aste, T. (2016).
#' Network filtering for big data: Triangulated maximally filtered graph.
#' \emph{Journal of Complex Networks}, \emph{5}, 161-178.
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
# TMFG Filtering Method----
# Updated 07.09.2023
TMFG <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    partial = FALSE, returnAllResults = FALSE,
    verbose = FALSE, 
    ...
)
{
  
  # Argument errors (return data in case of tibble)
  data <- TMFG_errors(data, n, partial, returnAllResults, verbose, ...)
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", TMFG)
  na.data <- set_default(na.data, "pairwise", TMFG)

  # Make sure there are variable names
  data <- ensure_dimension_names(data)
  
  # Generic function to get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = 1, # "n" is not used but input `1` to avoid error
    corr = corr, na.data = na.data, 
    verbose = verbose, needs_usable = FALSE, # skips usable data check
    ...
  )
  
  # Get correlations
  correlation_matrix <- output$correlation_matrix
  
  # Obtain number of nodes
  nodes <- dim(correlation_matrix)[2]
  
  # Set warning for fewer than 9 nodes
  if(nodes < 9 & isTRUE(verbose)){
    warning(
      "The TMFG method requires more than 9 nodes to obtain a chordal network. With fewer than 9 nodes, this property does not hold",
      call. = FALSE
    )
  }
  
  # For signed correlations, use absolute correlation matrix
  # and obtain element-wise inclusion of signed correlation matrix later
  absolute_matrix <- abs(correlation_matrix)
  
  # Initialize inserted nodes vector
  inserted <- numeric(nodes)
  
  # Separator rows
  separator_rows <- nodes - 4
  
  # Initialize triangles and separators matrix
  triangles <- matrix(nrow = 2 * nodes - 4, ncol = 3)
  separators <- matrix(nrow = separator_rows, ncol = 3)
  
  # Obtain four nodes with the largest strength
  # which is greater than the average strength
  ## Compute node strength
  node_strength <- colSums(absolute_matrix, na.rm = TRUE)
  ## Select the four nodes with the largest strength
  four_nodes <- colSums(
    absolute_matrix * (absolute_matrix > mean(absolute_matrix, na.rm = TRUE)),
    na.rm = TRUE
  )
  
  # First four nodes
  first_four <- seq_len(4)
  
  # Insert the top four nodes
  inserted[first_four] <- order(four_nodes, decreasing = TRUE)[first_four]
  
  # Set remaining nodes
  remaining <- setdiff(seq_len(nodes), inserted)
  
  # Build tetrahedron
  triangles[1,] <- inserted[seq_len(3)]; triangles[2,] <- inserted[2:4];
  triangles[3,] <- inserted[c(1, 2, 4)]; triangles[4,] <- inserted[c(1, 3, 4)];
  
  # Initialize network (correlations to retain)
  network <- diag(1, nrow = nodes, ncol = nodes)
  
  # Add nodes to network
  network[inserted[first_four], inserted[first_four]] <-
    correlation_matrix[inserted[first_four], inserted[first_four]]
  
  # Build gain table
  gain_columns <- 2 * (nodes - 2)
  gain <- matrix(-Inf, nrow = nodes, ncol = gain_columns)
  gain[remaining, 1] <- rowSums(absolute_matrix[remaining, triangles[1,]], na.rm = TRUE)
  gain[remaining, 2] <- rowSums(absolute_matrix[remaining, triangles[2,]], na.rm = TRUE)
  gain[remaining, 3] <- rowSums(absolute_matrix[remaining, triangles[3,]], na.rm = TRUE)
  gain[remaining, 4] <- rowSums(absolute_matrix[remaining, triangles[4,]], na.rm = TRUE)
  
  # Number of triangles
  triangle_count <- 4
  gain_vertex <- numeric(gain_columns)
  gain_weight <- numeric(gain_columns)
  
  # Loop over remaining edges
  for(i in 5:nodes){
  
    # Check for one remaining node
    if(length(remaining) == 1){
      
      # Last vertex to add
      add_vertex <- remaining; existing_vertex <- 1;
      gain_vertex <- 1; max_gain <- which.max(gain[remaining,]);
      
    }else{
      
      # Get max gains
      max_gains <- gain[remaining,] == max(gain[remaining,])
      
      # Get total gains
      total_gains <- sum(max_gains)

      # Check for more than 1
      if(total_gains > 1){
        max_gains[max_gains][2:total_gains] <- FALSE
      }
      
      # Vectorized solution (avoids nested loop)
      gain_location <- which(max_gains, arr.ind = TRUE)
      
      # Obtain existing vertex (already in network)
      # and vertex from remaining based on the maximum gain
      max_gain <- gain_location[,"col"]
      existing_vertex <- gain_location[,"row"]
      add_vertex <- remaining[existing_vertex]
      
    }
    
    # Update lists
    remaining <- remaining[-existing_vertex]
    inserted[i] <- add_vertex
    
    # Add edges to network
    network[add_vertex, triangles[max_gain,]] <-
      correlation_matrix[add_vertex, triangles[max_gain,]]
    
    # Add to other side
    network[triangles[max_gain,], add_vertex] <-
      correlation_matrix[triangles[max_gain,], add_vertex]
    
    # Update separators
    separators[i-4,] <- triangles[max_gain,]
    
    # Update triangles list
    ## Add two triangles
    triangles[triangle_count + 1,] <- c(triangles[max_gain, c(1, 3)], add_vertex)
    triangles[triangle_count + 2,] <- c(triangles[max_gain, c(2, 3)], add_vertex)
    ## Replace maximum gain triangle (no longer a triangle)
    triangles[max_gain,] <- c(triangles[max_gain, c(1, 2)], add_vertex)
    
    # Update gain table
    ## Set gain to zero for added node
    gain[add_vertex,] <- 0
    ## Maximum gain
    gain[remaining, max_gain] <- rowSums(
      absolute_matrix[remaining, triangles[max_gain,], drop = FALSE], 
      na.rm = TRUE
    )
    ## First new triangle
    gain[remaining, triangle_count + 1] <- rowSums(
      absolute_matrix[remaining, triangles[triangle_count + 1,], drop = FALSE], 
      na.rm = TRUE
    )
    ## Second new triangle
    gain[remaining, triangle_count + 2] <- rowSums(
      absolute_matrix[remaining, triangles[triangle_count + 2,], drop = FALSE], 
      na.rm = TRUE
    )
    
    ## Increase triangle count
    triangle_count <- triangle_count + 2
    
  }
  
  # Create cliques
  cliques <- rbind(
    inserted[first_four],
    cbind(separators, inserted[5:nodes])
  )
  
  # Check for whether partial correlation network should be computed
  ## An extension of the TMFG method using the LoGo method (Barfuss et al., 2016)
  if(partial){
    
    # Initialize partial correlation network 
    partial_network <- matrix(0, nrow = nodes, ncol = nodes)
    
    # Loop over cliques and separators
    for(i in seq_len(separator_rows)){
      
      # Obtain clique
      clique <- cliques[i,]
      
      # Add clique
      partial_network[clique, clique] <-
        partial_network[clique, clique] +
        solve(correlation_matrix[clique, clique])
      
      # Obtain separator
      separator <- separators[i,]
        
      # Subtract separator
      partial_network[separator, separator] <-
        partial_network[separator, separator] -
        solve(correlation_matrix[separator, separator])
      
    }
    
    # There is always one more clique than separator
    ## Obtain last clique
    clique <- cliques[separator_rows + 1,]
    
    ## Add last clique
    partial_network[clique, clique] <-
      partial_network[clique, clique] +
      solve(correlation_matrix[clique, clique])
    
    # Convert inverse covariance matrix to partial correlations
    ## Replaces the original network
    network <- -cov2cor(partial_network)
    diag(network) <- 0
    
  }
  
  # Transfer variable names
  network <- transfer_names(data, network)
  
  # Set methods attribute
  attr(network, "methods") <- list(
    corr = corr, partial = partial
  )
  
  # Set up return
  if(!returnAllResults){
    return(network)
  }else{
    
    # Set up return list
    return(
      list(
        network = network,
        separators = separators,
        cliques = cliques
      )
    )
    
  }
  
}

# Bug Checking ----
# ## Basic input
# data = wmt2[,7:24]; n = NULL
# corr = "auto"; na.data = "pairwise"
# partial = FALSE; returnAllResults = FALSE
# verbose = FALSE

#' @noRd
# Errors ----
# Updated 07.09.2023
TMFG_errors <- function(data, n, partial, returnAllResults, verbose, ...)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "TMFG")
  
  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }
  
  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "TMFG")
    typeof_error(n, "numeric", "TMFG")
  }
  
  # 'partial' errors
  length_error(partial, 1, "TMFG")
  typeof_error(partial, "logical", "TMFG")
  
  # 'returnAllResults' errors
  length_error(returnAllResults, 1, "TMFG")
  typeof_error(returnAllResults, "logical", "TMFG")
  
  # 'verbose' errors
  length_error(verbose, 1, "TMFG")
  typeof_error(verbose, "logical", "TMFG")
  
  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }
  
  # Return usable data in case of tibble
  return(data)
  
}

