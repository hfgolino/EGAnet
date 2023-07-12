#' Network Loadings
#'
#' @description Computes the between- and within-community
#' \code{strength} of each item
#' for each community. This function uses the
#' \code{comcat} and
#' \code{stable} functions to calculate
#' the between- and within-community strength of each item, respectively.
#'
#' @param A Matrix, data frame, or \code{\link[EGAnet]{EGA}} object.
#' A network adjacency matrix
#'
#' @param wc Numeric or character vector.
#' A vector of community assignments.
#' If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
#' then \code{wc} is automatically detected
#' 
#' @param rotation Character.
#' A rotation to use, like factor loadings, to obtain
#' a simple structure.
#' Defaults to \code{\link[GPArotation]{geominQ}}.
#' For a list of rotations, see \code{\link{GPArotation}}
#' 
#' @param min.load Numeric.
#' Sets the minimum loading allowed in the standardized
#' network loading matrix. Values equal or greater than
#' the minimum loading are kept in the output. Values
#' less than the minimum loading are removed. This matrix can
#' be viewed using \code{print()} or \code{summary()}
#' Defaults to \code{0}
#'
#' @return Returns a list containing:
#'
#' \item{unstd}{A matrix of the unstandardized within- and between-community
#' strength values for each node}
#'
#' \item{std}{A matrix of the standardized within- and between-community
#' strength values for each node}
#' 
#' \item{minLoad}{The minimum loading to appear in summary of network loadings.
#' Use \code{print()} or \code{summary()} to view}
#'
#' @details Simulation studies have demonstrated that a node's strength
#' centrality is roughly equivalent to factor loadings
#' (Christensen, Golino, & Silvia, 2019; Hallquist, Wright, & Molenaar, in press).
#' Hallquist and colleagues (in press) found that node strength represented a
#' combination of dominant and cross-factor loadings. This function computes
#' each node's strength within each specified dimension, providing a rough
#' equivalent to factor loadings (including cross-loadings).
#'
#' For more details, type \code{vignette("Network_Scores")}
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )}
#'
#' # Network loadings
#' net.loads(ega.wmt)
#' 
#' \dontrun{
#' # Produce Methods section
#' methods.section(
#'   ega.wmt,
#'   stats = "net.loads"
#' )}
#'
#' @references
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}, \emph{53}, 1563-1580.
#' 
#' Christensen, A. P., Golino, H., & Silvia, P. J. (2020).
#' A psychometric network perspective on the validity and validation of personality trait questionnaires.
#' \emph{European Journal of Personality}, \emph{34}, 1095-1108.
#'
#' Hallquist, M., Wright, A. C. G., & Molenaar, P. C. M. (2019).
#' Problems with centrality measures in psychopathology symptom networks: Why network psychometrics cannot escape psychometric theory.
#' \emph{Multivariate Behavioral Research}, 1-25.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Loadings
# Updated 12.07.2023
# Default = "BRM" or `net.loads` from version 1.2.3
# Experimental = new signs and cross-loading adjustment
net.loads <- function(
    A, wc, loading.method = c("BRM", "experimental"),
    rotation = NULL, ...
)
{
  
  # Check for missing arguments (argument, default, function)
  loading.method <- set_default(loading.method, "brm", net.loads)
  
  # Organize and extract input
  input <- organize_input(A, wc)
  A <- input$A; wc <- input$wc
  
  # Get number of nodes and their names
  nodes <- length(wc); node_names <- names(wc)
  
  # Get unique communities (`NA` is OK)
  unique_communities <- sort(unique(wc)) # put in order
  
  # Get number of communities
  communities <- length(unique_communities)
  
  # Get return order of node names and communities (without NA)
  return_node_order <- order(node_names)
  return_community_order <- order(
    unique_communities[unique_communities != "NA"]
  )
  
  # If all singleton communities, then send NA for all
  if(nodes == communities){
    
    # Send results
    return(
      list(
        unstd = unstandardized[return_node_order, return_community_order],
        std = unstandardized[return_node_order, return_community_order]
      )
    )
    
  }
  
  # Not all singleton dimensions, so carry on
  
  # Check for method
  if(loading.method == "brm"){
    
    # Compute unstandardized loadings (absolute sums)
    unstandardized <- absolute_weights(A, wc, nodes, unique_communities)
    
    # Add signs to the loadings
    unstandardized <- old_add_signs(unstandardized, A, wc, unique_communities)
    
    # Before rounding occurred prior to standardization (no rounding)
    standardized <- standardize(unstandardized)
    
    # Get descending order
    standardized <- descending_order(standardized, wc, unique_communities)
    
    # Set up results
    results <- list(
      unstd = unstandardized[dimnames(standardized)[[1]],],
      std = standardized
    )
    
    # Add attributes
    attr(results, "methods") <- list(
      loading.method = loading.method, rotation = rotation
    )
    
    # Add class
    class(results) <- "net.loads"
    
    # Return results
    return(results)
    
  }
  
  
  # If not "BRM", run experimental
  
  # Experimental unstandardized loadings
  unstandardized <- experimental(
    A, wc, nodes, node_names, communities, unique_communities
  )
  
  # Obtain standardized loadings
  standardized <- standardize(unstandardized)
  
  # Get descending order
  standardized <- descending_order(standardized, wc, unique_communities)
    
  # Check for rotation
  if(!is.null(rotation)){
    
    # Errors for...
    # Missing packages: {GPArotation} and {fungible}
    # Invalid rotations
    rotation_errors(rotation)
    
    # If rotation exists, then obtain it
    rotation_FUN <- get(rotation, envir = asNamespace("GPArotation"))
    
    # Get ellipse arguments
    ellipse <- list(...)
    
    # Get arguments for function
    rotation_ARGS <- obtain_arguments(rotation_FUN, ellipse)
    
    # Check for "NA" community
    if("NA" %in% wc){
      standardized <- standardized[, dimnames(standardized)[[2]] != "NA"]
      communities <- communities - 1
      unique_communities <- unique_communities[unique_communities != "NA"]
    }
    
    # Supply loadings
    rotation_ARGS$A <- standardized
    
    # Set default arguments for rotations
    rotation_ARGS <- rotation_defaults(rotation, rotation_ARGS, ellipse)
    
    # Perform rotations
    rotation_OUTPUT <- do.call(rotation_FUN, rotation_ARGS)
    
    # Align rotated loadings
    aligned_output <- fungible::faAlign(
      F1 = standardized,
      F2 = rotation_OUTPUT$loadings,
      Phi2 = rotation_OUTPUT$Phi
    )
    
    # Set rotated loadings objects
    ## Loadings
    rotated_loadings <- aligned_output$F2
    dimnames(rotated_loadings) <- dimnames(standardized)
    ## Phi
    rotated_Phi <- aligned_output$Phi2
    dimnames(rotated_Phi) <- list(unique_communities, unique_communities)
    
    # Make rotated results list
    rotated <- list(
      loadings = rotated_loadings,
      Phi = rotated_Phi
    )
    
  }else{ # If rotation is NULL, then rotated is NULL
    rotated <- NULL
  }
  
  # Set up results
  results <- list(
    unstd = unstandardized[dimnames(standardized)[[1]],],
    std = standardized,
    rotated = rotated
  )
  
  # Add "methods" attributes
  attr(results, "methods") <- list(
    loading.method = loading.method, rotation = rotation
  )
  
  # Set class
  class(results) <- "net.loads"
  
  # Return results
  return(results)
  
  
}


# Bug checking ----
# 
# set.seed(1234)
# 
# # Generate data
# sim_data <- latentFactoR::simulate_factors(
#   factors = 3,
#   variables = 10,
#   loadings = 0.60,
#   cross_loadings = 0.10,
#   correlations = 0.30,
#   sample_size = 1000,
#   variable_categories = 5,
#   skew_range = c(-1, 1)
# )
# 
# # Add wording effects (for negative loadings)
# sim_data <- latentFactoR::add_wording_effects(
#   sim_data, method = "mixed"
# )
# 
# # Estimate EGA
# ega <- EGA(sim_data$data, plot.EGA = FALSE)
# ega$wc[8] <- NA
# A = ega; loading.method = "brm"
# rotation = "geominq"

#' @noRd
# Organize input ----
# Updated 10.07.2023
organize_input <- function(A, wc)
{
  
  # Check for `EGA` object
  if(any(class(A) %in% c("EGA", "EGA.fit", "riEGA"))){
    
    # Get `EGA` object
    ega_object <- get_EGA_object(A)
    
    # Set network and memberships
    A <- ega_object$network
    wc <- ega_object$wc
    
  }else{
    
    # Produce errors for miss aligned data
    length_error(wc, dim(A)[2]) # length between network and memberships
    object_error(A, c("matrix", "data.frame")) # must be matrix or data frame
    object_error(wc, c("vector", "matrix", "data.frame")) # must be one of these
    
  }
  
  # Generally, good to proceed
  A <- as.matrix(A); wc <- force_vector(wc)
  
  # Set memberships as string
  wc <- paste(wc)
  
  # Ensure names
  A <- ensure_dimension_names(A)
  names(wc) <- dimnames(A)[[2]]
  
  # Set orders
  ordering <- order(wc)
  
  # Return ordered network and memberships
  return(
    list(A = A[ordering, ordering], wc = wc[ordering])
  )
  
}

#' @noRd
# Obtain signs ----
# Function to obtain signs on dominant community
# Updated 11.07.2023
obtain_signs <- function(target_network)
{
  
  # Initialize signs to all positive orientation
  signs <- rep(1, dim(target_network)[2])
  names(signs) <- dimnames(target_network)[[2]]
  
  # Initialize row sums and minimum value
  row_sums <- rowSums(target_network, na.rm = TRUE)
  minimum_value <- which.min(row_sums)
  
  # Set while loop
  while(sign(row_sums[minimum_value]) == -1){
    
    # Get minimum value name
    minimum_name <- names(minimum_value)
    
    # Flip variable
    target_network[minimum_name,] <- 
      target_network[,minimum_name] <-
      -target_network[,minimum_name]
    
    # Set sign as flipped
    signs[minimum_name] <- -signs[minimum_name]
    
    # Update row sums and minimum value
    row_sums <- rowSums(target_network, na.rm = TRUE)
    minimum_value <- which.min(row_sums)
    
  }
  
  # Add signs as an attribute to the target network
  attr(target_network, "signs") <- signs
  
  # Return results
  return(target_network)
  
}

#' @noRd
# Experimental loadings ----
# Updated 12.07.2023
experimental <- function(A, wc, nodes, node_names, communities, unique_communities)
{
  
  # Initialize loading matrix
  loading_matrix <- matrix(
    0, nrow = nodes, ncol = communities,
    dimnames = list(node_names, unique_communities)
  )
  
  # Initialize sign vector
  signs <- rep(1, nodes)
  names(signs) <- node_names
  
  # Populate loading matrix
  for(community in unique_communities){
    
    # Get community index
    community_index <- wc == community
    
    # Determine positive direction for dominant loadings
    target_network <- obtain_signs(A[community_index, community_index, drop = FALSE])
    
    # Compute absolute sum for dominant loadings
    loading_matrix[community_index, community] <- colSums(target_network, na.rm = TRUE)
    
    # Determine positive direction for dominant loadings
    signs[community_index] <- attr(target_network, "signs")
    
  }
  
  # Check for unidimensional structure
  if(communities > 1){
    
    # Check for any negative signs
    if(any(signs == -1)){
      
      # Make a copy of the network
      A_copy <- A
      
      # Flip signs
      A[signs == -1,] <- A[,signs == -1] <- A_copy[signs == -1,]
      
      
    }
    
    # Populate loading matrix with cross-loadings
    for(community in unique_communities){
      for(cross in unique_communities){
        
        # No need for same community loadings
        if(community != cross){
          
          # Get community index
          community_index <- wc == community
          
          # Compute algebraic sum for cross-loadings
          loading_matrix[community_index, cross] <- colSums(
            A[wc == cross, community_index, drop = FALSE], na.rm = TRUE
          )
          
        }
        
      }
    }
    
  }
  
  # Set signs
  return(loading_matrix * signs)
  
}

#' @noRd
# Standardize loadings ----
# Updated 12.07.2023
standardize <- function(unstandardized)
{
  return(t(t(unstandardized) / sqrt(colSums(abs(unstandardized), na.rm = TRUE))))
}

#' @noRd
# Descending order ----
# Updated 11.07.2023
descending_order <- function(standardized, wc, unique_communities) 
{
  
  # Initialize order names
  order_names <- character(dim(standardized)[1])
  
  # Loop over communities
  for(community in unique_communities){
    
    # Get community index
    community_index <- wc == community
    
    # Get order
    ordering <- order(standardized[community_index, community])
    
    # Input ordering into order names
    order_names[community_index] <- dimnames(standardized)[[1]][community_index]
    
  }
  
  
  
  # Return reordered results
  return(standardized[order_names,])
  
}

#' @noRd
# Rotation errors ----
# Updated 12.07.2023
rotation_errors <- function(rotation)
{
  
  # Check for packages
  ## Needs {GPArotation} and {fungible}
  check_package(c("GPArotation", "fungible"))
  
  # Get rotations available in {GPArotation}
  rotation_names <- ls(asNamespace("GPArotation"))
  
  # Check if rotation exists
  if(!rotation %in% rotation_names){
    
    # Send error that rotation is not found
    stop(
      paste0(
        "Invalid rotation: ", rotation, "\n\n",
        "The rotation \"", rotation, "\" is not available in the {GPArotation} package. ",
        "\n\nSee `?GPArotation::rotations` for the list of available rotations."
      )
    )
    
  }
  
}

#' @noRd
# Rotation default arguments ----
# Updated 12.07.2023
rotation_defaults <- function(rotation, rotation_ARGS, ellipse)
{
  
  # Check for "n.rotations" (used in {psych})
  if("n.rotations" %in% ellipse){
    rotation_ARGS$randomStarts <- ellipse$n.rotations
  }
  
  # Check for random starts
  if(!"randomStarts" %in% names(ellipse) & !"n.rotations" %in% names(ellipse)){
    rotation_ARGS$randomStarts <- 10
  }
  
  # Check for maximum iterations argument
  if(!"maxit" %in% names(ellipse)){
    rotation_ARGS$maxit <- 1000
  }
  
  # Check for epsilon argument
  if(!"eps" %in% names(ellipse) & grepl("geomin", rotation)){
    
    # Based on number of dimensions, switch epsilon
    rotation_ARGS$eps <- switch(
      as.character(dim(rotation_ARGS$A)[2]),
      "2" = 0.0001, # two dimensions
      "3" = 0.001, # three dimensions
      0.01 # four or more dimensions
    )
    
  }
  
  # Return arguments
  return(rotation_ARGS)
  
}

#%%%%%%%%%%%%%%%%%
# BRM Legacy ----
#%%%%%%%%%%%%%%%%%

#' @noRd
## Absolute weights ("BRM") ----
# Updated 10.07.2023
absolute_weights <- function(A, wc, nodes, unique_communities)
{
  
  # Ensure network is absolute
  A <- abs(A)
  
  # Loop over communities
  return(
    nvapply(
      unique_communities, function(community){
        colSums(A[wc == community,, drop = FALSE], na.rm = TRUE)
      }, LENGTH = nodes
    )
  )
  
}

#' @noRd
## Add signs ("BRM") ----
# From CRAN version 1.2.3
# Updated 10.07.2023
old_add_signs <- function(unstandardized, A, wc, unique_communities)
{
  
  # Loop over main loadings
  for(community in unique_communities){
    
    # Get community index
    community_index <- wc == community
    
    # Get number of nodes
    node_count <- sum(community_index)
    
    # Get community sub-network
    community_network <- A[community_index, community_index, drop = FALSE]
    
    # Initialize sign matrix
    community_signs <- sign(community_network)
    
    # Initialize signs to all positive
    signs <- rep(1, node_count)
    
    # Loop over nodes
    for(node in seq_len(node_count)){
      
      # Make copy of signs
      signs_copy <- community_signs
      
      # Get current maximum sum
      current_max <- sum(colSums(community_signs, na.rm = TRUE), na.rm = TRUE)
      
      # Flip sign of each node
      community_signs[node,] <- -community_signs[node,]
      
      # Get new maximum sum
      new_max <- sum(colSums(community_signs, na.rm = TRUE), na.rm = TRUE)
      
      # Check for increase
      if(new_max > current_max){
        signs[node] <- -1 # with increase, flip sign
      }else{ # otherwise, return sign matrix to original state
        community_signs <- signs_copy
      }
      
    }
    
    # Update signs in loadings
    unstandardized[community_index, community] <-
      unstandardized[community_index, community] * signs
    
    # Sweep across community
    A[, community_index] <- sweep(
      A[, community_index, drop = FALSE], MARGIN = 2, signs, `*`
    )
    
  }
  
  # Loop over communities
  for(community1 in unique_communities){
    
    # Get first community index
    community_index1 <- wc == community1
    
    # Get number of nodes
    node_count <- sum(community_index1)
    
    # Loop over other communities
    for(community2 in unique_communities){
      
      # Check for the same community
      if(community1 != community2){
        
        # Get second community index
        community_index2 <- wc == community2
        
        # Get community sub-network
        community_network <- A[community_index1, community_index2, drop = FALSE]
        
        # Initialize sign matrix
        community_signs <- sign(community_network)
        
        # Initialize signs to all positive
        signs <- rep(1, node_count)
        
        # Loop over nodes
        for(node in seq_len(node_count)){
          
          # Make copy of signs
          signs_copy <- community_signs
          
          # Get current maximum sum
          current_max <- sum(colSums(community_signs, na.rm = TRUE), na.rm = TRUE)
          
          # Flip sign of each node
          community_signs[node,] <- -community_signs[node,]
          
          # Get new maximum sum
          new_max <- sum(colSums(community_signs, na.rm = TRUE), na.rm = TRUE)
          
          # Check for increase
          if(new_max > current_max){
            signs[node] <- -1 # with increase, flip sign
          }else{ # otherwise, return sign matrix to original state
            community_signs <- signs_copy
          }
          
        }
        
        # Update signs in loadings
        unstandardized[community_index1, community2] <-
          unstandardized[community_index1, community2] * signs
        
      }
      
    }
    
  }
  
  # Flip direction of community with main loadings
  for(community in unique_communities){
    
    # Get community indices
    community_index <- wc == community
    
    # Determine direction with sign
    if(sign(sum(unstandardized[community_index, community])) != 1){
      unstandardized[community_index,] <- 
        -unstandardized[community_index,]
    }
    
  }
  
  # Return unstandardized loadings
  return(unstandardized)
  
}

# SAVE ----

#' Oblimin rotation (to be consistent in \code{\link[EGAnet]{net.loads}})
#' @noRd
# Updated 28.11.2022 -- Marcos
oblimin_rotate <- function(
    loadings,
    n.rotations = 10, # defaults to 10
    nfactors, # number of factors
    maxit = 1e4, # iterations
    eps = 1e-5, # convergence
    gamma = 0, # gamma in oblimin
    rotate = "oblimin" # rotation
)
{
  
  # Check for number of factors
  if(missing(nfactors)){
    nfactors <- ncol(loadings)
  }
  
  # Check for unidimensional structure (28.11.2022)
  if(nfactors == 1){
    
    # Return results
    return(
      list(
        loadings = loadings,
        Phi = 1
      )
    )
    
  }
  
  x <- list()
  fs <- vector(length = n.rotations)
  
  for(i in 1:n.rotations) {
    
    X <- replicate(nfactors, rnorm(nfactors))
    Q <- qr.Q(qr(X)) # Random orthogonal matrix (initial value)
    x[[i]] <- GPArotation::GPFoblq(A = loadings, method = rotate, 
                                   Tmat = Q, maxit = maxit, eps = eps) # Fit
    fs[i] <- oblimin(x[[i]]$loadings, gamma) # Fit values
    
  }
  
  # Results
  index <- which.min(fs) # Index pertaining to the minimum fit value
  loadings <- x[[i]]$loadings # Select the corresponding loadings
  Phi <- x[[i]]$Phi # Factor correlation matrix
  
  # Return results
  return(
    list(
      loadings = loadings,
      Phi = Phi
    )
  )
  
}

#' Oblimin rotation
#' @noRd
# Updated 21.11.2022 -- Marcos
oblimin <- function(L, gamma = 0) {
  
  # Fit value for oblimin
  # L = rotated loading matrix
  # gamma = fixed parameter
  
  nfactors <- ncol(L)
  p <- nrow(L)
  I <- diag(p)
  gC <- matrix(gamma/p, p, p)
  IgC <- I - gC
  
  N <- matrix(1, nfactors, nfactors)
  diag(N) <- 0
  L2 <- L*L
  f <- sum(diag(t(L2) %*% (IgC %*% L2 %*% N)))/4
  return(f)
  
}

#' Factor scores
#' @noRd
# Updated 28.11.2022 -- Marcos
fscores <- function(S, loadings, Phi, Shat, scores = NULL, method = "Thurstone") {
  
  # S = empirical correlation matrix
  # loadings = rotated loading matrix
  # Phi = estimated factor correlation matrix
  # Shat = model correlation matrix
  # scores = raw scores
  # method = factor.scores method
  
  if(is.null(scores)) stop("Please, provide the matrix of observed scores")
  
  uniquenesses <- 1 - diag(Shat)
  invS <- solve(S)
  n <- nrow(scores)
  p <- nrow(loadings)
  q <- ncol(loadings)
  z <- scale(scores)
  
  Lambda <- loadings
  Phi <- Phi
  LP <- Lambda %*% Phi # Correlations between factors and items
  
  # Find the weights:
  
  if(method == "regression" | method == "Thurstone") {
    
    weights <- solve(S, LP)
    
  } else if(method == "tenBerge") {
    
    SVD <- svd(Phi)
    Phi12 <- SVD$u %*% diag(sqrt(SVD$d)) %*% t(SVD$v)
    SVD <- svd(S)
    R12 <- SVD$u %*% diag(1/sqrt(SVD$d)) %*% t(SVD$v)
    L <- Lambda %*% Phi12
    SVD <- svd(t(L) %*% invS %*% L)
    LRL12 <- SVD$u %*% diag(1/sqrt(SVD$d)) %*% t(SVD$v)
    C <- R12 %*% L %*% LRL12
    weights <- R12 %*% C %*% Phi12
    
  } else if(method == "Bartlett") {
    
    U <- c(fit$efa$uniquenesses)
    U2 <- diag(1/(U*U))
    weights <- U2 %*% Lambda %*% solve(t(Lambda) %*% U2 %*% Lambda)
    
  } else if(method == "Harman") {
    
    weights <- solve(fit$efa$Rhat) %*% Lambda
    
  }
  
  fs <- z %*% weights # Factor scores
  
  # Validity coefficients:
  # invL <- diag(1/apply(fs, MARGIN = 2, FUN = sd))
  C <- t(weights) %*% S %*% weights
  # invL was updated (28.11.2022) to ensure matrix
  # when nfactors = 1
  invL <- diag(sqrt(diag(C)), ncol(loadings)) # Standard deviations of the factor scores
  validity_univocality <- t(LP) %*% weights %*% invL
  validity <- matrix(diag(validity_univocality), nrow = 1)
  rownames(validity) <- ""
  univocality <- validity_univocality
  diag(univocality) <- NA
  
  # Accuracy:
  accuracy <- stats::cor(fs)
  
  # Standard errors for factor scores:
  r <- matrix(diag(invS), ncol = 1)
  Rj <- matrix(1-c(validity^2), nrow = 1)
  se <- sqrt(r %*% Rj / (n-p-1))
  se <- matrix(se, nrow = p, ncol = q)
  
  colnames(fs) <- colnames(weights) <- colnames(validity) <-
    colnames(univocality) <- rownames(univocality) <-
    colnames(accuracy) <- rownames(accuracy) <- colnames(se) <-
    paste("F", sprintf(paste("%0", nchar(q), "d", sep = ""), 1:q), sep = "")
  
  result <- list(fscores = fs, weights = weights, validity = validity,
                 univocality = univocality, accuracy = accuracy, se = se)
  
  return(result)
  
}

#' Exploratory factor analysis
#' @noRd
# Updated 21.11.2022 -- Marcos
efa <- function(data, nfactors, fm = "minres", rotate = "oblimin",
                n.rotations = 10, maxit = 1e4, factor.scores = "Thurstone",
                gamma = 0, eps = 1e-5) {
  
  # Function to perform EFA
  # data = raw data
  # nfactors = number of factors to extract
  # fm = factor extraction method
  # rotate = rotation method (only oblimin is implemented right now)
  # n.rotations = number of rotations to perform with different initial values
  # maxit = maximum number of iterations for the rotation to converge
  # factor.scores = factor scores' method
  # gamma = fixed parameter for oblimin
  # eps = stop criteria for the rotation (relative tolerance value)
  
  if(rotate != "oblimin") stop("oblimin is the only rotation currently implemented.")
  
  # Compute the correlation matrix
  S <- qgraph::cor_auto(data, verbose = FALSE)
  
  # Factor extraction:
  fa <- psych::fa(S, nfactors = nfactors, fm = fm, rotate = "none")

  # Perform oblimin rotation
  rotation <- oblimin_rotate(
    loadings = fa$loadings,
    n.rotations = n.rotations,
    nfactors = nfactors,
    maxit = maxit,
    eps = eps,
    gamma = gamma,
    rotate = rotate
  )
  loadings <- rotation$loadings # Select the corresponding loadings
  Phi <- rotation$Phi # Factor correlation matrix
  Shat <- fa$model # Model correlation matrix
  uniquenesses <- fa$uniquenesses
  
  # Compute the factor scores
  fs <- fscores(S = S, loadings = loadings, Phi = Phi, 
                Shat = Shat, scores = data, method = factor.scores)
  
  result <- list(loadings = loadings, Phi = Phi, Shat = Shat,
                 uniquenesses = uniquenesses, factor.scores = fs$fscores)
  
}
