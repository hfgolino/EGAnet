#' @title Network Loadings
#'
#' @description Computes the between- and within-community
#' \code{strength} of each variable for each community
#'
#' @param A Network matrix, data frame, or \code{\link[EGAnet]{EGA}} object
#'
#' @param wc Numeric or character vector (length = \code{ncol(A)}).
#' A vector of community assignments.
#' If input into \code{A} is an \code{\link[EGAnet]{EGA}} object,
#' then \code{wc} is automatically detected
#'
#' @param loading.method Character (length = 1).
#' Sets network loading calculation based on implementation
#' described in \code{"original"} (Christensen & Golino, 2021) or
#' the \code{"revised"} (Christensen et al., 2024) implementation.
#' Defaults to \code{"revised"}
#'
#' @param scaling Numeric (length = 1).
#' Scaling factor for the magnitude of the \code{"experimental"} network loadings.
#' Defaults to \code{2}.
#' \code{10} makes loadings roughly the size of factor loadings when correlations
#' between factors are orthogonal
#'
#' @param rotation Character.
#' A rotation to use to obtain a simpler structure.
#' For a list of rotations, see \code{\link[GPArotation]{rotations}} for options.
#' Defaults to \code{NULL} or no rotation.
#' By setting a rotation, \code{scores} estimation will be
#' based on the rotated loadings rather than unrotated loadings
#'
#' @param ... Additional arguments to pass on to \code{\link[GPArotation]{rotations}}
#'
#' @return Returns a list containing:
#'
#' \item{unstd}{A matrix of the unstandardized within- and between-community
#' strength values for each node}
#'
#' \item{std}{A matrix of the standardized within- and between-community
#' strength values for each node}
#'
#' \item{rotated}{\code{NULL} if \code{rotation = NULL}; otherwise,
#' a list containing the rotated standardized network loadings
#' (\code{loadings}) and correlations between dimensions (\code{Phi})
#' from the rotation}
#'
#' @details Simulation studies have demonstrated that a node's strength
#' centrality is roughly equivalent to factor loadings
#' (Christensen & Golino, 2021; Hallquist, Wright, & Molenaar, 2019).
#' Hallquist and colleagues (2019) found that node strength represented a
#' combination of dominant and cross-factor loadings. This function computes
#' each node's strength within each specified dimension, providing a rough
#' equivalent to factor loadings (including cross-loadings; Christensen & Golino, 2021).
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#'
#' # Network loadings
#' net.loads(ega.wmt)
#'
#' @references
#' \strong{Original implementation and simulation} \cr
#' Christensen, A. P., & Golino, H. (2021).
#' On the equivalency of factor and network loadings.
#' \emph{Behavior Research Methods}, \emph{53}, 1563-1580.
#'
#' \strong{Demonstration of node strength similarity to CFA loadings} \cr
#' Hallquist, M., Wright, A. C. G., & Molenaar, P. C. M. (2019).
#' Problems with centrality measures in psychopathology symptom networks: Why network psychometrics cannot escape psychometric theory.
#' \emph{Multivariate Behavioral Research}, 1-25.
#'
#' \strong{Revised network loadings} \cr
#' Christensen, A. P., Golino, H., Abad, F. J., & Garrido, L. E. (2024).
#' Revised network loadings.
#' \emph{PsyArXiv}.
#'
#' @author Alexander P. Christensen <alexpaulchristensen@gmail.com> and Hudson Golino <hfg9s at virginia.edu>
#'
#' @export
#'
# Network Loadings ----
# Updated 06.03.2025
net.loads <- function(
    A, wc, loading.method = c("original", "revised"),
    scaling = 2, rotation = NULL, ...
)
{

  # Check for no input in 'loading.method'
  if(length(loading.method) > 1){
    loading.method <- "revised"
  }else{

    # Switch out old calls
    loading.method <- switch(
      tolower(loading.method),
      "brm" = "original",
      "experimental" = "revised",
      loading.method
    )

    # Check for missing arguments (argument, default, function)
    loading.method <- set_default(loading.method, "revised", net.loads)

  }

  # Organize and extract input (handles argument errors)
  # `wc` is made to be a character vector to allow `NA`
  input <- organize_input(A, wc)
  A <- input$A; wc <- input$wc

  # Get number of nodes and their names
  nodes <- length(wc); node_names <- names(wc)

  # Get unique communities (`NA` is OK)
  unique_communities <- sort(unique(wc)) # put in order

  # Get number of communities
  communities <- length(unique_communities)

  # If all singleton communities, then send NA for all
  if(nodes == communities){

    # Initialize loading matrix
    loading_matrix <- matrix(
      NA, nrow = nodes, ncol = communities,
      dimnames = list(node_names, unique_communities)
    )

    # Set up results
    results <- list(
      unstd = loading_matrix,
      std = loading_matrix
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

  # Not all singleton dimensions, so carry on

  # Check for method
  if(loading.method == "revised"){

    # Revised unstandardized loadings
    unstandardized <- revised_loadings(
      A, wc, nodes, node_names, communities, unique_communities, flip
    )

    # Store attributes
    unstd_attributes <- attr(unstandardized, "community")

  }else{

    # Compute unstandardized loadings (absolute sums)
    unstandardized <- absolute_weights(A, wc, nodes, unique_communities)

    # Add signs to the loadings
    unstandardized <- old_add_signs(unstandardized, A, wc, unique_communities)


  }

  # Obtain standardized loadings
  standardized <- standardize(unstandardized, loading.method, A, wc, scaling)

  # Get descending order
  standardized <- descending_order(standardized, wc, unique_communities, node_names)

  # Check for rotation
  if(!is.null(rotation)){

    # Errors for...
    # Missing packages: {GPArotation} and {fungible}
    # Invalid rotations
    # Returns: proper capitalization of rotation
    # For example: "geominq" returns "geominQ"
    rotation <- rotation_errors(rotation)

    # If rotation exists, then obtain it
    rotation_FUN <- get(rotation, envir = asNamespace("GPArotation"))

    # Get ellipse arguments
    ellipse <- list(...)

    # Get arguments for function
    rotation_ARGS <- obtain_arguments(rotation_FUN, ellipse)

    # Check for "NA" community
    NA_community <- unique_communities == "NA"

    # Check for "NA" community
    if(any(NA_community)){
      unique_communities <- unique_communities[!NA_community]
      standardized <- standardized[,unique_communities]
      communities <- communities - 1
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

  # Add "community" attributes
  if(loading.method == "revised"){
    attr(results, "methods")$scaling <- scaling
    attr(results, "community") <- unstd_attributes
  }

  # Add "membership" attributes for `net.scores`
  attr(results, "membership") <- list(
    wc = wc
  )

  # Set class
  class(results) <- "net.loads"

  # Send message about changed defaults
  if(loading.method == "revised"){
    message(
      paste(
        "The default 'loading.method' has changed to \"revised\" in {EGAnet} version >= 2.0.7.\n\n",
        "For the previous default (version <= 2.0.6), use `loading.method = \"original\"`"
      )
    )
  }

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
# A = ega; loading.method = "revised"
# rotation = "geominq"

#' @exportS3Method
# S3 Print Method ----
# Updated 12.08.2024
print.net.loads <- function(x, ...)
{

  # Get ellipse arguments
  ellipse <- list(...)

  # Get method attributes
  method_attributes <- attr(x, "methods")

  # Print method
  cat(
    paste0(
      "Loading Method: ", swiftelse(
        method_attributes$loading.method == "revised",
        "Revised", "Original"
      )
    )
  )

  # Check for rotation
  if(!is.null(method_attributes$rotation)){

    # Print rotation
    cat(
      paste0(
        "\nRotation: ", method_attributes$rotation
      )
    )

  }

  # Add breakspace
  cat("\n\n")

  # Get rounded loadings
  rounded_loadings <- round(x$std, 3)

  # Get minimum loadings value
  minimum <- swiftelse(
    "minimum" %in% names(ellipse),
    ellipse$minimum, 0.10
  )

  # Set loadings below minimum to empty string
  rounded_loadings[abs(x$std) < minimum] <- ""

  # Set loadings
  print(
    column_apply(
      as.data.frame(rounded_loadings),
      format_decimal, places = 3
    ), quote = FALSE
  )

  # Add message about minimum loadings
  cat(
    paste0(
      "Standardized loadings >= |", format_decimal(minimum, 2),
      "| are displayed. To change this 'minimum', use ",
      "`print(net.loads_object, minimum = 0.10)`"
    )
  )

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 12.07.2023
summary.net.loads <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @noRd
# Organize input ----
# Updated 13.08.2023
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
    length_error(wc, dim(A)[2], "net.loads") # length between network and memberships
    object_error(A, c("matrix", "data.frame"), "net.loads") # must be matrix or data frame
    object_error(wc, c("vector", "matrix", "data.frame"), "net.loads") # must be one of these

  }

  # Generally, good to proceed
  A <- as.matrix(A); wc <- force_vector(wc)

  # Set memberships as string
  if(is.numeric(wc)){
    wc <- format_integer(
      wc, places = digits(max(wc, na.rm = TRUE)) - 1
    )
  }

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

# Faster sign method (faster than eigenvectors but only marginally)
# Obtain signs ----
# Function to obtain signs on dominant community
# Updated 22.03.2024
obtain_signs <- function(target_network)
{

  # Initialize signs to all positive orientation
  signs <- rep(1, dim(target_network)[2])
  names(signs) <- dimnames(target_network)[[2]]

  # Initialize row sums and minimum index
  row_sums <- rowSums(target_network, na.rm = TRUE)
  minimum_index <- which.min(row_sums)

  # Set while loop
  while(sign(row_sums[minimum_index]) == -1){

    # Flip variable
    target_network[minimum_index,] <-
      target_network[,minimum_index] <-
        -target_network[minimum_index,]

    # Set sign as flipped
    signs[minimum_index] <- -signs[minimum_index]

    # Update row sums and minimum value
    row_sums <- rowSums(target_network, na.rm = TRUE)
    minimum_index <- which.min(row_sums)

  }

  # Determine whether signs should be flipped
  if(sum(signs) < 0){
    signs <- -signs
  }

  # Add signs as an attribute to the target network
  attr(target_network, "signs") <- signs

  # Return results
  return(target_network)

}

#' @noRd
# Revised loadings ----
# Updated 22.08.2024
revised_loadings <- function(
    A, wc, nodes, node_names,
    communities, unique_communities
)
{

  # Initialize loading matrix
  loading_matrix <- matrix(
    0, nrow = nodes, ncol = communities,
    dimnames = list(node_names, unique_communities)
  )

  # Initialize sign vector
  signs <- rep(1, nodes)
  names(signs) <- node_names

  # Get community numbers
  community_table <- fast_table(wc)

  # Populate loading matrix
  for(community in unique_communities){

    # Get community index
    community_index <- wc == community

    # Obtain target network
    target_network <- A[community_index, community_index, drop = FALSE]

    # Determine positive direction for dominant loadings
    target_network <- obtain_signs(
      A[community_index, community_index, drop = FALSE]
    )

    # Compute absolute sum for dominant loadings
    loading_matrix[community_index, community] <- colSums(
      target_network, na.rm = TRUE
    ) / (community_table[community] - 1)

    # Determine positive direction for dominant loadings
    signs[community_index] <- attr(target_network, "signs")


    # Revert back to original sign algorithm
    # Eigenvectors depend on the matrix manipulation to
    # orient variables in the proper direction to get the
    # appropriate signs

    # # Compute absolute sum for dominant loadings
    # loading_matrix[community_index, community] <- colSums(
    #   obtain_signs(target_network), na.rm = TRUE
    # ) / (community_table[community] - 1)
    #
    # # Obtain signs
    # target_signs <- sign(eigen(target_network, symmetric = TRUE)$vector[,1])
    # # Thank you to Sacha Epskamp for pointing out this simpler approach to us!
    #
    # # Determine positive direction for dominant loadings
    # signs[community_index] <- swiftelse(
    #   sum(target_signs) < 0, -target_signs, target_signs
    # )

  }

  # Take the average of the within-community values
  # and multiply them by the number of values
  loading_matrix <- sweep(
    loading_matrix, 2,
    STATS = community_table,
    FUN = "*"
  )

  # Compute sums
  community_sums <- colSums(abs(loading_matrix), na.rm = TRUE)

  # Check for unidimensional structure
  if(communities > 1){

    # Get negative sign indices
    negative_signs <- which(signs == -1)

    # Loop over negative signs
    if(length(negative_signs) != 0){

      # Flip signs
      for(negative in negative_signs){
        A[negative,] <- A[,negative] <- -A[,negative]
      }

    }

    # Populate loading matrix with cross-loadings
    for(community in unique_communities){

      # Get community index
      community_index <- wc == community

      # Loop across other communities
      for(cross in unique_communities){

        # No need for same community loadings
        if(community != cross){

          # Compute algebraic sum for cross-loadings
          loading_matrix[community_index, cross] <- colSums(
            A[wc == cross, community_index, drop = FALSE], na.rm = TRUE
          )

        }

      }
    }

  }

  # Set signs
  loading_matrix <- loading_matrix * signs

  # Add attributes
  attr(loading_matrix, "community") <- list(
    community_sums = community_sums,
    community_table = community_table
  )

  # Return loading matrix
  return(loading_matrix)

}

#' @noRd
# Standardize loadings ----
# Updated 12.08.2024
standardize <- function(unstandardized, loading.method, A, wc, scaling)
{

  # Check for loading method
  if(loading.method == "revised"){

    # Get attributes
    community <- attr(unstandardized, "community")

    # Return loadings
    return(
      t(t(unstandardized) / (community$community_sums^(1 / log(scaling * community$community_table))))
    )

  }else if(loading.method == "original"){
    return(t(t(unstandardized) / sqrt(colSums(abs(unstandardized), na.rm = TRUE))))
  }

}

#' @noRd
# Convert scale of revised loadings ----
# Updated 26.09.2024
scaling_conversion <- function(standardized, community, original_scaling, new_scaling)
{

  # Needs 'community' attribute from the unstandardized loadings
  # This attribute is attached to the overall loading return when 'loading.method = "revised"'

  # Compute change in scaling
  delta <- community$community_sums^(1 / log(original_scaling * community$community_table)) /
    community$community_sums^(1 / log(new_scaling * community$community_table))

  # Return adjustment
  return(t(t(standardized) * delta))

}

#' @noRd
# Descending order ----
# Updated 24.07.2023
descending_order <- function(standardized, wc, unique_communities, node_names)
{

  # Initialize order names
  order_names <- character(length(node_names))

  # Get order names
  order_names <- ulapply(
    unique_communities, function(community){

      # Get community index
      community_index <- wc == community

      # Get order
      ordering <- order(standardized[community_index, community], decreasing = TRUE)

      # Return order
      return(node_names[community_index][ordering])

    }
  )

  # Return reordered results
  return(standardized[order_names,, drop = FALSE])

}

#' @noRd
# Rotation errors ----
# Updated 13.08.2023
rotation_errors <- function(rotation)
{

  # Check for packages
  ## Needs {GPArotation} and {fungible}
  check_package(c("GPArotation", "fungible"))

  # Get rotations available in {GPArotation}
  rotation_names <- ls(asNamespace("GPArotation"))

  # Get lowercase names
  rotation_lower <- tolower(rotation)
  rotation_names_lower <- tolower(rotation_names)

  # Check if rotation exists
  if(!rotation_lower %in% rotation_names_lower){

    # Send error that rotation is not found
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Invalid rotation: ", rotation, "\n\n",
        "The rotation \"", rotation, "\" is not available in the {GPArotation} package. ",
        "\n\nSee `?GPArotation::rotations` for the list of available rotations."
      ),
      call = "net.loads"
    )

  }

  # If rotation exists, then return proper name
  return(rotation_names[rotation_lower == rotation_names_lower])

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

#%%%%%%%%%%%%%%%%%%%%%
# Original Legacy ----
#%%%%%%%%%%%%%%%%%%%%%

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
# Updated 13.07.2023
old_add_signs <- function(unstandardized, A, wc, unique_communities)
{

  # Loop over main loadings
  for(community in unique_communities){

    # Get community index
    community_index <- wc == community

    # Get number of nodes
    node_count <- sum(community_index)

    # Initialize sign matrix
    community_signs <- sign(A[community_index, community_index, drop = FALSE])

    # Initialize signs to all positive
    signs <- rep(1, node_count)

    # Loop over nodes
    for(node in seq_len(node_count)){

      # Make copy of signs
      signs_copy <- community_signs

      # Get current maximum sum
      current_max <- sum(community_signs, na.rm = TRUE)

      # Flip sign of each node
      community_signs[node,] <- -community_signs[node,]

      # Get new maximum sum
      new_max <- sum(community_signs, na.rm = TRUE)

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

        # Initialize sign matrix
        community_signs <- sign(A[community_index1, community_index2, drop = FALSE])

        # Initialize signs to all positive
        signs <- rep(1, node_count)

        # Loop over nodes
        for(node in seq_len(node_count)){

          # Make copy of signs
          signs_copy <- community_signs

          # Get current maximum sum
          current_max <- sum(community_signs, na.rm = TRUE)

          # Flip sign of each node
          community_signs[node,] <- -community_signs[node,]

          # Get new maximum sum
          new_max <- sum(community_signs, na.rm = TRUE)

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

