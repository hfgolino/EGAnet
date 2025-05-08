#' @title Total Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit (TEFI) of a dimensionality structure using Von Neumman's entropy
#' when the input is a correlation matrix. Lower values suggest better fit of a structure to the data.
#'
#' @param data Matrix, data frame, or \code{*EGA} class object.
#' Matrix or data frame can be raw data or a correlation matrix.
#' All \code{*EGA} objects are accepted. \code{\link[EGAnet]{hierEGA}}
#' input will produced the Generalized TEFI (see \code{\link[EGAnet]{genTEFI}})
#'
#' @param structure Numeric or character vector (length = \code{ncol(data)}).
#' Can be theoretical factors or the structure detected by \code{\link{EGA}}
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{TRUE} to see all messages and warnings for every
#' function call.
#' Set to \code{FALSE} to ignore messages and warnings
#'
#' @return Returns a data frame with columns:
#'
#' \strong{Non-hierarchical Structure}
#'
#' \item{VN.Entropy.Fit}{The Total Entropy Fit Index using Von Neumman's entropy}
#'
#' \item{Total.Correlation}{The total correlation of the dataset}
#'
#' \item{Average.Entropy}{The average entropy of the dataset}
#'
#' \strong{Hierarchical Structure}
#'
#' \item{VN.Entropy.Fit}{The Generalized Total Entropy Fit Index using Von Neumman's entropy}
#'
#' \item{Level_#_VN}{An individual level's Von Neumann's entropy}
#'
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' # Estimate EGA model
#' ega.wmt <- EGA(
#'   data = wmt, model = "glasso",
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#'
#' # Compute entropy indices for empirical EGA
#' tefi(ega.wmt)
#'
#' # User-defined structure (with `EGA` object)
#' tefi(ega.wmt, structure = c(rep(1, 5), rep(2, 5), rep(3, 8)))
#'
#' @references
#' \strong{Initial formalization and simulation} \cr
#' Golino, H., Moulder, R. G., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Nesselroade, J., Sadana, R., Thiyagarajan, J. A., & Boker, S. M. (2020).
#' Entropy fit indices: New fit measures for assessing the structure and dimensionality of multiple latent variables.
#' \emph{Multivariate Behavioral Research}.
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen@gmail.com>, and Robert Moulder <rgm4fd@virginia.edu>
#'
#' @export
# Total Entropy Fit Index Function (for correlation matrices)
# Updated 09.08.2023
tefi <- function(data, structure = NULL, verbose = TRUE)
{

  # Check for errors (returns 'data' and 'ega_class' flag)
  error_return <- tefi_errors(data, verbose)

  # Get 'data' and 'ega_class' flag
  data <- error_return$data; ega_class <- error_return$ega_class

  # Branch for `EGA` class
  if(any(ega_class)){

    # Get `EGA` object
    ega_object <- get_EGA_object(data)

    # Get structure
    structure <- get_tefi_structure(data, structure, ega_object)

    # Get correlation matrix based on EGA
    if(is(data, "hierEGA")){
      correlation_matrix <- ega_object$lower_order$correlation
    }else{
      correlation_matrix <- ega_object$correlation
    }

  }else{ # Non-EGA objects

    # Get structure
    structure <- get_tefi_structure(data, structure, NULL)

    # Generic function to get necessary inputs
    output <- obtain_sample_correlations(
      data = data, n = 1L, # set to 1 to avoid error
      corr = "auto", na.data = "pairwise",
      verbose = verbose
    )

    # Get correlation matrix
    correlation_matrix <- output$correlation_matrix

  }

  # Get absolute correlation matrix
  correlation_matrix <- abs(correlation_matrix)

  # Branch based on hierarchical structure
  return(
    swiftelse( # hierarchical will be a list
      get_object_type(structure) == "list",
      generalized_tefi(correlation_matrix, structure, verbose),
      tefi_standard(correlation_matrix, structure, verbose)
    )
  )

}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]; ega.wmt <- EGA(data, plot.EGA = FALSE)
# data <- ega.wmt$correlation
# structure <- ega.wmt$wc

#' @noRd
# Argument errors
# Updated 03.05.2025
tefi_errors <- function(data, verbose)
{

  # Get `EGA` class
  ega_class <- grepl("EGA", class(data))

  # 'verbose' errors
  length_error(verbose, 1, "tefi")
  typeof_error(verbose, "logical", "tefi")

  # 'data' errors
  if(any(!ega_class)){

    # Check for appropriate data
    object_error(data, c("matrix", "data.frame", "tibble"), "tefi")

    # Check for tibble
    if(get_object_type(data) == "tibble"){
      data <- as.data.frame(data)
    }

    # Ensure usable data
    data <- usable_data(data, verbose)

  }

  # Return data and `EGA` classes
  return(list(data = data, ega_class = ega_class))

}

#' @noRd
# Handle structure input ----
get_tefi_structure <- function(data, structure, ega_object = NULL)
{

  # Check for whether `EGA` object is NULL
  if(is.null(ega_object)){

    # Get number of variables
    variables <- dim(data)[2L]

    # Determine if NULL
    if(is.null(structure)){

      stop(
        paste(
          "Input to 'structure' was `NULL` and 'data' was not identified as",
          "an `EGA` type object. Input 'data' or an `EGA` type object."
        ),
        call. = FALSE
      )

    }else if(get_object_type(structure) == "list"){

      # Determine if list (for hierarchical structures)

      # Check all levels
      level_lengths <- nvapply(structure, length)

      # Check if first level length matches number of variables
      length_error(structure[[1]], variables, "tefi")

      # Iterate through levels to check consistency
      for(level in seq_len(length(structure) - 1)){

        # Set higher index
        higher_index <- level + 1

        # Get length of higher meberships
        n_higher <- length(higher_memberships)

        # Collect memberships
        lower_memberships <- structure[[level]]
        higher_memberships <- structure[[higher_index]]

        # Get unique lower memberships
        lower_communities <- unique_length(lower_memberships)

        # Length should be equal to lower_communities or number of variables
        if(!(n_higher %in% c(lower_communities, variables))){
          stop(
            paste0(
              "Mismatch in hierarchical structure levels: level ", level,
              " has ", lower_communities, " communities but level ", higher_index,
              " does not have matching length (", n_higher, ")."
            ),
            call. = FALSE
          )
        }

        # If higher level only has communities (not expanded to all variables), expand
        if(n_higher == lower_communities){
          structure[[higher_index]] <- single_revalue_memberships(lower_memberships, higher_memberships)
        }

      }

    }else{ # Simple flat structure

      # Perform length check
      length_error(structure, variables, "tefi")

    }

  }else{

    # Get flag for hierarchical
    if(is(data, "hierEGA")){  # Use internal `hierEGA_structure` from `itemStability`
      structure <- hierEGA_structure(ega_object, structure)
    }else if(is.null(structure)){ # Use EGA memberships
      structure <- ega_object$wc
    }else{ # Ensure proper length
      length_error(structure, length(ega_object$wc), "tefi")
    }

  }

  # Convert if string
  if(is.list(structure)){

    # Check for characters
    structure <- lapply(
      structure, function(x){

        # If characters, then convert to numeric
        swiftelse(is.character(x), as.numeric(reindex_memberships(x)), x)

      }
    )

  }else if(is.character(structure)){
    structure <- as.numeric(reindex_memberships(structure))
  }

  # Return structure
  return(structure)

}

#' @noRd
# `tefi` standard function ----
# Updated 07.08.2023
tefi_standard <- function(correlation_matrix, structure, verbose)
{

  # Check structure
  if(anyNA(structure)){

    # Determine variables that are NA
    rm.vars <- is.na(structure)

    # Send warning message
    if(verbose){
      warning(
        paste(
          "Some variables did not belong to a dimension:",
          paste0(dimnames(correlation_matrix)[[2]][rm.vars], collapse = ", "),
          "\n\nUse caution: These variables have been removed from the TEFI calculation"
        ), call. = FALSE
      )
    }

    # Keep variables
    keep_vars <- !rm.vars

    # Keep available variables
    correlation_matrix <- correlation_matrix[keep_vars, keep_vars]

    # Remove NAs from structure
    structure <- structure[keep_vars]

  }

  # Obtain Von Neumann's entropy of density matrix
  H_vn <- matrix_entropy(correlation_matrix / dim(correlation_matrix)[2L])

  # Obtain communities
  communities <- unique_length(structure)

  # Get Von Neumman entropy by community
  H_vn_wc <- nvapply(seq_len(communities), function(community){

    # Get indices
    indices <- structure == community

    # Get community matrix
    community_matrix <- correlation_matrix[indices, indices]

    # Return Von Neumann entropy
    return(matrix_entropy(community_matrix / dim(community_matrix)[2L]))

  })

  # Pre-compute values
  ## Mean of community Von Neumann
  mean_H_vn_wc <- mean(H_vn_wc, na.rm = TRUE)
  ## Sum of community Von Neumann
  sum_H_vn_wc <- mean_H_vn_wc * communities
  ## Difference between total and total community
  # H_diff <- H_vn - sum_H_vn_wc
  ## Average entropy
  mean_H_vn <- mean_H_vn_wc - H_vn

  # Set up results
  return(
    fast.data.frame(
      data = c(
        mean_H_vn + ((H_vn - sum_H_vn_wc) * sqrt(communities)),
        sum_H_vn_wc - H_vn,
        mean_H_vn
      ), ncol = 3,
      colnames = c(
        "VN.Entropy.Fit", "Total.Correlation", "Average.Entropy"
      )
    )
  )

}

#' @noRd
# `tefi` generalized function ----
# Updated 08.05.2025
generalized_tefi <- function(correlation_matrix, structure, verbose = FALSE)
{

  # Get variables
  variables <- dim(correlation_matrix)[2L]

  # Determine NAs across all structure levels
  NA_memberships <- rowSums(lvapply(structure, is.na, LENGTH = variables))

  # Determine variables that are NA
  rm.vars <- NA_memberships != 0

  # Check structure
  if(any(rm.vars)){

    if(verbose){
      warning(
        paste(
          "Some variables did not belong to a dimension:",
          paste0(dimnames(correlation_matrix)[[2]][rm.vars], collapse = ", "),
          "\n\nUse caution: These variables have been removed from the TEFI calculation"
        ), call. = FALSE
      )
    }

    # Keep only available variables
    keep_vars <- !rm.vars
    correlation_matrix <- correlation_matrix[keep_vars, keep_vars]

    # Remove NAs from structure
    structure <- lapply(structure, function(x){x[keep_vars]})

  }

  # Von Neumann entropy of total structure
  H_vn <- matrix_entropy(correlation_matrix / dim(correlation_matrix)[2L])

  # Get number of levels
  n_levels <- length(structure)
  level_sequence <- seq_len(n_levels)

  # Prepare storage for each level's summed entropy
  summed_entropies <- num_communities <- numeric(n_levels)

  # Loop over each level of structure
  for(level in level_sequence){

    # Set index
    index <- structure[[level]]

    # Unique communities at this level
    communities <- unique_length(index)
    num_communities[[level]] <- communities

    # Von Neumann entropy for each community
    H_vn_communities <- nvapply(seq_len(communities), function(community){

      indices <- index == community
      community_matrix <- correlation_matrix[indices, indices]

      return(matrix_entropy(community_matrix / dim(community_matrix)[2L]))

    })

    # Compute sum of the entropies
    summed_entropies[[level]] <- sum(H_vn_communities, na.rm = TRUE)

  }

  # Compute generalized TEFI

  # Total sum of entropies across all levels
  total_entropy_sum <- sum(summed_entropies)

  # Scaling factor (sqrt of first-order number of communities)
  sqrt_first_order <- sqrt(num_communities[1L])

  # Pre-multiply levels by entropy
  total_entropy <- n_levels * H_vn

  # Return results
  return(
    fast.data.frame(
      c(
        # Generalized TEFI
        (total_entropy_sum / num_communities[1L]) - total_entropy +
        (total_entropy - total_entropy_sum) * sqrt_first_order,
        # Individual TEFIs
        (summed_entropies / num_communities[1L]) - H_vn +
        (H_vn - summed_entropies) * sqrt_first_order
      ),
      ncol = n_levels + 1,
      colnames = c("VN.Entropy.Fit", paste0("Level_", level_sequence, "_VN"))
    )
  )

}
