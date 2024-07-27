#' @title Jensen-Shannon Distance
#'
#' @description Computes the Jensen-Shannon Distance between two networks
#'
#' @param network1 Matrix or data frame.
#' Network to be compared
#'
#' @param network2 Matrix or data frame.
#' Second network to be compared
#'
#' @param method Character (length = 1).
#' Method to compute Jensen-Shannon Distance.
#' Defaults to \code{"spectral"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"kld"} --- Uses Kullback-Leibler Divergence
#'
#' \item \code{"spectral"} --- Uses eigenvalues of combinatorial Laplacian matrix to compute
#' Von Neumann entropy
#'
#' }
#'
#' @param signed Boolean. (length = 1).
#' Should networks be remain signed?
#' Defaults to \code{TRUE}
#'
#' @examples
#' # Obtain wmt2 data
#' wmt <- wmt2[,7:24]
#'
#' # Set seed (for reproducibility)
#' set.seed(1234)
#'
#' # Split data
#' split1 <- sample(
#'   1:nrow(wmt), floor(nrow(wmt) / 2)
#' )
#' split2 <- setdiff(1:nrow(wmt), split1)
#'
#' # Obtain split data
#' data1 <- wmt[split1,]
#' data2 <- wmt[split2,]
#'
#' # Perform EBICglasso
#' glas1 <- EBICglasso.qgraph(data1)
#' glas2 <- EBICglasso.qgraph(data2)
#'
#' # Spectral JSD
#' jsd(glas1, glas2)
#' # 0.1595893
#'
#' # Spectral JSS (similarity)
#' 1 - jsd(glas1, glas2)
#' # 0.8404107
#'
#' # Jensen-Shannon Divergence
#' jsd(glas1, glas2, method = "kld")
#' # 0.1393621
#'
#' @return Returns Jensen-Shannon Distance
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
#'
# Jensen-Shannon Distance
# Updated 27.07.2024
jsd <- function(
    network1, network2,
    method = c("kld", "spectral"),
    signed = TRUE
)
{

  # Check for missing arguments (argument, default, function)
  method <- set_default(method, "spectral", jsd)

  # Argument errors (send back networks in case of tibble)
  error_return <- jsd_errors(network1, network2, signed)

  # Get networks
  network1 <- error_return$network1; network2 <- error_return$network2

  # Check for method
  if(method == "spectral"){

    # Obtain rescaled Laplacian matrices
    laplacian1 <- rescaled_laplacian(network1)
    laplacian2 <- rescaled_laplacian(network2)

    # Obtain individual VN entropies
    lentropy1 <- entropy_laplacian(laplacian1)
    lentropy2 <- entropy_laplacian(laplacian2)

    # Obtain combined VN entropy
    lentropy_combined <- entropy_laplacian(0.50 * (laplacian1 + laplacian2))

    # Compute JSD
    JSD <- sqrt(abs(lentropy_combined - (0.50 * (lentropy1 + lentropy2))))

  }else if(method == "kld"){

    # Pre-compute inverse covariance matrix of combined networks
    inverse_combined <- pcor2inv(0.50 * (network1 + network2))

    # Compute KLDs
    kld1 <- kld(inverse_combined, pcor2inv(network1))
    kld2 <- kld(inverse_combined, pcor2inv(network2))

    # Compute JSD
    JSD <- 0.50 * kld1 + 0.50 * kld2

  }

  # Return (ensure real numbers)
  return(Re(JSD))

}

#' @noRd
# Argument errors ----
# Updated 27.07.2024
jsd_errors <- function(network1, network2, signed)
{

  # 'network1' errors
  object_error(network1, c("matrix", "data.frame", "tibble"), "jsd")

  # Check for tibble
  if(!is(network1, "matrix")){
    network1 <- as.matrix(network1)
  }

  # 'network2' errors
  object_error(network2, c("matrix", "data.frame"), "jsd")

  # Check for tibble
  if(!is(network2, "matrix")){
    network2 <- as.matrix(network2)
  }

  # 'Check for 'signed' errors
  length_error(signed, 1, "jsd")
  typeof_error(signed, "logical", "jsd")

  # Return networks
  return(
    list(
      network1 = swiftelse(signed, network1, abs(network1)),
      network2 = swiftelse(signed, network2, abs(network2))
    )
  )

}

#' @noRd
# Rescaled Laplacian matrix ----
# Updated 10.07.2023
rescaled_laplacian <- function(network)
{

  # Ensure diagonal is zero
  diag(network) <- 0

  # Get node strength
  node_strength <- colSums(network, na.rm = TRUE)

  # Return
  return((diag(node_strength) - network) / sum(node_strength, na.rm = TRUE))

}

#' @noRd
# Von Neumann Entropy ----
# Called "entropy_laplacian" to avoid conflict with `vn.entropy`
# Updated 11.10.2023
entropy_laplacian <- function(laplacian_matrix)
{

  # Set NAs to zero
  if(anyNA(laplacian_matrix)){
    laplacian_matrix[is.na(laplacian_matrix)] <- 0
  }

  # Get eigenvalues
  eigenvalues <- eigen(laplacian_matrix, symmetric = TRUE, only.values = TRUE)$values

  # Return entropy
  return(silent_call(-sum(eigenvalues * log2(eigenvalues), na.rm = TRUE)))

}

#' @noRd
# Kullback-Leibler Divergence ----
# Updated 03.07.2023
# Compute Kullback-Leibler Divergence
kld <- function(network1, network2)
{

  # network1 = P
  # network2 = Q
  # KLD(P || Q)

  # Pre-compute matrix multiplication
  combined_network <- crossprod(solve(network1), network2)

  # Return KLD
  return(
    trace(combined_network) -
    log2(det(combined_network)) -
    dim(network1)[2]
  )

}

#' @noRd
# Faster pairwise spectral JSD
# Updated 27.07.2024
pairwise_spectral_JSD <- function(network_list, ...)
{

  # Get ellipse
  ellipse <- list(...)

  # Get length of list
  network_length <- length(network_list)

  # Get ID names
  ID_names <- names(network_list)

  # Initialize matrix
  jsd_matrix <- matrix(
    nrow = network_length, ncol = network_length,
    dimnames = list(ID_names, ID_names)
  )

  # Check for signed
  if("signed" %in% names(ellipse) && !ellipse$signed){

    # Get absolute networks
    network_list <- lapply(network_list, abs)

  }

  # Pre-compute rescaled Laplacian matrices
  rescaled_L <- lapply(network_list, rescaled_laplacian)

  # Pre-compute entropies
  H <- nvapply(rescaled_L, entropy_laplacian)

  # Loop over networks
  for(i in seq_len(network_length)){
    for(j in i:network_length){

      # Obtain combined VN entropy
      lentropy_combined <- entropy_laplacian(
        0.50 * (rescaled_L[[i]] + rescaled_L[[j]])
      )

      # Compute JSD
      jsd_matrix[i,j] <- jsd_matrix[j,i] <- Re(
        sqrt(abs(lentropy_combined - (0.50 * (H[i] + H[j]))))
      )

    }
  }

  # Return matrix
  return(jsd_matrix)

}


#' @noRd
# Faster comparison spectral JSD
# Updated 27.07.2024
comparison_spectral_JSD <- function(base, network_list, ...)
{

  # Get ellipse
  ellipse <- list(...)

  # Check for signed
  if("signed" %in% names(ellipse) && !ellipse$signed){

    # Get absolute networks
    base <- abs(base)
    network_list <- lapply(network_list, abs)

  }

  # Pre-compute rescaled Laplacian and entropy for base
  rescaled_base <- rescaled_laplacian(base)
  H_base <- entropy_laplacian(rescaled_base)

  # Return JSDs with comparison network list
  return(
    nvapply(network_list, function(network){

      # Get rescaled network and entropy
      rescaled_network <- rescaled_laplacian(network)
      H_network <- entropy_laplacian(rescaled_network)

      # Obtain combined VN entropy
      lentropy_combined <- entropy_laplacian(
        0.50 * (rescaled_base + rescaled_network)
      )

      # Return JSD
      return(
        Re(sqrt(abs(lentropy_combined - (0.50 * (H_base + H_network)))))
      )

    })
  )

}
