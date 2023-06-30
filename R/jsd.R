#' Jensen-Shannon Distance
#' 
#' @description Computes the Jensen-Shannon Distance between two networks
#'
#' @param network1 Matrix or data frame.
#' Network to be compared
#' 
#' @param network2 Matrix or data frame.
#' Second network to be compared
#'
#' @param method Character.
#' Method to compute Jensen-Shannon Distance.
#' Defaults to \code{"spectral"}.
#' Options:
#' 
#' \itemize{
#' 
#' \item{\code{"kld"}}
#' {Uses Kullback-Leibler Divergence}
#' 
#' \item{\code{"spectral"}}
#' {Uses eigenvalues of combinatiorial Laplacian matrix to compute
#' Von Neumann entropy}
#' 
#' }
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
#' jsd(glas1, glas2) # 0.1618195
#' 
#' # Spectral JSS (similarity)
#' 1 - jsd(glas1, glas2) # 0.8381805
#' 
#' # Jensen-Shannon Divergence
#' jsd(glas1, glas2, method = "kld") # 0.1923636
#'
#' @return Returns Jensen-Shannon Distance
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
#' 
# Jensen-Shannon Distance
# Updated 29.06.2023
jsd <- function(
    network1, network2,
    method = c("kld", "spectral")
)
{
  
  # Check for missing arguments (argument, default, function)
  method <- set_default(method, "spectral", jsd)
  
  # Check for method
  if(method == "spectral"){
    
    # Obtain rescaled Laplacian matrices
    laplacian1 <- rescaled_laplacian(network1)
    laplacian2 <- rescaled_laplacian(network2)
    
    # Obtain individual VN entropies
    lentropy1 <- silent_call(entropy_laplacian(laplacian1))
    lentropy2 <- silent_call(entropy_laplacian(laplacian2))
    
    # Obtain combined VN entropy
    laplacian_combined <- 0.5 * (laplacian1 + laplacian2)
    lentropy_combined <- silent_call(entropy_laplacian(laplacian_combined))
    
    # Compute JSD
    JSD <- sqrt(abs(lentropy_combined - (0.5 * (lentropy1 + lentropy2))))
    
  }else if(method == "kld"){
    
    # Combine networks
    network_combined <- 0.5 * (network1 + network2)
    
    # Pre-compute inverse covariance matrix of combined networks
    inverse_combined <- pcor2inv(network_combined)
    
    # Compute KLDs
    kld1 <- kld(inverse_combined, pcor2inv(network1))
    kld2 <- kld(inverse_combined, pcor2inv(network2))
    
    # Compute JSD
    JSD <- 0.5 * kld1 + 0.5 * kld2
    
  }
  
  # Return (ensure real numbers)
  return(Re(JSD))
  
}

#' @noRd
# Rescaled Laplacian matrix ----
# Updated 29.06.2023
rescaled_laplacian <- function(network)
{
  # Ensure diagonal is zero
  diag(network) <- 0
  
  # Return
  return((diag(colSums(network)) - network) / sum(network))
  
}

#' @noRd
# Von Neumann Entropy ----
# Called "entropy_laplacian" to avoid conflict with `vn.entropy`
# Updated 29.06.2023
entropy_laplacian <- function(laplacian_matrix)
{
  return(entropy(matrix_eigenvalues(laplacian_matrix), base = 2))
}

#' @noRd
# Kullback-Leibler Divergence ----
# Updated 29.06.2023
# Compute Kullback-Leibler Divergence
kld <- function(network1, network2)
{
  
  # network1 = P
  # network2 = Q
  # KLD(P || Q)
  
  # Compute inverse of first network
  inverse_network1 <- solve(network1)
  
  # Pre-compute matrix multiplication
  combined_network <- inverse_network1 %*% network2
  
  # Return KLD
  return(
    trace(combined_network) -
    log2(det(combined_network)) -
    dim(network1)[2]
  )
  
}