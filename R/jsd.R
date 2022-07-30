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
#' # Compute Jensen-Shannon Distance
#' jsd(glas1, glas2) # 0.1832602
#' 
#' # Compute Jensen-Shannon Similarity
#' 1 - jsd(glas1, glas2) # 0.8167398
#'
#' @return Returns Jensen-Shannon Distance
#'
#' @author Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
#'
#' @export
#' 
# Jensen-Shannon Distance
# Updated 30.07.2022
jsd <- function(
    network1, network2,
    method = c("kld", "spectral")
)
{
  
  # Missing method
  if(missing(method)){
    method <- "spectral"
  }else{method <- tolower(match.arg(method))}
  
  # Check for method
  if(method == "spectral"){
    
    # Obtain rescaled Laplacian matrices
    rL1 <- rescaled_laplacian(network1)
    rL2 <- rescaled_laplacian(network2)
    
    # Obtain individual VN entropies
    vn1 <- vn_entropy(rL1)
    vn2 <- vn_entropy(rL2)
    
    # Obtain combined VN entropy
    rL_comb <- 0.5 * (rL1 + rL2)
    vn_comb <- vn_entropy(rL_comb)
    
    # Compute JSD
    JSD <- sqrt(
      vn_comb - (0.5 * (vn1 + vn2))
    )
    
  }else if(method == "kld"){
    
    # Compute Kullback-Leibler Divergence
    kld <- function(network1, network2){
      sum(diag(network1 %*% solve(network2))) -
        log2(det(network1 %*% solve(network2))) -
        ncol(network1)
    }
    
    # Combine networks
    network_comb <- 0.5 * (network1 + network2)
    
    # Compute KLDs
    kld1 <- kld(network1, network_comb)
    kld2 <- kld(network2, network_comb)
    
    # Compute JSD
    JSD <- sqrt(
      0.5 * kld1 + 0.5 * kld2
    )
    
  }
  
  # Return
  return(Re(JSD))
  
}
