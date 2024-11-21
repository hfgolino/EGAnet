#' @title Entropy Fit Index using Von Neumman's entropy (Quantum Information Theory) for correlation matrices
#'
#' @description Computes the fit of a dimensionality structure using Von Neumman's 
#' entropy when the input is a correlation matrix. Lower values suggest better 
#' fit of a structure to the data
#'
#' @param data Matrix or data frame.
#' Contains variables to be used in the analysis
#'
#' @param structure Numeric or character vector (length = \code{ncol(data)}).
#' A vector representing the structure (numbers or labels for each item).
#' Can be theoretical factors or the structure detected by \code{\link[EGAnet]{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{VN.Entropy.Fit}{The Entropy Fit Index using Von Neumman's entropy}
#'
#' \item{Total.Correlation}{The total correlation of the dataset}
#'
#' \item{Average.Entropy}{The average entropy of the dataset}
#'
#' @examples
#' # Get EGA result
#' ega.wmt <- EGA(
#'   data = wmt2[,7:24], model = "glasso",
#'   plot.EGA = FALSE # no plot for CRAN checks
#' )
#' 
#' # Compute Von Neumman entropy
#' vn.entropy(ega.wmt$correlation, ega.wmt$wc)
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
#' 
# VN Entropy Function ----
# Updated 07.08.2023
vn.entropy <- function(data, structure)
{
  
  # Argument errors (return data in case of tibble)
  data <- vn.entropy_errors(data, structure)
  
  # Generic function to get necessary inputs
  output <- obtain_sample_correlations(
    data = data, n = 1, # set to 1 to avoid error
    corr = "auto", na.data = "pairwise", 
    verbose = FALSE
  )
  
  # Get absolute correlations
  data <- abs(output$correlation_matrix)
  
  # Check structure
  if(anyNA(structure)){
    
    # Determine variables that are NA
    rm.vars <- is.na(structure)
    
    # Send warning message
    warning(
      paste(
        "Some variables did not belong to a dimension:",
        paste0(dimnames(data)[[2]][rm.vars], collapse = ", ")
      ), call. = FALSE
    )
    message("Use caution: These variables have been removed from the TEFI calculation")
    
    # Keep available variables
    data <- data[!rm.vars, !rm.vars]
    
  }
  
  # Obtain communities
  communities <- unique_length(structure)
  
  # Obtain eigenvalues of density matrix
  eigenvalues <- matrix_eigenvalues(data / dim(data)[2])
 
  # Obtain Von Neumann's entropy
  H_vn <- entropy(eigenvalues)
  
  # Obtain eigenvalues by community
  eigenvalues_wc <- lapply(seq_len(communities), function(community){
    
    # Get indices
    indices <- structure == community
    
    # Get community matrix
    community_matrix <- data[indices, indices]

    # Return eigenvalues
    return(matrix_eigenvalues(community_matrix / dim(community_matrix)[2]))
    
  })
  
  # Get Von Neumann's entropy by community
  H_vn_wc <- nvapply(eigenvalues_wc, entropy)
  
  # Get eigenvalue Kronecker product
  eigenvalue_kronecker <- Reduce("%x%", eigenvalues_wc)
  
  # Get Von Neumman's joint entropy
  H_vn_joint <- entropy(eigenvalue_kronecker)
  
  # Pre-compute values
  ## Get average entropy
  mean_H_vn <- mean(H_vn_wc, na.rm = TRUE)
  ## Compute average entropy
  H_average <- mean_H_vn - H_vn_joint
  
  # Set up results
  return(
    fast.data.frame(
      c(
        H_average - ((H_vn - mean_H_vn) * sqrt(communities)),
        mean_H_vn * communities - H_vn_joint,
        # uses `mean_H_vn * communities` which is faster than `sum(H_vn_wc)`
        H_average
      ), ncol = 3,
      colnames = c(
        "VN.Entropy.Fit", "Total.Correlation", "Average.Entropy"
      )
    )
  )
  
}

#' @noRd
# Argument errors ----
# Updated 13.08.2023
vn.entropy_errors <- function(data, structure)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "vn.entropy")
  
  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }
  
  # 'structure' errors
  object_error(structure, "vector", "vn.entropy")
  length_error(structure, dim(data)[2], "vn.entropy")
  
  # Return usable data in case of tibble
  return(usable_data(data, verbose = TRUE))
  
}

# Bug Checking ----
# ## Basic input
# data <- wmt2[,7:24]; ega.wmt <- EGA(data, plot.EGA = FALSE)
# data <- ega.wmt$correlation
# structure <- ega.wmt$wc

# # Old Code (possible errors)
# vn.entropy <- function(data, structure){
#   if(!ncol(data)==nrow(data)){
#     data <- qgraph::cor_auto(data)
#   }
#   uniq <- unique(structure)
#   num.comm <- structure
#   data <- abs(data)
#     cor1 <- data/ncol(data)
#     eigen1 <- eigen(cor1)$values
#     h.vn <- -sum(eigen1*log(eigen1))
# 
#     n <- max(structure)
#     cor.fact <- vector("list")
#     eigen.fact <- vector("list")
#     l.eigen.fact <- vector("list")
#     h.vn.fact <- vector("list")
#     for(i in 1:n){
#       cor.fact[[i]] <- data[which(structure==unique(structure)[i]),which(structure==unique(structure)[i])]
#       cor.fact[[i]] <- cor.fact[[i]]/ncol(cor.fact[[i]])
#       eigen.fact[[i]] <- eigen(cor.fact[[i]])$values
#       l.eigen.fact[[i]] <- eigen.fact[[i]]*log(eigen.fact[[i]])
#       h.vn.fact[[i]] <- -sum(l.eigen.fact[[i]])
#     }
#     # Joint entropy using the eigenvalues of a Kronecker product of a list of matrices
#     cor.joint <- vector("list", n)
#     for(i in 1:n){
#       cor.joint[[i]] <- cor1[which(num.comm==uniq[i]),which(num.comm==uniq[i])]/table(num.comm)[[i]]
#       # ERROR HERE? POTENTIALLY TWO?
#       #
#       # Error 1: `cor1` is already a density matrix; shouldn't it be `data` instead?
#       # `data` would be the correlation matrix that is then divided by the number of
#       # variables to get the density matrix (which what I believe we need)
#       #
#       # Error 2: `table(num.comm)[[i]]` does not necessarily always correspond with `uniq[i]`
#       # and instead should be `table(num.comm)[[uniq[i]]]` to ensure correspondence
#     }
# 
#     eigen.val <- lapply(cor.joint, function(x) eigen(x)$values)
#     eigen.kronecker <- eigen.val[[1]]
#     for (i in 2:n){
#       eigen.kronecker <- eigen.kronecker %x% eigen.val[[i]]
#     }
# 
#     h.vn.joint <- -sum(eigen.kronecker*log(eigen.kronecker))
# 
#   h.vn.fact2 <- unlist(h.vn.fact)
# 
#   # Difference between Max the sum of the factor entropies:
#   Hdiff <- h.vn-mean(h.vn.fact2)
#   results <- data.frame(matrix(NA, nrow = 1, ncol = 3))
#   colnames(results) <- c("VN.Entropy.Fit", "Total.Correlation","Average.Entropy")
#   results$VN.Entropy.Fit <- ((mean(h.vn.fact2)-h.vn.joint))-((Hdiff-(mean(h.vn.fact2))*sqrt(n)))
#
#   ERROR HERE?
#   Based on equation from paper: 
#   (mean(h.vn.fact2) - h.vn.joint) - ((h.vn - mean(h.vn.fact2)) * sqrt(n))
#
#   Expansion on what's given in the original code
#   (mean(h.vn.fact2)-h.vn.joint) - ((h.vn - mean(h.vn.fact2)) - mean(h.vn.fact2) * sqrt(n))
#
#   results$Total.Correlation <- sum(h.vn.fact2)-h.vn.joint
#   results$Average.Entropy <- mean(h.vn.fact2)-h.vn.joint
#   return(results)
# }
