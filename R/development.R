#' @noRd
# Polytomous IRT parameters ----
# Updated 26.07.2023
poly.irt <- function(loadings, data)
{
  
  # Check for {EGAnet} network loadings
  if(is(loadings, "net.loads")){
    
    # Check for rotated loadings
    if(is.null(loadings$rotated)){
      loadings <- loadings$std
    }else{
      loadings <- loadings$rotated
    }
    
  }
  
  # Ensure matrix (and in same order as data)
  loadings <- as.matrix(loadings)[dimnames(data)[[2]],]
  
  # Unique variance
  s <- sqrt(1 - rowSums(loadings^2))
  
  # Estimate discrimination parameters
  est_a <- sweep(loadings, 1, s, "/")
  
  # Estimate threshold parameters
  thresholds <- lapply(as.data.frame(data), obtain_thresholds)
  
  # Estimate location parameters
  est_d <- lapply(seq_along(thresholds), function(variable){
    -thresholds[[variable]] / s[variable]
  })
  names(est_d) <- names(thresholds)
  
  # Return results
  return(
    list(
      discrimination = est_a * 1.702,
      location = est_d
    )
  )
  
}
