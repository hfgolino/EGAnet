#%%%%%%%%%%%%%%%%%%%%
# MATH FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Trace of matrix
# Updated 27.05.2023
trace <- function(object)
{
  
  # Ensure matrix
  if(!is(object, "matrix") | !is(object, "table")){
    object <- as.matrix(object)
  }

  # Return trace
  return(sum(diag(object)))
  
}
