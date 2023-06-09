#%%%%%%%%%%%%%%%%%%%%%%%
# GENERAL FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Determine number of digits in a number ----
# Updated 27.05.2023
digits <- function(number)
{
  # Obtain the lowest value of log base 10 and add 1
  return(floor(log10(number)) + 1)
  
}

#' @noRd
# Force vector ----
# (usually for data frame rows/columns)
# Updated 27.05.2023
force_vector <- function(desired_vector)
{
  # Convert to matrix then vector
  new_vector <- as.vector(as.matrix(desired_vector))
  
  # Keep names (if any)
  if(!is.null(colnames(desired_vector))){
    names(new_vector) <- colnames(desired_vector)
  }
  
  return(new_vector)
  
}

#' @noRd
# Force matrix ----
# (usually for vectors)
# Updated 09.06.2023
force_matrix <- function(desired_matrix)
{
  
  # Check for matrix form already
  if(is(desired_matrix, "matrix")){
    return(desired_matrix)
  }else if(is(desired_matrix, "data.frame")){
    return(as.matrix(desired_matrix))
  }else{# Convert vector to matrix
    return(matrix(desired_matrix, ncol = 1))
  }
  
}

#' @noRd
# Determine number of categories in data ----
# Updated 02.02.2023
data_categories <- function(data)
{
  
  # Ensure data is matrix
  data <- as.matrix(data)
  
  # Loop over columns
  categories <- apply(
    data, 2, function(x){
      length(na.omit(unique(x)))
    }
  )
  
  # Return categories
  return(categories)
  
}

#' @noRd
# Convert version to number ----
# Updated 02.02.2023
version_conversion <- function(version)
{
  
  # Convert to character
  version <- as.character(version)
  
  # Remove periods
  version <- gsub("\\.", "", version)
  
  # Convert to numeric
  version <- as.numeric(version)
  
  # Return version
  return(version)
  
}

#' @noRd
# All-purpose symmetric checker ----
# Updated 03.02.2023
is_symmetric <- function(data){
  
  # Check for whether rows equal columns
  if(nrow(data) == ncol(data)){
    
    # Convert to matrix
    data_matrix <- as.matrix(data)
    
    # Remove names
    data_matrix <- unname(data_matrix)
    
    # Check that lower triangle equal upper triangle
    lower_triangle <- data_matrix[lower.tri(data_matrix)]
    transpose_matrix <- t(data_matrix) # ensures similar orientation
    upper_triangle <- transpose_matrix[lower.tri(transpose_matrix)]
    
    # Check that all are equal
    all_equal <- all(lower_triangle == upper_triangle, na.rm = TRUE)
    
  }else{
    
    # Not a matrix
    return(FALSE)
    
  }
  
  # Return whether all are equal
  return(all_equal)
  
}

#' @noRd
# Ensure data has dimension names ----
# Updated 03.02.2023
ensure_dimension_names <- function(data)
{
  
  # Check for column names
  if(is.null(colnames(data))){
    
    # Standardize names
    colnames(data) <- paste0(
      "V", formatC(
        x = 1:ncol(data),
        digits = (digits(ncol(data)) - 1),
        format = "d", flag = "0"
      )
    )
    
  }
  
  # Check for matrix
  if(nrow(data) == ncol(data)){
    
    # Check for row names
    if(is.null(data) | any(row.names(data) != colnames(data))){
      
      # Assign column names to row names
      row.names(data) <- colnames(data)
      
    }
    
  }
  
  # Return named data
  return(data)
  
}

#' @noRd
# No names print ----
# Updated 03.02.2023
no_name_print <- function(object){
  
  # Convert object to data frame
  df <- as.data.frame(object)
  
  # Remove column names
  colnames(df) <- NULL
  
  # Print with no quotes or row names
  print(df, quote = FALSE, row.names = FALSE)

}

#' @noRd
#'
# General function to check for packages
# Updated 12.04.2023
check_package <- function(packages)
{
  
  # Check for packages
  installed <- packages %in% row.names(installed.packages())
  
  # Determine which packages are not installed
  not_installed <- packages[!installed]
  
  # Print error with missing packages
  if(length(not_installed) != 0){
    
    # Organize packages error output
    if(length(not_installed) > 1){
      
      # Get missing packages
      missing_packages <- paste0("{", packages , "}", collapse = ", ")
      packages <- paste0("\"", packages, "\"", collapse = ", ")
    
      # Stop and tell user to install package
      stop(
        paste0(
          missing_packages, 
          " are not installed but are required for this function. ",
          "Please run \n\n",
          "install.packages(c(", packages, "))",
          "\n\nOnce installed, re-run this function (you may need to restart R/RStudio)."
        )
      )
      
    }else{
      
      # Get missing packages
      missing_packages <- paste0("{", packages, "}")
      packages <- paste0("\"", packages, "\"")
      
      # Stop and tell user to install package
      stop(
        paste0(
          missing_packages, 
          " is not installed but is required for this function. ",
          "Please run \n\n",
          "install.packages(c(", packages, "))",
          "\n\nOnce installed, re-run this function (you may need to restart R/RStudio)."
        )
      )
      
    }
    
  }
  
}

#' @noRd
#'
# General function to silently obtain output
# Updated 11.05.2023
silent_call <- function(...){
  
  # Make call
  sink <- capture.output(
    result <- suppressWarnings(
      suppressMessages(
        ...
      )
    )
  )
  
  # Return result
  return(result)
  
}

















