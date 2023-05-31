#%%%%%%%%%%%%%%%%%%%%%%%%
# FUNCTION FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Standard EGA arguments
# Updated 02.02.2023
ega_arguments <- function(arguments)
{
  
  # Set default arguments
  default_arguments <- list(
    data = NULL, n = NULL,
    corr = "cor_auto", uni.method = "louvain",
    model = "glasso", model.args = list(),
    algorithm = "walktrap", algorithm.args = list(),
    consensus.method = "most_common_tefi", consensus.iter = 100,
    plot.EGA = TRUE, plot.args = list()
  )
  
  # Match arguments with input
  matched_arguments <- match(names(arguments), names(default_arguments))
  
  # Check for whether any are not missing
  if(any(!is.na(matched_arguments))){
    
    # Remove missing
    replace_arguments <- na.omit(matched_arguments)
    
    # Replace arguments
    default_arguments[matched_arguments] <- arguments
    
  }
  
  # Return arguments
  return(default_arguments)
  
}

#' @noRd
# Standard {lavaan} arguments
# Updated 02.02.2023
lavaan_arguments <- function(arguments)
{
  
  # Set default arguments
  default_arguments <- list(
    model = NULL, data = NULL, ordered = NULL,
    sampling.weights = NULL, sample.cov = NULL,
    sample.mean = NULL, sample.th = NULL,
    sample.nobs = NULL, group = NULL, cluster = NULL,
    constraints = "", WLS.V = NULL, NACOV = NULL,
    std.lv = TRUE
  )
  
  # Match arguments with input
  matched_arguments <- match(names(arguments), names(default_arguments))
  
  # Check for whether any are not missing
  if(any(!is.na(matched_arguments))){
    
    # Remove missing
    replace_arguments <- na.omit(matched_arguments)
    
    # Replace arguments
    default_arguments[matched_arguments] <- arguments
    
  }
  
  # Return arguments
  return(default_arguments)
  
}

#' @noRd
# Make unidimensional CFA model
# Updated 02.02.2023
make_unidimensional_cfa <- function(variable_names)
{
 
  # Check for number of variables
  if(length(variable_names) == 2){
    
    # Make latent factor model
    # Fix all loadings to 1
    model <- paste0(
      "LF =~ ",
      paste0(
        "1*", variable_names, collapse = " + "
      )
    )
    
    # Fix residual variances to same value
    model <- paste(
      model,
      paste0(
        variable_names,
        "~~theta*",
        variable_names,
        collapse = " \n "
      ),
      sep = " \n "
    )
    
  }else{
    
    # Make latent factor model
    model <- paste0(
      "LF =~ ",
      paste(
        variable_names, collapse = " + "
      )
    )
    
  }
  
  # Return model
  return(model)
  
}

#' @noRd
# Determine estimator arguments
# Updated 03.23.2023
estimator_arguments <- function(lavaan_ARGS)
{
  
  # Obtain categories
  categories <- data_categories(
    data = lavaan_ARGS$data
  )
  
  # Check for categories
  if(any(categories <= 7)){
    
    # Set arguments
    lavaan_ARGS$estimator <- "WLSMV"
    lavaan_ARGS$missing <- "pairwise"
    lavaan_ARGS$ordered <- names(categories)[
      categories <= 7
    ]
    
  }else{
    
    # Set arguments
    lavaan_ARGS$estimator <- "MLR"
    lavaan_ARGS$missing <- "fiml"
    lavaan_ARGS$ordered <- FALSE
    
  }
  
  # Return arguments
  return(lavaan_ARGS)
  
  
}

#' @noRd
# Removes MASS package dependency (from version 7.3.54)
# Updated 23.12.2021
MASS_mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE, EISPACK = FALSE)
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  if (EISPACK)
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    drop(X)
  else t(X)
}


