#' Exploratory Graph Model
#'
#' @description Function to fit the Exploratory Graph Model
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' Can be raw data or a correlation matrix
#'
#' @param EGM.model Character vector (length = 1).
#' Sets the procedure to conduct \code{EGM}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"explore"} --- A two-step procedure where the first step estimates
#' network loadings for Walktrap community structures up to \code{communities} using the
#' standard beta-min condition to threshold the network. The second step computes networks
#' using the community-aware beta-min criterion based on the
#' loading-implied partial correlation matrices from the previous step. The
#' network (corresponding to communities and a set of loadings) with the lowest
#' BIC is selected as the model. An optional (but recommended and default) procedure is to
#' optimize the network toward the model-implied zero-order correlation matrix
#' (`network.optimize = TRUE`); otherwise, the community-aware beta-min criterion is used
#' on the model-implied partial correlation matrix.
#'
#' \item \code{"EGA"} (default) --- Applies \code{\link[EGAnet]{EGA}} to obtain the
#' (sparse) regularized network structure, communities, and memberships
#'
#' \item \code{"probability"} --- Estimates communities based on the non-regularized
#' empirical partial correlation matrix and sparsity is set using probabilities for the
#' within- (\code{p.in}) and between-community (\code{p.out}) connectivity
#'
#' }
#'
#' @param communities Numeric vector (length = 1).
#' Number of communities to use for the \code{"standard"} type of EGM.
#' Defaults to \code{NULL}.
#' Providing no input will use the communities and memberships output
#' from the Walktrap algorithm (\code{\link[igraph]{cluster_walktrap}}) based
#' on the empirical non-regularized partial correlation matrix
#'
#' @param structure Numeric or character vector (length = \code{ncol(data)}).
#' Can be theoretical factors or the structure detected by \code{\link[EGAnet]{EGA}}.
#' Defaults to \code{NULL}
#'
#' @param search Boolean (length = 1).
#' Whether a search over parameters should be conducted.
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} to select a model over a variety of parameters that
#' optimizes the \code{opt} objective
#'
#' @param p.in Numeric vector (length = 1).
#' Probability that a node is randomly linked to other nodes in the same community.
#' Within community edges are set to zero based on \code{quantile(x, prob = 1 - p.in)}
#' ensuring the lowest edge values are set to zero (i.e., most probable to \emph{not}
#' be randomly connected).
#' Only used for \code{EGM.type = "standard"}.
#' Defaults to \code{NULL} but must be set
#'
#' @param p.out Numeric vector (length = 1).
#' Probability that a node is randomly linked to other nodes \emph{not} in the same community.
#' Between community edges are set to zero based on \code{quantile(x, prob = 1 - p.out)}
#' ensuring the lowest edge values are set to zero (i.e., most probable to \emph{not}
#' be randomly connected).
#' Only used for \code{EGM.type = "standard"} and \code{search = FALSE}.
#' Defaults to \code{NULL} but must be set
#'
#' @param opt Character vector (length = 1).
#' Fit index used to select from when searching over models
#' (only applies to \code{search = TRUE}).
#' Available options include maximum likelihood (default; \code{"logLik"}) and
#' standardized root mean residual (\code{"SRMR"})
#'
#' @param model.select Character vector (length = 1).
#' Criterion to select the best fitting model.
#' Defaults to \code{"AICq"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"logLik"} --- log-likelihood; not recommended as more communities
#' will \emph{almost} always have better likelihood; \eqn{\log{(\hat{L})}}
#'
#' \item \code{"AIC"} --- Akaike's information criterion;
#' \eqn{-2\log{(\hat{L})} \cdot 2k} where
#' \emph{k} is the number of parameters (non-zero edges plus non-loadings)
#'
#' \item \code{"AICc"} --- Corrected AIC (adjusts for sample size);
#' \eqn{AIC + \frac{2k^2 + 2k}{n - k - 1}} where \emph{n} is sample size
#'
#' \item \code{"AICq"} (default)--- Modularity adjusted AIC;
#' \eqn{AIC + 6Q\log{n}} where
#' \emph{Q} is modularity (see below). The constant of 6 adjust modularity
#' such that it contributes roughly 4-6 AIC for each increase in 0.10 modularity
#' across most usual sample sizes (between 500-10000)
#'
#' \item \code{"BIC"} --- Bayesian information criterion;
#' \eqn{k\log{n} - 2\log{(\hat{L})}}
#'
#' \item \code{"Q"} --- Modularity;
#' \eqn{\frac{1}{2m} \sum_{ij} (A_{ij} - \frac{k_i k_j}{2m}) \delta(c_i, c_j)} where
#' \emph{m} is the sum total edges of the network, \emph{A}, \emph{k} is the
#' node strength, and \emph{c} is the membership of nodes \emph{i} and \emph{j}
#'
#' }
#'
#' @param constrain.structure Boolean (length = 1).
#' Whether memberships of the communities should
#' be added as a constraint when optimizing the network loadings.
#' Defaults to \code{TRUE} which ensures assigned loadings are
#' guaranteed to never be smaller than any cross-loadings.
#' Set to \code{FALSE} to freely estimate each loading similar to exploratory factor analysis
#'
#' @param constrain.zeros Boolean (length = 1).
#' Whether zeros in the estimated network loading matrix should
#' be retained when optimizing the network loadings.
#' Defaults to \code{TRUE} which ensures that zero networks loadings are retained.
#' Set to \code{FALSE} to freely estimate each loading similar to exploratory factor analysis
#'
#' @param verbose Boolean (length = 1).
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not display progress
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}},
#' \code{\link[EGAnet]{community.unidimensional}},
#' \code{\link[EGAnet]{EGA}}, and
#' \code{\link[EGAnet]{net.loads}}
#'
#' @examples
#' # Get depression data
#' data <- depression[,24:44]
#'
#' \dontrun{
#' # Estimate EGM (using "explore")
#' egm <- EGM(data)
#'
#' # Constrain EGM to only estimate 3 communities
#' egm_communities <- EGM(data, communities = 3)
#'
#' # Estimate EGM (using EGA) specifying structure
#' egm_ega_structure <- EGM(
#'   data, structure = c(
#'     1, 1, 1, 2, 1, 1, 1,
#'     1, 1, 1, 3, 2, 2, 2,
#'     2, 3, 3, 3, 3, 3, 2
#'   ), EGM.model = "EGA"
#' )
#'
#'
#' # Estimate EGM (using EGA search)
#' egm_ega_search <- EGM(data, EGM.model = "EGA", search = TRUE)
#'
#' # Estimate EGM (using search)
#' egm_prob_search <- EGM(
#'   data, EGM.model = "probability", search = TRUE,
#'   communities = 3, # need communities or structure
#'   p.in = 0.95 # only need 'p.in'
#' )}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Estimate EGM ----
# Updated 21.06.2025
EGM <- function(
    data, EGM.model = c("explore", "EGA", "probability"),
    communities = NULL, structure = NULL, search = FALSE,
    p.in = NULL, p.out = NULL, opt = c("logLik", "SRMR"),
    model.select = c("logLik", "AIC", "AICc", "AICq", "BIC", "Q"),
    constrain.structure = TRUE, constrain.zeros = TRUE,
    verbose = TRUE, ...
)
{

  # Set default
  EGM.model <- set_default(EGM.model, "explore", EGM)
  opt <- set_default(opt, "loglik", EGM)
  model.select <- set_default(model.select, "aicq", EGM)

  # Set up EGM type internally
  EGM.type <- switch(
    EGM.model,
    "explore" = "explore",
    "ega" = swiftelse(search, "ega.search", "ega"),
    "probability" = swiftelse(search, "search", "probability")
  )

  # Check for communities if model is 'explore'
  if(EGM.type == "explore" & is.null(communities)){

    # Ensure 'search' and 'communities' are set
    search <- TRUE; communities <- 8

  }

  # Check data and structure
  data <- EGM_errors(
    data, EGM.type, communities, search, p.in, p.out,
    constrain.structure, constrain.zeros,
    verbose, ...
  )

  # Switch and return results based on type
  return(
    switch(
      EGM.type,
      "explore" = EGM.explore(data, communities, search, opt, model.select, ...),
      "ega" = EGM.EGA(data, structure, opt, constrain.structure, constrain.zeros, ...),
      "ega.search" = EGM.EGA.search(data, communities, structure, opt, constrain.structure, constrain.zeros, verbose, ...),
      "probability" = EGM.probability(data, communities, structure, p.in, p.out, opt, constrain.structure, constrain.zeros, ...),
      "search" = EGM.search(data, communities, structure, p.in, opt, constrain.structure, constrain.zeros, verbose, ...)
    )
  )

}

#' @noRd
# EGM Errors ----
# Updated 20.06.2025
EGM_errors <- function(
    data, EGM.type, communities, search, p.in, p.out,
    constrain.structure, constrain.zeros,
    verbose, ...
)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "EGM")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # Check for NULL communities
  if(!is.null(communities) & EGM.type != "explore"){

    # If not NULL, check 'communities' errors
    typeof_error(communities, "numeric", "EGM")
    object_error(communities, "vector", "EGM")
    length_error(communities, 1, "EGM")

  }

  # 'search' errors
  length_error(search, 1, "EGM")
  typeof_error(search, "logical", "EGM")

  # Check first for parameters involved in both search and standard
  if(grepl("probability", EGM.type)){

    # Check for NULL in 'p.in'
    if(is.null(p.in)){
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "Input for 'p.in' is `NULL`. A value between 0 and 1 for 'p.in' must be provided."
        ),
        call = "EGM"
      )
    }

    # Check 'p.in' errors
    typeof_error(p.in, "numeric", "EGM")
    range_error(p.in, c(0, 1), "EGM")
    length_error(p.in, c(1, communities), "EGM")

  }

  # 'constrain.structure' errors
  length_error(constrain.structure, 1, "EGM")
  typeof_error(constrain.structure, "logical", "EGM")

  # 'constrain.zeros' errors
  length_error(constrain.zeros, 1, "EGM")
  typeof_error(constrain.zeros, "logical", "EGM")

  # Check for EGM type
  if(EGM.type == "search"){

    # 'verbose' errors
    length_error(verbose, 1, "EGM")
    typeof_error(verbose, "logical", "EGM")

  }else if(EGM.type == "standard"){

    # Check for NULL in 'p.out'
    if(is.null(p.out)){
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "Input for 'p.out' is `NULL`. A value between 0 and 1 for 'p.out' must be provided."
        ),
        call = "EGM"
      )
    }

    # Check 'p.out' errors
    typeof_error(p.out, "numeric", "EGM")
    range_error(p.out, c(0, 1), "EGM")
    length_error(p.out, c(1, communities), "EGM")

  }

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, FALSE)
  }

  # Return data in case of tibble
  return(data)

}

#' @exportS3Method
# S3 Print Method ----
# Updated 09.11.2024
print.EGM <- function(x, digits = 3, ...)
{

  # Return EGA results
  print(x$EGA)

  # Add break
  cat("\n\n----\n\n")

  # Send title
  cat("Fit based on Optimized Loadings\n\n")

  # Add fit metrics
  print(format_decimal(x$model$optimized$fit, digits), quote = FALSE)

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 09.11.2024
summary.EGM <- function(object, digits = 3, ...)
{
  print(object) # same as print
}

#' @exportS3Method
# S3 Plot Method ----
# Updated 12.11.2024
plot.EGM <- function(x, ...)
{

  # Return EGA plot
  plot(x$EGA, ...)

}

# LOADINGS AND CORRELATIONS ----

#' @noRd
# Community correlations from loadings
# Updated 21.07.2025
community_correlations <- function(simple_structure, loading_structure)
{
  return(cov2cor(crossprod(simple_structure, nload2cor(loading_structure)) %*% simple_structure))
}

#' @noRd
# Network loadings to partial correlations ----
# Updated 20.03.2025
nload2pcor <- function(loadings)
{

  # Compute covariance matrix
  P <- tcrossprod(loadings)

  # Compute interdependence
  diag(P) <- interdependence <- sqrt(diag(P))

  # Compute matrix I
  I <- diag(sqrt(1 / interdependence))

  # Compute zero-order correlations
  R <- I %*% P %*% I

  # Obtain inverse covariance matrix
  INV <- solve(R)

  # Compute matrix D
  D <- diag(sqrt(1 / diag(INV)))

  # Compute partial correlations
  P <- -D %*% INV %*% D

  # Set diagonal to zero
  diag(P) <- 0

  # Return partial correlation
  return(P)

}

#' @noRd
# Network loadings to correlations ----
# Updated 20.03.2025
nload2cor <- function(loadings)
{

  # Compute partial correlation
  P <- tcrossprod(loadings)

  # Compute interdependence
  diag(P) <- interdependence <- sqrt(diag(P))

  # Compute matrix I
  I <- diag(sqrt(1 / interdependence))

  # Return correlation
  return(I %*% P %*% I)

}

#' @noRd
# Compute model parameters for TEFI adjustment ----
# Updated 22.10.2024
compute_tefi_adjustment <- function(loadings, correlations)
{

  # Obtain dimensions
  dimensions <- dim(loadings)

  # Total number of parameters
  parameters <- (dimensions[1] * dimensions[2]) + dimensions[1] +
    ((dimensions[2] * (dimensions[2] - 1)) / 2)

  # Obtain model parameters
  model_parameters <- parameters - sum(loadings == 0)

  # Return adjustment
  return(-2 * log(model_parameters) + mean(abs(correlations)))

}

# EGM TYPES ----

#' @noRd
# EGM | Exploratory ----
# Updated 20.06.2025
EGM.explore <- function(data, communities, search, opt, model.select, ...)
{

  # Obtain data dimensions
  data_dimensions <- dim(data)

  # Estimate correlations
  empirical_R <- auto.correlate(data, ...)
  empirical_P <- cor2pcor(empirical_R)
  diag(empirical_P) <- 0

  # Obtain variable names from the correlations
  variable_names <- dimnames(empirical_R)[[2]]

  # Need at least two variables per community
  at_least_two <- round(data_dimensions[2] / 2) - 1
  max_communities <- max(communities)

  # Determine communities
  community_sequence <- swiftelse(
    search, seq_len(
      swiftelse(
        max_communities > at_least_two,
        at_least_two, max_communities
      )
    ), max_communities
  )

  # Obtain expected edges
  EE <- expected_edges(empirical_P, data_dimensions)

  # Obtain partial correlations that are greater than chance
  null_P <- empirical_P * (abs(empirical_P) > (EE - 1.96 * attr(EE, "SE")))

  # Collect results
  results <- lapply(
    community_sequence, EGM.explore.core, null_P = null_P,
    cluster = hclust(d = as.dist(1 - sqrt(abs(null_P))), method = "average"),
    variable_names = variable_names, data_dimensions = data_dimensions,
    empirical_R = empirical_R, opt = opt
  )

  # Obtain fits
  fits <- do.call(rbind.data.frame, lapply(results, function(x){x$fit}))

  # Get converged fits
  converged <- grepl("0", fits$converged)
  good_fits <- fits[converged,, drop = FALSE] # ensure matrix
  results <- results[converged]

  # Extract optimal results
  optimal_results <- results[[which.min(good_fits[[model.select]])]]

  # Obtain structure
  structure <- structure(optimal_results$structure, names = variable_names)

  # Set communities
  communities <- unique_length(structure)
  community_sequence <- seq_len(communities)

  # Set simple loading structure
  simple_loadings <- optimal_results$loadings
  for(i in community_sequence){
    simple_loadings[structure != i, i] <- 0
  }

  # Compute network scores
  optimized_scores <- scale(data) %*% simple_loadings

  # Compute community correlations
  optimized_correlations <- cor(optimized_scores, use = "pairwise")

  # Obtain model-implied correlations
  optimized_R <- pcor2cor(optimal_results$network)
  optimized_P <- cor2pcor(optimized_R)

  # Set up EGA
  ega_list <- list(
    dim.variables = data.frame(
      items = variable_names,
      dimension = structure
    ),
    network = optimal_results$network, wc = structure,
    n.dim = communities, correlation = empirical_R,
    n = data_dimensions[1], TEFI = tefi(empirical_R, structure)$VN.Entropy.Fit
  ); class(ega_list) <- "EGA"

  # Attach methods to network
  attr(ega_list$network, which = "methods") <- list(
    model = "egm", communities = dim(optimal_results$loadings)[2],
    lambda = optimal_results$lambda
  )

  # Attach class to memberships
  class(ega_list$wc) <- "EGA.community"

  # Attach methods to memberships
  attr(ega_list$wc, which = "methods") <- list(
    algorithm = swiftelse(is.null(structure), NULL, "EGM Explore"),
    objective_function = NULL
  )

  # Set up results
  results <- list(
    EGA = ega_list, structure = structure,
    model = list( # general list
      standard = list( # using standard parameters
        loadings = NULL,
        scores = NULL,
        correlations = NULL,
        fit = NULL,
        implied = list(R = NULL, P = NULL)
      ),
      optimized = list( # using optimized parameters
        loadings = optimal_results$loadings,
        scores = optimized_scores,
        correlations = optimized_correlations,
        fit = fit(
          n = data_dimensions[1], p = data_dimensions[2],
          R = optimized_R, S = empirical_R,
          loadings = optimal_results$loadings,
          correlations = optimized_correlations,
          structure = structure,
          ci = 0.95, remove_correlations = FALSE
        ),
        implied = list(R = optimized_R, P = optimized_P),
        fits = fits
      )
    )
  )

  # Set class
  class(results) <- "EGM"

  # Return results
  return(results)

}

#' @noRd
# EGM | EGA ----
# Updated 02.04.2024
EGM.EGA <- function(data, structure, opt, constrain.structure, constrain.zeros, ...)
{

  # Obtain data dimensions
  data_dimensions <- dim(data)

  # Estimate EGA
  ega <- EGA(data, plot.EGA = FALSE, ...)
  empirical_P <- cor2pcor(ega$correlation)

  # Update number of rows based on EGA
  data_dimensions[1] <- ega$n

  # Obtain variable names from the network
  variable_names <- dimnames(ega$network)[[2]]

  # Set memberships based on structure
  if(!is.null(structure)){

    # Update EGA features
    ega$wc[] <- structure
    ega$dim.variables$dimension <- structure

  }else{
    structure <- ega$wc
  }

  # Set communities
  ega$n.dim <- communities <- unique_length(structure)

  # Obtain standard network loadings
  output <- silent_call(net.loads(A = ega$network, wc = ega$wc, ...))

  # Obtain standard loadings
  standard_loadings <- output$std[variable_names,, drop = FALSE]

  # Compute network scores
  standard_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  standard_correlations <- cor(standard_scores, use = "pairwise")

  # Obtain model-implied correlations
  standard_R <- nload2cor(standard_loadings)
  standard_P <- cor2pcor(standard_R)

  # Get loading dimensions
  dimensions <- dim(standard_loadings)
  dimension_names <- dimnames(standard_loadings)

  # Obtain loadings vector and get bounds
  loadings_vector <- as.vector(standard_loadings)
  loadings_length <- length(loadings_vector)
  zeros <- swiftelse(constrain.zeros, loadings_vector != 0, rep(1, loadings_length))

  # Set up loading structure
  loading_structure <- matrix(
    FALSE, nrow = dimensions[1],
    ncol = dimensions[2],
    dimnames = list(dimension_names[[1]], dimension_names[[2]])
  )

  # Fill structure
  for(i in seq_along(structure)){
    loading_structure[i, structure[i]] <- TRUE
  }

  # Optimize over loadings
  result <- egm_optimize(
    loadings_vector = loadings_vector,
    loadings_length = loadings_length,
    zeros = zeros, R = ega$correlation,
    loading_structure = loading_structure,
    rows = communities, n = data_dimensions[1],
    v = data_dimensions[2],
    constrained = constrain.structure,
    lower_triangle = lower.tri(ega$correlation),
    opt = opt
  )

  # Extract optimized loadings
  optimized_loadings <- matrix(
    result$par, nrow = dimensions[1],
    dimnames = dimension_names
  )

  # Replace loadings output with optimized loadings
  output$loadings$std <- optimized_loadings

  # Compute network scores
  optimized_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  optimized_correlations <- cor(optimized_scores, use = "pairwise")

  # Obtain model-implied correlations
  optimized_R <- nload2cor(optimized_loadings)
  optimized_P <- cor2pcor(optimized_R)

  # Set up results
  results <- list(
    EGA = ega, structure = structure,
    model = list( # general list
      standard = list( # using standard parameters
        loadings = standard_loadings,
        scores = standard_scores,
        correlations = standard_correlations,
        fit = fit(
          n = data_dimensions[1], p = data_dimensions[2],
          R = standard_R, S = ega$correlation,
          loadings = output$std, correlations = standard_correlations,
          structure = ega$wc, ci = 0.95, remove_correlations = FALSE
        ),
        implied = list(R = standard_R, P = standard_P)
      ),
      optimized = list( # using optimized parameters
        loadings = optimized_loadings,
        scores = optimized_scores,
        correlations = optimized_correlations,
        fit = fit(
          n = data_dimensions[1], p = data_dimensions[2],
          R = optimized_R, S = ega$correlation,
          loadings = optimized_loadings,
          correlations = optimized_correlations,
          structure = ega$wc, ci = 0.95, remove_correlations = FALSE
        ),
        implied = list(R = optimized_R, P = optimized_P)
      )
    )
  )

  # Set class
  class(results) <- "EGM"

  # Return results
  return(results)

}

# EGM | EGA search ----
# Updated 20.03.2025
EGM.EGA.search <- function(data, communities, structure, opt, constrain.structure, constrain.zeros, verbose, ...)
{

  # Get number of dimensions
  data_dimensions <- dim(data)

  # Put data through `EBICglasso.qgraph`
  glasso_output <- network.estimation(data, network.only = FALSE, ...)

  # Obtain attributes
  glasso_attr <- attributes(glasso_output$estimated_network)

  # Obtain fits
  fits <- lapply(
    seq_len(glasso_attr$methods$nlambda), function(i){

      # Convert to partial correlations
      P <- wi2net(glasso_output$output$results$wi[,,i])
      dimnames(P) <- glasso_attr$dimnames

      # Obtain organized output
      output <- try(
        silent_call(
          glasso_fit(
            data = data, S = glasso_output$output$S,
            glasso_attr = glasso_attr, P = P,
            data_dimensions = data_dimensions,
            constrain.structure = constrain.structure,
            constrain.zeros = constrain.zeros, opt = opt,
            ...
          )
        ), silent = TRUE
      )

      # Return output
      return(swiftelse(is(output, "try-error"), NULL, output))

    }
  )

  # Add names to fits
  names(fits) <- glasso_output$output$lambda

  # Fit name based on fit
  fit_name <- switch(
    opt,
    "loglik" = "logLik",
    "srmr" = "SRMR"
  )

  # Obtain fit index
  fit_index <- nvapply(fits[!lvapply(fits, is.null)], function(x){
    return(x$model$optimized$fit[[fit_name]])
  })

  # Index function based on fit
  index_FUN <- switch(
    opt,
    "loglik" = max,
    "srmr" = min
  )

  # Obtain optimum value
  opt_value <- index_FUN(fit_index)

  # Obtain index
  opt_index <- max(which(fit_index == opt_value))

  # Obtain lambda (as character)
  lambda <- names(fit_index)[[opt_index]]

  # Obtain results
  results <- fits[[lambda]]

  # Set up 'dim.variables' for EGA
  dim.variables <- fast.data.frame(
    c(dimnames(results$EGA$network)[[2]], structure),
    nrow = data_dimensions[2], ncol = 2,
    colnames = c("items", "dimension")
  )

  # Add to EGA results
  results$EGA$dim.variables <- dim.variables[order(dim.variables$dimension),]

  # Compute TEFI for EGA
  results$EGA$TEFI <- tefi(glasso_output$output$S, results$EGA$wc)$VN.Entropy.Fit

  # Set up network methods
  attributes(results$EGA$network) <- glasso_attr
  attributes(results$EGA$network)$methods[
    c("model.selection", "lambda", "criterion")
  ] <- list(
    model.selection = fit_name, lambda = as.numeric(lambda), criterion = opt_value
  )

  # Overwrite class
  class(results) <- "EGM"

  # Return results
  return(results)

}

#' @noRd
# EGM | Probability ----
# Updated 04.04.2025
EGM.probability <- function(data, communities, structure, p.in, p.out, opt, constrain.structure, constrain.zeros, ...)
{

  # Get ellipse
  ellipse <- list(...)

  # Get dimensions
  data_dimensions <- dim(data)
  dimension_names <- dimnames(data)

  # Determine if sample size was provided
  if(data_dimensions[1] == data_dimensions[2]){

    # Check if sample size was provided
    if("n" %in% names(ellipse)){

      # Set sample size
      data_dimensions[1] <- ellipse$n

      # Set empirical zero-order correlations
      empirical_R <- data

    }else{

      # Stop and send error
      .handleSimpleError(
        h = stop,
        msg = paste0(
          "A symmetric (", data_dimensions[1], " x ",
          data_dimensions[2], ") matrix was input but sample size was not provided.\n\n",
          "If you'd like to use a correlation matrix, please set the argument 'n' to your sample size"
        ),
        call = "EGM"
      )

    }

  }else{

    # Estimate zero-order correlations
    empirical_R <- auto.correlate(data, ...)

  }

  # Obtain partial correlations
  empirical_P <- cor2pcor(empirical_R)

  # Check for whether structure is provided
  if(is.null(structure)){

    # Obtain structure if not provided
    walktrap <- igraph::cluster_walktrap(convert2igraph(abs(empirical_P)))

    # Obtain memberships based on number of communities
    structure <- swiftelse(
      is.null(communities) || unique_length(walktrap$membership) == communities,
      walktrap$membership, cutree(as.hclust(walktrap), k = communities)
    )

  }

  # Set number of communities
  communities <- unique_length(structure)

  # Set up community variables
  community_variables <- lapply(
    seq_len(communities), function(community){structure == community}
  )

  # Currently, sets up for `p.in` and `p.out` to correspond to each
  # individual community
  #
  # This might make sense for particularly simple structures but other
  # options might exist
  #
  # For example:
  #
  # # Set up community variables
  # community_variables <- do.call(
  #   cbind, lapply(
  #     seq_len(communities), function(community){structure == community}
  #   )
  # )
  #
  # # Create "on" and "off" block structure
  # on_block <- tcrossprod(community_variables)
  # off_block <- !on_block
  #
  # This approach allows the block structure to use `p.in` as a general
  # thresholding scheme across communities rather than for *each* community
  #
  # Similarly, the off-block sets `p.out` to threshold across all non-community
  # edges rather than on *each* non-community block

  # Obtain community structure
  community_P <- create_community_structure(
    P = empirical_P, total_variables = data_dimensions[2],
    communities = communities,
    community_variables = community_variables,
    p.in = p.in, p.out = p.out
  )

  # Update loadings
  output <- silent_call(net.loads(community_P, structure, ...))
  output$std <- output$std[dimension_names[[2]],, drop = FALSE]

  # Obtain dimension names
  loading_names <- dimnames(output$std)

  # Compute network scores
  standard_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  standard_correlations <- cor(standard_scores, use = "pairwise")

  # Compute model-implied correlations
  standard_R <- nload2cor(output$std)
  standard_P <- cor2pcor(standard_R)

  # Obtain loadings vector and get bounds
  loadings_vector <- as.vector(output$std)
  loadings_length <- length(loadings_vector)
  zeros <- swiftelse(constrain.zeros, loadings_vector != 0, rep(1, loadings_length))

  # Set up loading structure
  # Uses transpose for 2x speed up in optimization
  loading_structure <- matrix(FALSE, nrow = communities, ncol = data_dimensions[2])

  # Fill structure
  for(i in seq_along(structure)){
    loading_structure[structure[i], i] <- TRUE
  }

  # Optimize over loadings
  result <- egm_optimize(
    loadings_vector = loadings_vector,
    loadings_length = loadings_length,
    zeros = zeros, R = empirical_R,
    loading_structure = loading_structure,
    rows = communities, n = data_dimensions[1],
    v = data_dimensions[2],
    constrained = constrain.structure,
    lower_triangle = lower.tri(empirical_R),
    opt = opt
  )

  # Extract optimized loadings
  optimized_loadings <- matrix(
    result$par, nrow = data_dimensions[2],
    dimnames = loading_names
  )

  # Replace loadings output with optimized loadings
  output$std <- optimized_loadings

  # Compute network scores
  optimized_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Compute community correlations
  optimized_correlations <- cor(optimized_scores, use = "pairwise")

  # Compute model-implied correlations
  optimized_R <- nload2cor(optimized_loadings)
  optimized_P <- cor2pcor(optimized_R)

  # Set up EGA
  ega_list <- list(
    dim.variables = data.frame(
      items = dimension_names[[2]],
      dimension = structure
    ),
    network = community_P, wc = structure,
    n.dim = unique_length(structure), correlation = empirical_R,
    n = data_dimensions[1], TEFI = tefi(empirical_R, structure)$VN.Entropy.Fit
  ); class(ega_list) <- "EGA"

  # Attach methods to network
  attr(ega_list$network, which = "methods") <- list(
    model = "egm", communities = communities, p.in = p.in, p.out = p.out
  )

  # Attach class to memberships
  class(ega_list$wc) <- "EGA.community"

  # Attach methods to memberships
  names(ega_list$wc) <- dimension_names[[2]]
  attr(ega_list$wc, which = "methods") <- list(
    algorithm = swiftelse(is.null(structure), NULL, "Walktrap"),
    objective_function = NULL
  )

  # Set up results
  results <- list(
    EGA = ega_list, structure = structure,
    model = list( # general list
      standard = list( # using standard parameters
        loadings = output$std,
        scores = standard_scores,
        correlations = standard_correlations,
        fit = fit(
          n = data_dimensions[1], p = data_dimensions[2],
          R = standard_R, S = empirical_R,
          loadings = output$std, correlations = standard_correlations,
          structure = structure, ci = 0.95, remove_correlations = FALSE
        ),
        implied = list(R = standard_R, P = standard_P)
      ),
      optimized = list( # using optimized parameters
        loadings = optimized_loadings,
        scores = optimized_scores,
        correlations = optimized_correlations,
        fit = fit(
          n = data_dimensions[1], p = data_dimensions[2],
          R = optimized_R, S = empirical_R,
          loadings = optimized_loadings, correlations = optimized_correlations,
          structure = structure, ci = 0.95, remove_correlations = FALSE
        ),
        implied = list(R = optimized_R, P = optimized_P)
      )
    )
  )

  # Set class
  class(results) <- "EGM"

  # Return results
  return(results)

}

#' @noRd
# EGM | Search ----
# Updated 20.03.2025
EGM.search <- function(data, communities, structure, p.in, opt, constrain.structure, constrain.zeros, verbose, ...)
{

  # Perform search based on 'p.in'
  p_grid <- expand.grid(
    p_in = seq(p.in, 1, 0.05), p_out = seq(0.00, p.in, 0.05)
  )

  # Get number of search
  p_number <- dim(p_grid)[1]

  # Loop over grid search
  grid_search <- parallel_process(
    iterations = p_number, datalist = seq_len(p_number),
    FUN = function(i, data, p_grid, communities, structure, opt, constrain.structure, constrain.zeros, ...){

      # First, try
      output <- silent_call(
        try(
          EGM.probability(
            data = data, communities = communities,
            structure = structure, p.in = p_grid$p_in[i],
            p.out = p_grid$p_out[i], opt = opt,
            constrain.structure = constrain.structure,
            constrain.zeros = constrain.zeros, ...
          ), silent = TRUE
        )
      )

      # Return result
      return(swiftelse(is(output, "try-error"), NULL, output))

    },
    # `EGM.probability` arguments
    data = data, p_grid = p_grid,
    communities = communities, structure = structure,
    opt = opt, constrain.structure = constrain.structure,
    constrain.zeros = constrain.zeros, ...,
    # `parallel_process` arguments
    ncores = 1, progress = verbose
  )

  # Fit name based on fit
  fit_name <- switch(
    opt,
    "aic" = "AIC", "bic" = "BIC",
    "loglik" = "logLik", "srmr" = "SRMR"
  )

  # Obtain fits
  optimized_fits <- nvapply(
    grid_search, function(x){
      swiftelse(is.null(x), NA, x$model$optimized$fit[[fit_name]])
    }
  )

  # Check for all NA
  if(all(is.na(optimized_fits))){

    # Set error if no solution is found
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Search could not converge to any solutions. \n\n",
        "Try:\n",
        "  +  increasing 'p.in'\n",
        "  +  changing 'communities'\n",
        "  +  changing 'structure'"
      ),
      call = "EGM"
    )

  }

  # Index function based on fit
  index_FUN <- switch(
    opt,
    "loglik" = which.max,
    "srmr" = which.min
  )

  # Obtain index
  opt_index <- index_FUN(optimized_fits)

  # Set up final model
  results <- grid_search[[opt_index]]

  # Add 'p.in' and 'p.out' parameters
  results$search <- c(
    p.in = p_grid[opt_index, "p_in"],
    p.out = p_grid[opt_index, "p_out"]
  )

  # Overwrite class
  class(results) <- "EGM"

  # Return results
  return(results)

}

# # EGM | Confirmatory ----
# # Updated 11.04.2025
# EGM.confirm <- function(data, structure, opt, ...)
# {
#
#   # Obtain data dimensions
#   data_dimensions <- dim(data)
#
#   # Set data names
#   colnames(data) <- paste0("V", seq_len(ncol(data)))
#
#   # Compute correlations
#   empirical_R <- auto.correlate(data, ...)
#   empirical_P <- cor2pcor(empirical_R)
#
#   # Obtain variable names from the network
#   variable_names <- dimnames(data)[[2]]
#
#   # Set up loading structure
#   loading_structure <- matrix(
#     FALSE, nrow = data_dimensions[2], ncol = 3,
#     dimnames = list(variable_names, format_integer(1:3, 0))
#   )
#
#   # Get dimension names
#   dimension_names <- dimnames(loading_structure)
#
#   # Fill structure
#   for(i in seq_along(structure)){
#     loading_structure[i, structure[i]] <- TRUE
#   }
#
#   structure <- loading_structure
#
#   # Assume structure is a matrix
#   communities <- dim(structure)[2]
#
#   # Set up structure as vector
#   structure_vector <- as.vector(structure)
#
#   # Obtain get bounds
#   structure_length <- length(structure_vector)
#
#   # Optimize over loadings
#   result <- egm_optimize(
#     loadings_vector = runif_xoshiro(structure_length, min = 0.00, max = 0.01),
#     loadings_length = structure_length,
#     zeros = structure_vector, R = empirical_R,
#     loading_structure = structure,
#     rows = communities, n = data_dimensions[1],
#     v = data_dimensions[2], constrained = TRUE,
#     lower_triangle = lower.tri(empirical_R),
#     opt = opt
#   )
#
#   # Extract optimized loadings
#   optimized_loadings <- matrix(
#     result$par, nrow = data_dimensions[2],
#     dimnames = dimension_names
#   )
#
#   # Replace loadings output with optimized loadings
#   ega <- EGA(data, plot.EGA = FALSE)
#   ega$wc <- max.col(structure)
#   output <- silent_call(net.loads(A = ega$network, wc = ega$wc))
#   output$loadings$std <- optimized_loadings
#
#   # Compute network scores
#   optimized_scores <- compute_scores(output, data, "network", "simple")$std.scores
#
#   # Compute community correlations
#   optimized_correlations <- cor(optimized_scores, use = "pairwise")
#
#   struct <- max.col(structure)
#
#   # Add cross-loadings based on correlations
#   for(i in 1:3)
#     for(j in 1:3){
#       if(i != j){
#
#         optimized_loadings[struct == i, j] <- 0.25 * optimized_correlations[i,j] /
#           ((1 - optimized_loadings[struct == i, i]) * sqrt(log(data_dimensions[2])))
#
#
#       }
#     }
#
#   # Obtain model-implied correlations
#   optimized_R <- nload2cor(optimized_loadings)
#   optimized_P <- cor2pcor(optimized_R)
#
#   # Set up results
#   results <- list(
#     EGA = ega, structure = structure,
#     model = list( # general list
#       standard = list( # using standard parameters
#         loadings = standard_loadings,
#         scores = standard_scores,
#         correlations = standard_correlations,
#         fit = fit(
#           n = data_dimensions[1], p = data_dimensions[2],
#           R = standard_R, S = ega$correlation,
#           loadings = output$std, correlations = standard_correlations,
#           structure = ega$wc, ci = 0.95, remove_correlations = FALSE
#         ),
#         implied = list(R = standard_R, P = standard_P)
#       ),
#       optimized = list( # using optimized parameters
#         loadings = optimized_loadings,
#         scores = optimized_scores,
#         correlations = optimized_correlations,
#         fit = fit(
#           n = data_dimensions[1], p = data_dimensions[2],
#           R = optimized_R, S = empirical_R,
#           loadings = structure,
#           correlations = optimized_correlations,
#           structure = max.col(structure), ci = 0.95,
#           remove_correlations = FALSE
#         ),
#         implied = list(R = optimized_R, P = optimized_P)
#       )
#     )
#   )
#
#   # Set class
#   class(results) <- "EGM"
#
#   # Return results
#   return(results)
#
# }

# EGM UTILITIES ----

#' @noRd
# EGM | Core Exploration ----
# Updated 13.07.2025
EGM.explore.core <- function(
    communities, null_P, cluster, variable_names,
    data_dimensions, empirical_R, opt, ...
)
{

  # Set bad fit from the git
  bad_fit <- data.frame(
    parameters = NA, loglik = NA, aic = NA, aicc = NA,
    bic = NA, q = NA, aicq = NA, min_eigenvalue = -Inf,
    converged = "1: error in loading estimation"
  )

  # Set memberships
  membership <- cutree(cluster, communities)

  # Compute loadings
  loadings <- silent_call(
    net.loads(A = null_P, wc = membership)$std[variable_names,, drop = FALSE]
  )

  # Check for single memberships
  singletons <- fast_table(membership) == 1

  # If singletons, replace loadings
  if(any(singletons)){

    # Loop over
    for(i in which(singletons)){

      # Set index
      index <- membership == i

      # Set singleton to max of overall connections (allow drop to vector)
      for(i in which(singletons)){
        loadings[,i] <- 0
        loadings[index, i] <- max(null_P[index,])
      }

    }

  }

  # Set up loadings vector
  loadings_vector <- as.vector(loadings)

  # Get loading dimensions
  dimensions <- dim(loadings)
  dimension_names <- dimnames(loadings)

  # Set up loading structure
  loading_structure <- matrix(
    TRUE, nrow = dimensions[1], ncol = dimensions[2],
    dimnames = dimension_names
  )

  # Allow zeros to be estimated
  zeros <- rep(1, length(loadings_vector))

  # Set lower triangle
  lower_triangle <- lower.tri(empirical_R)

  # Optimize for initial loadings
  initial_loadings <- loadings_optimization(
    iter = 10, loadings_vector = loadings_vector, zeros = zeros,
    R = empirical_R, loading_structure = loading_structure,
    communities = communities, data_dimensions = data_dimensions,
    lower_triangle = lower_triangle, opt = opt
  )

  # Check for bad result
  if(is.null(initial_loadings)){
    return(list(loadings = loadings, fit = bad_fit))
  }

  # Compute minimum eigenvalue of hessian
  min_eigenvalue <- min(matrix_eigenvalues(initial_loadings$hessian))

  # Format loadings
  loadings[] <- initial_loadings$par

  # Set lambda sequence
  # (reverse order ensures higher lambda is used when there is a tie)
  lambdas <- seq(0.20, 0.01, -0.01)

  # Optimize SCAD soft threshold
  soft_threshold <- nvapply(
    lambdas, scad_threshold_cost,
    loadings = loadings, empirical_R = empirical_R,
    lower_triangle = lower_triangle,
    data_dimensions = data_dimensions
  )

  # Get updated loadings
  loadings <- scad_threshold(loadings, lambdas[[which.min(soft_threshold)]])

  # Obtain solution
  membership <- max.col(abs(loadings))

  # Threshold network
  P <- set_network(loadings, membership, data_dimensions)

  # Get quality flags
  negative_flag <- min_eigenvalue < 0 & initial_loadings$convergence == 1
  meaningful_flag <- unique_length(membership) != dimensions[2]
  singleton_flag <- any(fast_table(membership) < 2)
  PD_flag <- anyNA(P) || !is_positive_definite(pcor2cor(P))

  # Check overall quality
  if(negative_flag || meaningful_flag || singleton_flag || PD_flag){

    # Check down reasons
    if(negative_flag){
      converged <- "1: negative hessian eigenvalue"
    }else if(meaningful_flag){
      converged <- "1: unique membership < communities"
    }else if(singleton_flag){
      converged <- "1: singleton communities"
    }else if(PD_flag){
      converged <- "1: not positive definite"
    }

    # Updated bad fit
    bad_fit[,c("min_eigenvalue", "converged")] <- c(min_eigenvalue, converged)

    # Return bad result
    return(list(loadings = loadings, fit = bad_fit))

  }

  # Add dimension names
  dimnames(P) <- dimnames(empirical_R)

  # Get implied correlations
  implied_R <- pcor2cor(P)

  # Compute log-likelihood
  logLik <- silent_call(
    -log_likelihood(
      n = data_dimensions[1], p = data_dimensions[2],
      R = implied_R, S = empirical_R, type = "zero"
    )
  )

  # Obtain loading parameters
  parameters <- sum(loadings != 0)

  # Compute negative 2 times log-likelihood
  logLik2 <- 2 * logLik

  # Compute 2 times parameters
  parameters2 <- 2 * parameters

  # Compute AIC
  aic <- logLik2 + parameters2

  # Compute modularity
  q <- obtain_modularity(P, membership)

  # Collect fit indices
  fit <- data.frame(
    parameters = parameters,
    loglik = logLik,
    aic = aic,
    aicc = aic + (parameters2 * (parameters + 1)) /
      (data_dimensions[2] - parameters2 - 1),
    bic = logLik2 + parameters * log(data_dimensions[1]),
    q = -q,
    aicq = aic - q * 40,
    min_eigenvalue = min_eigenvalue,
    converged = "0: successful"
  )

  # Return result
  return(list(network = P, structure = membership, fit = fit, loadings = loadings))

}

#' @noRd
# Expected edge values ----
# Updated 26.06.2025
expected_edges <- function(network, data_dimensions = NULL)
{

  # Compute node strength
  strength <- colSums(network)

  # Obtain the normalized cross-product
  EE <- tcrossprod(strength) / sum(strength)

  # Obtain differences as SE
  if(!is.null(data_dimensions)){

    # Delta method for the standard error

    # Ensure absolute
    network <- abs(network)

    # Convert network to Fisher's z
    z_network <- r2z(network)

    # Compute node strength
    strength <- colSums(z_network)
    total_strength <- sum(strength)
    standard_strength <- strength / total_strength
    total_squared <- total_strength^2
    cross_strength <- tcrossprod(strength)
    dEE <- cross_strength / total_squared

    # Get variables minus one
    p_minus_one <- data_dimensions[2] - 1

    # Because of Chung-Lu configuration model, all nodes are
    # assumed to have similar variability
    variance <- p_minus_one * (1 / (data_dimensions[1] - data_dimensions[2] - 1))

    # Compute gradient
    gradient_j <- (matrix(standard_strength, nrow = data_dimensions[2], ncol = data_dimensions[2]) - dEE)^2

    # Get SE (still in Fisher's z)
    # t(gradient_j) = shorthand for gradient i
    # dEE^2 = shorthand for gradient k
    SE <- sqrt(variance * (t(gradient_j) + gradient_j + (p_minus_one - 1) * dEE^2))

    # Absolute partial correlations (assumes similar transformation deviations)
    lower_triangle <- lower.tri(EE)
    average <- sum(EE[lower_triangle] * network[lower_triangle]) / sum(EE[lower_triangle])

    # Calculate the Jacobian
    attr(EE, "SE") <- SE * sqrt((1 - average^2)^2)

  }

  # Return result
  return(EE)

}

#' @noRd
# Obtain modularity edge values ----
# Updated 13.07.2025
obtain_modularity <- function(network, membership = NULL)
{

  # Ensure absolute
  network <- abs(network)

  # Set modularity value
  return(
    swiftelse(
      is.null(membership), 1,
      swiftelse(
        unique_length(membership) == 1,
        mean(network[lower.tri(network)]),
        modularity(network, membership)
      )
    )
  )

}

#' @noRd
# l1 soft threshold ----
# Updated 11.07.2025
l1_threshold <- function(loadings, lambda)
{
  return(sign(loadings) * pmax(abs(loadings) - lambda, 0))
}

#' @noRd
# SCAD soft threshold ----
# Updated 11.07.2025
scad_threshold <- function(loadings, lambda)
{

  # gamma = 3.7
  # tau = 1

  # Set absolute
  L <- abs(loadings)

  # Pre-compute components
  gamma_lambda <- 3.7 * lambda

  # Return values
  return(
    ifelse(
      L <= 2 * lambda,
      l1_threshold(loadings, lambda),
      ifelse(
        L <= gamma_lambda,
        1.588235 * l1_threshold(loadings, gamma_lambda / 2.7),
        loadings
      )
    )
  )

}

#' @noRd
# Threshold for network ----
# Updated 13.07.2025
set_network <- function(loadings, membership, data_dimensions)
{

  # Obtain partial correlations
  P <- nload2pcor(loadings)

  # Set edges to zero for zero loadings
  for(i in seq_len(data_dimensions[2])){

    # Check for zero loadings
    zero_index <- which(loadings[i,] == 0)

    # Identify whether zero loadings exist
    if(length(zero_index) != 0){

      # Loop over zero loading memberships
      for(community in zero_index){

        # Get index
        index <- membership == community

        # Set values to zero
        P[index, i] <- P[i, index] <- 0

      }

    }

  }

  # Return network
  return(P)

}

#' @noRd
# SCAD soft threshold ----
# Updated 13.07.2025
scad_threshold_cost <- function(lambda, loadings, empirical_R, lower_triangle, data_dimensions)
{

  # Apply SCAD threshold
  loadings <- scad_threshold(loadings, lambda)

  # Set membership
  membership <- max.col(abs(loadings))

  # Obtain model-implied partial correlations with SCAD
  P <- try(
    set_network(
      loadings = loadings, membership = membership,
      data_dimensions = data_dimensions
    ), silent = TRUE
  )

  # Check first for missing issues (all zero loadings)
  if(is(P, "try-error") || anyNA(P)){
    return(Inf)
  }

  # Convert partial correlations to zero-order correlations
  R <- try(pcor2cor(P), silent = TRUE)

  # Check for positive definite
  if(is(R, "try-error") || anyNA(R) || !is_positive_definite(R)){
    return(Inf)
  }else{

    # Set parameters
    parameters <- sum(loadings != 0)

    # Compute log-likelihood
    loglik <- log_likelihood(
      n = data_dimensions[1], p = data_dimensions[2],
      R = pcor2cor(P), S = empirical_R, type = "zero"
    )

    # Send result
    return(-2 * loglik + 2 * parameters - obtain_modularity(P, membership) * 40)

  }

}

#' @noRd
# Creates community structure ----
# Updated 08.10.2024
create_community_structure <- function(
    P, total_variables, communities, community_variables, p.in, p.out
)
{

  # Set vectors of 'p.in' and 'p.out'
  p.in <- swiftelse(length(p.in) == 1, rep(p.in, communities), p.in)
  p.out <- swiftelse(length(p.out) == 1, rep(p.out, communities), p.out)

  # Set diagonal to missing
  diag(P) <- NA

  # Set community blocks
  for(i in seq_len(communities)){

    # Randomly set zero in block
    indices <- P[community_variables[[i]], community_variables[[i]]]

    # Check for current sparsity
    if(compute_density(indices) > p.in[i]){

      # Get lower triangle
      lower_triangle <- lower.tri(indices)

      # Sample to set to zero
      indices[
        abs(indices) < quantile(
          abs(indices[lower_triangle]), probs = 1 - p.in[i], na.rm = TRUE
        )
      ] <- 0

      # Set back into block
      P[community_variables[[i]], community_variables[[i]]] <- indices

    }

    # Set between-community indices
    indices <- P[community_variables[[i]], -unlist(community_variables[-i])]

    # Check for current sparsity
    if(compute_density(indices) > p.out[i]){

      # Get threshold value
      threshold_value <- quantile(abs(indices), probs = 1 - p.out[i], na.rm = TRUE)

      # Set below to zero
      indices[abs(indices) < threshold_value] <- 0

      # Set back into block
      P[community_variables[[i]], -unlist(community_variables[-i])] <- indices

      # Do other side
      indices <- P[-unlist(community_variables[-i]), community_variables[[i]]]

      # Set below to zero
      indices[abs(indices) < threshold_value] <- 0

      # Set back into block
      P[-unlist(community_variables[-i]), community_variables[[i]]] <- indices

    }

  }

  # Set diagonal to zero
  diag(P) <- 0

  # Return community network
  return(P)

}

#' @noRd
# Known Graph ----
# Follows pages 631--634:
# Hastie, T., Tibshirani, R., & Friedman, J. (2008).
# The elements of statistical learning: Data mining, inference, and prediction (2nd ed.).
# New York, NY: Springer.
# Updated 24.12.2024
known_graph <- function(S, A, tol = 1e-06, max_iter = 100)
{

  # Obtain number of nodes
  nodes <- dim(S)[2]

  # Get node sequence
  node_sequence <- seq_len(nodes)

  # Initialize matrix
  W <- S

  # Initialize zero betas
  zero_beta <- matrix(0, nrow = nodes - 1, ncol = 1)

  # Obtain non-zero edge list
  non_zero_list <- lapply(node_sequence, function(i){
    A[-i, i] != 0
  })

  # Initialize convergence criteria
  iteration <- 0; difference <- Inf

  # Loop until criterion is met or max iterations is reached
  while(difference > tol && iteration < max_iter){

    # Increase iterations
    iteration <- iteration + 1

    # Set old W
    W_old <- W

    # Iterate over each node
    for(j in node_sequence){

      # Identify edges in the adjacency matrix for node j
      edges <- non_zero_list[[j]]

      # Initialize beta
      beta <- zero_beta

      # Partition for W11
      W11 <- W[-j, -j, drop = FALSE]

      # Extract the corresponding rows/columns in W11 and elements of S12
      if(any(edges)){

        # Compute S12
        S12 <- S[-j, j, drop = FALSE]

        # Solve for beta* = (W11*)^{-1} * S12*
        beta[edges,] <- solve(
          W11[edges, edges, drop = FALSE], S12[edges, , drop = FALSE]
        )

        # Update W12 = W11 * beta
        W12_new <- W11 %*% beta

        # Update the corresponding parts of W
        W[j, -j] <- W[-j, j] <- W12_new
        W[j, j] <- S[j, j] - sum(beta * W12_new)

      }

    }

    # Compute maximum absolute difference in W for convergence check
    difference <- max(abs(W - W_old))

  }

  # Obtain Theta
  Theta <- solve(W)

  # Obtain P
  P <- -cov2cor(Theta); diag(P) <- 0

  # Return all output
  return(
    list(
      W = W, Theta = Theta, P = P,
      iterations = iteration,
      converged = difference <= tol
    )
  )

}

#' @noRd
# Computes density ----
# Updated 08.10.2024
compute_density <- function(network)
{

  # Get dimensions
  dimensions <- dim(network)

  # Obtain total number of edges
  edges <- prod(dimensions)

  # Check for on-diagonal
  if(dimensions[1] == dimensions[2]){

    # Total possible edges
    total_possible <- (edges - dimensions[2]) / 2

    # Compute sparsity
    return((total_possible - (sum(network == 0, na.rm = TRUE) / 2)) / total_possible)

  }else{ # Otherwise, treat as off-diagonal

    # Compute sparsity
    return((edges - sum(network == 0, na.rm = TRUE)) / edges)

  }

}

#' @noRd
# GLASSO fit ----
# Updated 20.03.2025
glasso_fit <- function(data, S, glasso_attr, P, data_dimensions, constrain.structure, constrain.zeros, opt, ...)
{

  # Obtain structure based on community detection
  structure <- community.detection(P, ...)

  # Obtain communities
  communities <- unique_length(structure)

  # Initialize loadings
  output <- silent_call(net.loads(A = P, wc = structure, ...))
  standard_loadings <- output$std[
    glasso_attr$dimnames[[2]], seq_len(communities), drop = FALSE
  ]

  # Obtain scores
  standard_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Extract standard correlations
  standard_correlations <- cor(standard_scores, use = "pairwise")

  # Obtain model-implied correlations
  standard_R <- nload2cor(standard_loadings)
  standard_P <- cor2pcor(standard_R)

  # Get dimension names
  dimensions <- dim(standard_loadings)
  dimension_names <- dimnames(standard_loadings)

  # Obtain loadings vector and get bounds
  loadings_vector <- as.vector(standard_loadings)
  loadings_length <- length(loadings_vector)
  zeros <- swiftelse(constrain.zeros, loadings_vector != 0, rep(1, loadings_length))

  # Set up loading structure
  # Uses transpose for 2x speed up in optimization
  loading_structure <- matrix(
    FALSE, nrow = dimensions[2],
    ncol = dimensions[1],
    dimnames = list(dimension_names[[2]], dimension_names[[1]])
  )

  # Fill structure
  for(i in seq_along(structure)){
    loading_structure[structure[i], i] <- TRUE
  }

  # Optimize over loadings
  result <- try(
    egm_optimize(
      loadings_vector = loadings_vector,
      loadings_length = loadings_length,
      zeros = zeros, R = S,
      loading_structure = loading_structure,
      rows = communities, n = data_dimensions[1],
      v = data_dimensions[2],
      constrained = constrain.structure,
      lower_triangle = lower.tri(S),
      opt = opt
    ), silent = TRUE
  )

  # Check for error
  if(is(result, "try-error")){
    stop("bad result")
  }

  # Extract optimized loadings
  optimized_loadings <- matrix(
    result$par, nrow = dimensions[1],
    dimnames = dimension_names
  )

  # Update output
  output$std <- optimized_loadings

  # Extract optimized scores
  optimized_scores <- compute_scores(output, data, "network", "simple")$std.scores

  # Extract optimized correlations
  optimized_correlations <- cor(optimized_scores, use = "pairwise")

  # Obtain model-implied correlations
  optimized_R <- nload2cor(optimized_loadings)
  optimized_P <- cor2pcor(optimized_R)

  # Set up EGA object
  ega <- list(
    network = P, wc = structure, n.dim = communities,
    correlation = S, n = data_dimensions[1]
  ); class(ega) <- "EGA"

  # Set up results
  return(
    list(
      EGA = ega, structure = structure,
      model = list( # general list
        standard = list( # using standard parameters
          loadings = standard_loadings,
          scores = standard_scores,
          correlations = standard_correlations,
          fit = fit(
            n = data_dimensions[1], p = data_dimensions[2],
            R = standard_R, S = ega$correlation,
            loadings = output$std, correlations = standard_correlations,
            structure = ega$wc, ci = 0.95
          ),
          implied = list(R = standard_R, P = standard_P)
        ),
        optimized = list( # using optimized parameters
          loadings = optimized_loadings,
          scores = optimized_scores,
          correlations = optimized_correlations,
          fit = fit(
            n = data_dimensions[1], p = data_dimensions[2],
            R = optimized_R, S = ega$correlation,
            loadings = optimized_loadings,
            correlations = optimized_correlations,
            structure = ega$wc, ci = 0.95
          ),
          implied = list(R = optimized_R, P = optimized_P)
        )
      )
    )
  )

}

