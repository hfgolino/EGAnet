#' @title Bootstrap Exploratory Graph Analysis
#'
#' @description \code{bootEGA} Estimates the number of dimensions of \code{iter} bootstraps
#' using the empirical zero-order correlation matrix (\code{"parametric"}) or 
#' \code{"resampling"} from the empirical dataset (non-parametric). \code{bootEGA}
#' estimates a typical median network structure, which is formed by the median or 
#' mean pairwise (partial) correlations over the \emph{iter} bootstraps (see
#' \strong{Details} for information about the typical median network structure).
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis
#'
#' @param n Numeric (length = 1).
#' Sample size if \code{data} provided is a correlation matrix
#' 
#' @param corr Character (length = 1).
#' Method to compute correlations.
#' Defaults to \code{"auto"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"auto"} --- }
#' {Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)}
#' 
#' \item{\code{"cor_auto"} --- }
#' {Uses \code{\link[qgraph]{cor_auto}} to compute correlations. Arguments
#' can be passed along to the function}
#' 
#' \item{\code{"pearson"} --- }
#' {Pearson's correlation is computed for all variables regardless of
#' categories}
#' 
#' \item{\code{"spearman"} --- }
#' {Spearman's rank-order correlation is computed for all variables
#' regardless of categories}
#' 
#' }
#' 
#' For other similarity measures, compute them first and input them
#' into \code{data} with the sample size (\code{n})
#' 
#' @param na.data Character (length = 1).
#' How should missing data be handled?
#' Defaults to \code{"pairwise"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"} --- }
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"} --- }
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param model Character (length = 1).
#' Defaults to \code{"glasso"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"BGGM"} --- }
#' {Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{\link[BGGM]{estimate}} for more details}
#' 
#' \item{\code{"glasso"} --- }
#' {Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details}
#' 
#' \item{\code{"TMFG"} --- }
#' {Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details}
#' 
#' }
#' 
#' @param algorithm Character or 
#' \code{\link{igraph}} \code{cluster_*} function (length = 1).
#' Defaults to \code{"walktrap"}.
#' Three options are listed below but all are available
#' (see \code{\link[EGAnet]{community.detection}} for other options):
#' 
#' \itemize{
#'
#' \item{\code{"leiden"} --- }
#' {See \code{\link[igraph]{cluster_leiden}} for more details}
#' 
#' \item{\code{"louvain"} --- }
#' {By default, \code{"louvain"} will implement the Louvain algorithm using 
#' the consensus clustering method (see \code{\link[EGAnet]{community.consensus}} 
#' for more information). This function will implement
#' \code{consensus.method = "most_common"} and \code{consensus.iter = 1000} 
#' unless specified otherwise}
#' 
#' \item{\code{"walktrap"} --- }
#' {See \code{\link[igraph]{cluster_walktrap}} for more details}
#' 
#' }
#'
#' @param uni.method Character (length = 1).
#' What unidimensionality method should be used? 
#' Defaults to \code{"louvain"}.
#' Available options:
#' 
#' \itemize{
#'
#' \item{\code{"expand"} --- }
#' {Expands the correlation matrix with four variables correlated 0.50.
#' If number of dimension returns 2 or less in check, then the data 
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This method was used in the Golino et al.'s (2020)
#' \emph{Psychological Methods} simulation}
#'
#' \item{\code{"LE"} --- }
#' {Applies the Leading Eigenvector algorithm
#' (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvector solution is used; otherwise, regular EGA
#' is used. This method was used in the Christensen et al.'s (2023) 
#' \emph{Behavior Research Methods} simulation}
#' 
#' \item{\code{"louvain"} --- }
#' {Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix. If the number of dimensions is 1, 
#' then the Louvain solution is used; otherwise, regular EGA is used. 
#' This method was validated Christensen's (2022) \emph{PsyArXiv} simulation.
#' Consensus clustering can be used by specifying either
#' \code{"consensus.method"} or \code{"consensus.iter"}}
#' 
#' }
#'
#' @param iter Numeric (length = 1).
#' Number of replica samples to generate from the bootstrap analysis.
#' Defaults to \code{500} (recommended)
#'
#' @param type Character (length = 1).
#' What type of bootstrap should be performed?
#' Defaults to \code{"parametric"}.
#' Available options:
#'
#' \itemize{
#'
#' \item{\code{"parametric"} --- }
#' {Generates \code{iter} new datasets from
#' (multivariate normal random distributions) based on the
#' original dataset using \code{\link[MASS]{mvrnorm}}}
#'
#' \item{\code{"resampling"} --- }
#' {Generates \code{iter} new datasets from random subsamples 
#' of the original data}
#'
#' }
#' 
#' @param ncores Numeric (length = 1).
#' Number of cores to use in computing results.
#' Defaults to \code{ceiling(parallel::detectCores() / 2)} or half of your
#' computer's processing power.
#' Set to \code{1} to not use parallel computing
#'
#' If you're unsure how many cores your computer has,
#' then type: \code{parallel::detectCores()}
#' 
#' @param EGA.type Character (length = 1).
#' Type of EGA model to use.
#' Defaults to \code{"EGA"}
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"\link[EGAnet]{EGA}"} --- }
#' {Uses standard exploratory graph analysis}
#' 
#' \item{\code{"\link[EGAnet]{EGA.fit}"} --- }
#' {Uses \code{\link[EGAnet]{tefi}} to determine best fit of
#' \code{\link[EGAnet]{EGA}}}
#' 
#' \item{\code{"\link[EGAnet]{hierEGA}"} --- }
#' {Uses hierarhical exploratory graph analysis}
#' 
#' \item{\code{"\link[EGAnet]{riEGA}"} --- }
#' {Uses random-intercept exploratory graph analysis}
#' 
#' }
#' 
#' Arguments for \code{EGA.type} can be added (see links
#' for details on specific function arguments)
#'
#' @param typicalStructure Boolean (length = 1).
#' If \code{TRUE}, returns the median (\code{"glasso"} or \code{"BGGM"}) or
#' mean (\code{"TMFG"}) network structure and estimates its dimensions 
#' (see \strong{Details} for more information).
#' Defaults to \code{TRUE}
#'
#' @param plot.typicalStructure Boolean (length = 1).
#' If \code{TRUE}, returns a plot of the typical network structure.
#' Defaults to \code{TRUE}
#' 
#' @param seed Numeric (length = 1).
#' Defaults to \code{NULL} or random results.
#' Set for reproducible results.
#' See \href{https://github.com/hfgolino/EGAnet/wiki/Reproducibility-and-PRNG}{Reproducibility and PRNG}
#' for more details on random number generation in \code{\link{EGAnet}}
#' 
#' @param verbose Boolean (length = 1).
#' Should progress be displayed?
#' Defaults to \code{TRUE}.
#' Set to \code{FALSE} to not display progress
#'
#' @param ... Additional arguments that can be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}},
#' \code{\link[EGAnet]{EGA}},
#' \code{\link[EGAnet]{EGA.fit}},
#' \code{\link[EGAnet]{hierEGA}}, and
#' \code{\link[EGAnet]{riEGA}}
#'
#' @details The typical network structure is derived from the median (or mean) value
#' of each pairwise relationship. These values tend to reflect the
#' "typical" value taken by an edge across the bootstrap networks. Afterward,
#' the same community detection algorithm is applied to the typical network as the
#' bootstrap networks. 
#' 
#' Because the community detection algorithm is applied to the typical network structure,
#' there is a possibility that the community algorithm determines
#' a different number of dimensions than the median number derived from the bootstraps.
#' The typical network structure (and number of dimensions) may \emph{not}
#' match the empirical \code{\link[EGAnet]{EGA}} number of dimensions or
#' the median number of dimensions from the bootstrap. This result is known
#' and \emph{not} a bug. 
#'
#' @return Returns a list containing:
#'
#' \item{iter}{Number of replica samples in bootstrap}
#' 
#' \item{bootGraphs}{A list containing the networks of each replica sample}
#' 
#' \item{boot.wc}{A matrix of membership assignments for each replica network
#' with variables down the columns and replicas across the rows}
#'
#' \item{boot.ndim}{Number of dimensions identified in each replica sample}
#'
#' \item{summary.table}{A data frame containing number of replica samples, 
#' median, standard deviation, standard error, 95\% confidence intervals, and 
#' quantiles (lower = 2.5\% and upper = 97.5\%)}
#'
#' \item{frequency}{A data frame containing the proportion of times the number of dimensions was identified
#' (e.g., .85 of 1,000 = 850 times that specific number of dimensions was found)}
#'
#' \item{TEFI}{\code{\link[EGAnet]{tefi}} value for each replica sample}
#' 
#' \item{type}{Type of bootstrap used}
#'
#' \item{EGA}{Output of the empirical EGA results
#' (output will vary based on \code{EGA.type})}
#' 
#' \item{EGA.type}{Type of \code{*EGA} function used}
#'
#' \item{typicalGraph}{A list containing:
#'
#' \itemize{
#'
#' \item{\code{graph} --- }
#' {Network matrix of the median network structure}
#'
#' \item{\code{typical.dim.variables} --- }
#' {An ordered matrix of item allocation}
#'
#' \item{\code{wc} --- }
#' {Membership assignments of the median network}
#' 
#' }
#' 
#' }
#' 
#' \item{plot.typical.ega}{Plot output if \code{plot.typicalStructure = TRUE}}
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' \dontrun{
#' # Standard EGA parametric example
#' boot.wmt <- bootEGA(
#'   data = wmt, iter = 500,
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Standard resampling example
#' boot.wmt <- bootEGA(
#'   data = wmt, iter = 500,
#'   type = "resampling", ncores = 2
#' )
#' 
#' # Example using {igraph} `cluster_*` function
#' boot.wmt.spinglass <- bootEGA(
#'   data = wmt, iter = 500,
#'   algorithm = igraph::cluster_spinglass,
#'   # use any function from {igraph}
#'   type = "parametric", ncores = 2
#' )
#'
#' # EGA fit example
#' boot.wmt.fit <- bootEGA(
#'   data = wmt, iter = 500,
#'   EGA.type = "EGA.fit",
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Hierarchical EGA example
#' boot.wmt.hier <- bootEGA(
#'   data = wmt, iter = 500,
#'   EGA.type = "hierEGA",
#'   type = "parametric", ncores = 2
#' )
#' 
#' # Random-intercept EGA example
#' boot.wmt.ri <- bootEGA(
#'   data = wmt, iter = 500,
#'   EGA.type = "riEGA",
#'   type = "parametric", ncores = 2
#' )}
#'
#' @references
#' \strong{Original implementation of bootEGA} \cr
#' Christensen, A. P., & Golino, H. (2021).
#' Estimating the stability of the number of factors via Bootstrap Exploratory Graph Analysis: A tutorial.
#' \emph{Psych}, \emph{3}(3), 479-500.
#'
#' @seealso \code{\link[EGAnet]{itemStability}} to estimate the stability of
#' the variables in the empirical dimensions and
#' \code{\link[EGAnet]{dimensionStability}} to estimate the stability of
#' the dimensions (structural consistency)
#'
#' @export
#'
# Bootstrap EGA ----
# Updated 04.09.2023
bootEGA <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("BGGM", "glasso", "TMFG"),  
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    iter = 500, type = c("parametric", "resampling"),
    ncores, EGA.type = c("EGA", "EGA.fit", "hierEGA", "riEGA"),
    typicalStructure = TRUE, plot.typicalStructure = TRUE,
    seed = NULL, verbose = TRUE, 
    ...
) 
{
  
  # Store random state (if there is one)
  store_state()
  
  # Get ellipse arguments
  ellipse <- list(needs_usable = FALSE, ...)
  
  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", bootEGA)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)
  type <- set_default(type, "parametric", bootEGA)
  EGA.type <- set_default(EGA.type, "ega", bootEGA)
  
  # Set cores
  if(missing(ncores)){ncores <- ceiling(parallel::detectCores() / 2)}
  
  # Argument errors (return data in case of tibble)
  data <- bootEGA_errors(
    data, n, iter, ncores, typicalStructure,
    plot.typicalStructure, seed, verbose, ...
  )
  
  # `EGA.estimate` will handle legacy arguments and data processing 
  
  # Ensure matrix
  data <- as.matrix(data)
  
  # Ensure variable names
  data <- ensure_dimension_names(data)
  
  # Get data dimensions
  dimensions <- dim(data)

  # Get variable names
  variable_names <- dimnames(data)[[2]]

  # Obtain EGA function
  ega_function <- switch(
    EGA.type,
    "ega" = EGA,
    "ega.fit" = EGA.fit,
    "riega" = riEGA,
    "hierega" = hierEGA,
    stop(
      paste0(
        "EGA.type = \"", EGA.type, "\" is not supported. Please use one of the default options: ", 
        paste0("\"", as.character(formals(bootEGA)$EGA.type)[-1], "\"", collapse = ", ")
      ), call. = FALSE
    )
  )
  
  # Set up default EGA arguments
  EGA_ARGS <- list(
    data = data, n = n, corr = corr,
    na.data = na.data, model = model,
    algorithm = algorithm, uni.method = uni.method,
    plot.EGA = FALSE, verbose = FALSE
  )
  
  # Set flag for hierarchical EGA
  hierarchical <- EGA.type == "hierega"
  
  # Branch for hierarchical
  if(hierarchical){
    
    # Handle hierarchical EGA arguments
    ellipse <- handle_hierEGA_ARGS(algorithm, ellipse, names(ellipse))
    
    # For `hierEGA`, remove "algorithm" from EGA arguments
    EGA_ARGS <- EGA_ARGS[names(EGA_ARGS) != "algorithm"]
    
  }
  
  # Estimate empirical EGA
  empirical_EGA <- do.call(
    what = ega_function,
    args = c(EGA_ARGS, ellipse)
  )
  # If "n" is NULL and a correlation matrix is supplied,
  # then an error will be thrown within `EGA.estimate`,
  # so no need to check for "n"
  
  # Obtain EGA output
  empirical_EGA_output <- get_EGA_object(empirical_EGA)
  
  # Check for hierarchical EGA (if so, get lowest level)
  if(hierarchical){
    empirical_EGA_output <- empirical_EGA_output$lower_order
  }
  
  # Check for seed
  if(!is.null(seed)){
    seeds <- reproducible_seeds(iter, seed)
  }else{ 
    
    # Set all seeds to zero (or random)
    seeds <- rep(0, iter)
    
    # Check for external suppression (from `invariance`)
    if(!"suppress" %in% names(ellipse) || !ellipse$suppress){
      message("Argument 'seed' is set to `NULL`. Results will not be reproducible. Set 'seed' for reproducible results")
    }

  }
  
  # Check for parametric (pre-compute values)
  if(type == "parametric"){
    
    # Get parameters
    mvrnorm_parameters <- mvrnorm_precompute(
      cases = empirical_EGA_output$n,
      Sigma = empirical_EGA_output$correlation
    )
    
    # Set case sequence to be NULL
    case_sequence <- NULL
    
  }else if(type == "resampling"){
    
    # Get case sequence
    case_sequence <- seq_len(empirical_EGA_output$n)
    
    # Set parameters to NULL
    mvrnorm_parameters <- NULL
    
  }
  
  # Perform bootstrap using parallel processing
  boots <- parallel_process(
    # Parallel processing arguments
    iterations = iter,
    datalist = seeds, # data is generated in each iteration
    ncores = ncores, progress = verbose,
    clear = swiftelse("clear" %in% names(ellipse), ellipse$clear, FALSE), # from `invariance`
    # Standard EGA arguments
    FUN = function(
      seed_value, data, type, case_sequence, 
      mvrnorm_parameters, EGA_ARGS, ellipse
    ){
      
      # Replace data in EGA arguments
      EGA_ARGS$data <- reproducible_bootstrap(
        seed = seed_value, data = data,
        case_sequence = case_sequence,
        mvrnorm_parameters = mvrnorm_parameters,
        type = type
      )
      
      # Estimate EGA
      return(
        empirical_EGA <- do.call(
          what = ega_function,
          args = c(EGA_ARGS, ellipse)
        )
      )

    }, # Send all argument variables
    data, type, case_sequence, 
    mvrnorm_parameters, EGA_ARGS, ellipse
  )
  
  # Obtain bootstrap EGA output
  bootstrap_EGA_output <- lapply(boots, get_EGA_object)
  
  # Branch based on hierarchical EGA
  if(hierarchical){
    
    # Prepare results
    results <- list(
      lower_order = prepare_bootEGA_results(
        lapply(bootstrap_EGA_output, function(x){x$lower_order}),
        iter
      ),
      # To avoid issues with differing communities,
      # re-value each variable's membership to their
      # higher order community
      higher_order = prepare_bootEGA_results(
        revalue_memberships(bootstrap_EGA_output),
        iter
      )
    )
    
  }else{
    
    # Get results
    results <- prepare_bootEGA_results(bootstrap_EGA_output, iter)
    
  }
  
  # Add additional results
  results[c("type", "iter", "EGA", "EGA.type")] <- list(
    type, iter, empirical_EGA, EGA.type
  ) # `iter` is redundant except for `hierEGA`
  
  # No attributes needed (all information is contained in output)
  
  # Set class
  class(results) <- "bootEGA"
  
  # Check for typical structure results
  if(typicalStructure){
    
    # Obtain results
    results$typicalGraph <- do.call(
      what = estimate_typicalStructure,
      args = c(
        list(data = data, results = results, verbose = verbose),
        ellipse
      )
    )

    # Check for plot
    if(plot.typicalStructure){
      
      # Get plot
      results$plot.typical.ega <- plot(results, ...)
      
      # Actually send plot
      silent_plot(results$plot.typical.ega)
      
    }
    
  }
  
  # Restore random state (if there is one)
  restore_state()
  
  # Return result
  return(results)
  
}

# Bug checking ----
# data = NetworkToolbox::neoOpen; n = NULL; corr = "auto"; na.data = "pairwise"
# model = "glasso"; algorithm = "walktrap"; uni.method = "louvain"
# iter = 100; type = "parametric"; ncores = 8; EGA.type = "EGA"
# typicalStructure = TRUE; plot.typicalStructure = FALSE;
# seed = 1234; verbose = TRUE

#' @noRd
# Errors ----
# Updated 07.09.2023
bootEGA_errors <- function(
    data, n, iter, ncores, typicalStructure,
    plot.typicalStructure, seed, verbose, ...
)
{
  
  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "bootEGA")
  
  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }
  
  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "bootEGA")
    typeof_error(n, "numeric", "bootEGA")
  }
  
  # 'iter' errors
  length_error(iter, 1, "bootEGA")
  typeof_error(iter, "numeric", "bootEGA")
  range_error(iter, c(1, Inf), "bootEGA")
  
  # 'ncores' errors
  length_error(ncores, 1, "bootEGA")
  typeof_error(ncores, "numeric", "bootEGA")
  range_error(ncores, c(1, parallel::detectCores()), "bootEGA")
   
  # 'typicalStructure' errors
  length_error(typicalStructure, 1, "bootEGA")
  typeof_error(typicalStructure, "logical", "bootEGA")
  
  # 'plot.typicalStructure' errors
  length_error(plot.typicalStructure, 1, "bootEGA")
  typeof_error(plot.typicalStructure, "logical", "bootEGA")
  
  # 'seed' errors
  if(!is.null(seed)){
    length_error(seed, 1, "bootEGA")
    typeof_error(seed, "numeric", "bootEGA")
    range_error(seed,  c(0, Inf), "bootEGA")
  }
  
  # 'verbose' errors
  length_error(verbose, 1, "bootEGA")
  typeof_error(verbose, "logical", "bootEGA")
  
  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }
  
  # Return usable data (in case of tibble)
  return(data)
  
}

#' @exportS3Method 
# S3 Print Method ----
# Updated 21.07.2023
print.bootEGA <- function(x, ...)
{
  
  # Ensure proper EGA object
  ega_object <- get_EGA_object(x)
  
  # Branch for hierarchical EGA
  if(is(ega_object, "hierEGA")){
    
    # Set proper EGA type name
    ega_type <- switch(
      x$EGA.type,
      "ega" = "EGA",
      "ega.fit" = "EGA.fit",
      "hierega" = "hierEGA",
      "riega" = "riEGA"
    )
    
    # Print EGA type
    cat(paste0("EGA Type: ", ega_type), "\n")
    
    # Set up methods
    cat(
      paste0(
        "Bootstrap Samples: ",
        x$iter, " (", totitle(x$type), ")"
      )
    )
    
    # Add breakspace
    cat("\n\n------------\n\n")
    
    # Print level
    cat(
      styletext(
        text = styletext(
          text =  "Lower Order\n\n", 
          defaults = "underline"
        ),
        defaults = "bold"
      )
    )
    
    # Print network information
    send_network_methods(ega_object$lower_order$network, boot = TRUE)
    
    # Add line break
    cat("\n")
    
    # Print community detection
    print(ega_object$lower_order$wc, boot = TRUE)
    
    # Add line break
    cat("\n")
    
    # Get unidimensional attributes
    unidimensional_attributes <- attr(ega_object$lower_order, "unidimensional")
    
    # Obtain unidimensional method
    unidimensional_method <- switch(
      unidimensional_attributes$uni.method,
      "expand" = "Expand",
      "le" = "Leading Eigenvector",
      "louvain" = "Louvain"
    )
    
    # Set up unidimensional print
    if(
      unidimensional_method == "Louvain" &
      "consensus.iter" %in% names(unidimensional_attributes$consensus)
    ){
      
      # Set up consensus attributes
      consensus_attributes <- unidimensional_attributes$consensus
      
      # Obtain consensus name
      consensus_name <- switch(
        consensus_attributes$consensus.method,
        "highest_modularity" = "Highest Modularity",
        "iterative" = "Iterative",
        "most_common" = "Most Common",
        "lowest_tefi" = "Lowest TEFI"
      )
      
      # Update unidimensional method text
      unidimensional_method <- paste0(
        unidimensional_method, " (", consensus_name,
        " for ", consensus_attributes$consensus.iter,
        " iterations)"
      )
      
    }
    
    # Print unidimensional
    cat("Unidimensional Method: ", unidimensional_method)
    
    # Add break space
    cat("\n\n----\n")
    
    # Print frequency table
    frequency_df <- as.data.frame(
      do.call(rbind, lapply(x$lower_order$frequency, as.character))
    )

    # Adjust dimension names (`dimnames` doesn't work)
    colnames(frequency_df) <- NULL
    row.names(frequency_df) <- c("", "Frequency: ")
    # Finally, print
    print(frequency_df)
    
    # Print summary table
    cat(
      paste0(
        "\nMedian dimensions: ", x$lower_order$summary.table$median.dim,
        " [", round(x$lower_order$summary.table$Lower.CI, 2), ", ",
        round(x$lower_order$summary.table$Upper.CI, 2), "] 95% CI"
      )
    )
    
    # Add breakspace
    cat("\n\n------------\n\n")
    
    # Print level
    cat(
      styletext(
        text = styletext(
          text =  "Higher Order\n\n", 
          defaults = "underline"
        ),
        defaults = "bold"
      )
    )
    
    # Print community detection
    print(ega_object$higher_order$wc, boot = TRUE)

    # Add break space
    cat("\n\n----\n")
    
    # Print frequency table
    frequency_df <- as.data.frame(
      do.call(rbind, lapply(x$higher_order$frequency, as.character))
    )
    
    # Adjust dimension names (`dimnames` doesn't work)
    colnames(frequency_df) <- NULL
    row.names(frequency_df) <- c("", "Frequency: ")
    # Finally, print
    print(frequency_df)
    
    # Print summary table
    cat(
      paste0(
        "\nMedian dimensions: ", x$higher_order$summary.table$median.dim,
        " [", round(x$higher_order$summary.table$Lower.CI, 2), ", ",
        round(x$higher_order$summary.table$Upper.CI, 2), "] 95% CI"
      )
    )
    
  }else{
    
    # Print network information
    send_network_methods(ega_object$network, boot = TRUE)
    
    # Add line break
    cat("\n")
    
    # Print community detection
    print(ega_object$wc, boot = TRUE)
    
    # Add line break
    cat("\n")
    
    # Do not print unidimensional for `EGA.fit`
    if(x$EGA.type != "ega.fit"){
      
      # Get unidimensional attributes
      unidimensional_attributes <- attr(ega_object, "unidimensional")
      
      # Obtain unidimensional method
      unidimensional_method <- switch(
        unidimensional_attributes$uni.method,
        "expand" = "Expand",
        "le" = "Leading Eigenvector",
        "louvain" = "Louvain"
      )
      
      # Set up unidimensional print
      if(
        unidimensional_method == "Louvain" &
        "consensus.iter" %in% names(unidimensional_attributes$consensus)
      ){
        
        # Set up consensus attributes
        consensus_attributes <- unidimensional_attributes$consensus
        
        # Obtain consensus name
        consensus_name <- switch(
          consensus_attributes$consensus.method,
          "highest_modularity" = "Highest Modularity",
          "iterative" = "Iterative",
          "most_common" = "Most Common",
          "lowest_tefi" = "Lowest TEFI"
        )
        
        # Update unidimensional method text
        unidimensional_method <- paste0(
          unidimensional_method, " (", consensus_name,
          " for ", consensus_attributes$consensus.iter,
          " iterations)"
        )
        
      }
      
      # Print unidimensional
      cat("Unidimensional Method: ", unidimensional_method)
      
    }
    
    # Add break space
    cat("\n\n----\n\n")
    
    # Set proper EGA type name
    ega_type <- switch(
      x$EGA.type,
      "ega" = "EGA",
      "ega.fit" = "EGA.fit",
      "hierega" = "hierEGA",
      "riega" = "riEGA"
    )
    
    # Print EGA type
    cat(paste0("EGA Type: ", ega_type), "\n")
    
    # Set up methods
    cat(
      paste0(
        "Bootstrap Samples: ",
        x$iter, " (", totitle(x$type), ")\n"
      )
    )
    
    # Print frequency table
    frequency_df <- as.data.frame(
      do.call(rbind, lapply(x$frequency, as.character))
    )

    # Adjust dimension names (`dimnames` doesn't work)
    colnames(frequency_df) <- NULL
    row.names(frequency_df) <- c("", "Frequency: ")
    # Finally, print
    print(frequency_df)
    
    # Print summary table
    cat(
      paste0(
        "\nMedian dimensions: ", x$summary.table$median.dim,
        " [", round(x$summary.table$Lower.CI, 2), ", ",
        round(x$summary.table$Upper.CI, 2), "] 95% CI"
      )
    )
    
  }

}

#' @exportS3Method 
# S3 Summary Method ----
# Updated 05.07.2023
summary.bootEGA <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method 
# S3 Plot Method ----
# Updated 21.07.2023
plot.bootEGA <- function(x, ...)
{
  
  # Check for hierarchical EGA
  if(x$EGA.type == "hierega"){
    
    # Set up results to be similar to `hierEGA`
    ## Lower order
    lower_order_result <- list(
      network = x$typicalGraph$lower_order$graph,
      wc = x$typicalGraph$lower_order$wc
    )
    class(lower_order_result) <- "EGA"
    
    ## Higher order
    higher_order_result <- list(
      network = x$typicalGraph$higher_order$graph,
      wc = x$typicalGraph$higher_order$wc
    )
    class(higher_order_result) <- "EGA"
    
    # Set up results
    results <- list(
      lower_order = lower_order_result,
      higher_order = higher_order_result,
      parameters = x$typicalGraph$parameters
    )
    
    # Transfer "methods" attributes
    attr(results, "methods") <- attr(x$EGA, "methods")
    
    # Set class
    class(results) <- "hierEGA"
    
    # Return plot as hierarchical EGA
    return(plot(results, ...))
    
  }else{
    
    # Return plot
    return(
      single_plot(
        network = x$typicalGraph$graph,
        wc = x$typicalGraph$wc,
        ...
      ) 
    )
    
  }
  
}

#' @noRd
# Handle `hierEGA` arguments ----
# Updated 31.07.2023
handle_hierEGA_ARGS <- function(algorithm, ellipse, ellipse_names)
{
  
  # Determine whether "lower" and "higher" algorithms are used
  if("lower.algorithm" %in% names(ellipse) & !"higher.algorithm" %in% ellipse_names){
    
    # Use "algorithm" as "higher.algorithm"
    ellipse$higher.algorithm <- algorithm
    
    # Send message
    warning(
      paste0(
        "Argument 'lower.algorithm' was set to \"", ellipse$lower.algorithm,
        "\" but 'higher.algorithm' was not set ",
        "with 'EGA.type = \"hierEGA\"'.", 
        "\n\nSetting 'higher.algorithm = \"", algorithm, "\"'"
      ), call. = FALSE
    )
    
  }else if(!"lower.algorithm" %in% names(ellipse) & "higher.algorithm" %in% ellipse_names){
    
    # Set default "louvain"
    ellipse$lower.algorithm <- "louvain"
    ellipse$consensus.method <- "most_common"
    ellipse$consensus.iter <- 1000
    
    # Send message
    warning(
      paste0(
        "Argument 'higher.algorithm' was set to \"", ellipse$higher.algorithm,
        "\" but 'lower.algorithm' was not set ",
        "with 'EGA.type = \"hierEGA\"'.",
        "\n\nSetting 'lower.algorithm' to the default: ",
        "\"louvain\" with most common consensus clustering (1000 iterations)"
      ), call. = FALSE
    )
    
  }else if(!any(c("lower.algorithm", "higher.algorithm") %in% ellipse_names)){
    
    # Set to defaults but use "algorithm" as higher order algorithm
    ellipse$lower.algorithm <- "louvain"
    ellipse$consensus.method <- "most_common"
    ellipse$consensus.iter <- 1000
    ellipse$higher.algorithm <- algorithm
    
  }
  
  # Return ellipse
  return(ellipse)

}

#' @noRd
# Revalue single membership ----
# Also used in `hierEGA`
# Updated 22.07.2023
single_revalue_memberships <- function(lower, higher)
{
  
  # Assign new memberships
  lower[] <- higher[lower]
  
  # Return lower
  return(lower)
  
}

#' @noRd
# Revalue higher order results ----
# Updated 22.07.2023
revalue_memberships <- function(bootstrap_EGA_output)
{
  
  # Return revalued memberships
  return(
    
    # Loop over iterations
    lapply(bootstrap_EGA_output, function(output){

      # Re-assign to the higher order output
      output$higher_order$wc <- single_revalue_memberships(
        output$lower_order$wc, output$higher_order$wc
      )
      
      # Return higher order output
      return(output$higher_order)
      
    })
    
  )

}

#' @noRd
# Prepare `bootEGA` results ----
# Self-contained to work on `EGA` bootstraps
# Updated 31.07.2023
prepare_bootEGA_results <- function(boot_object, iter)
{
  
  # Get networks
  boot_networks <- lapply(boot_object, function(x){x$network})
  
  # Get memberships
  boot_memberships <- t(
    nvapply(boot_object, function(x){x$wc}, LENGTH = length(boot_object[[1]]$wc))
  )
  
  # Get bootstrap dimensions
  boot_n.dim <- nvapply(boot_object, function(x){x$n.dim})
  
  # Set up data frame
  boot_dimensions <- fast.data.frame(
    c(seq_len(iter), boot_n.dim),
    nrow = iter, ncol = 2,
    colnames = c("Boot.Number", "N.Dim")
  )
  
  # Pre-compute standard deviation and confidence interval
  median_boot <- median(boot_n.dim, na.rm = TRUE)
  se_boot <- sd(boot_n.dim, na.rm = TRUE)
  ci_boot <- se_boot * qt(0.95 / 2 + 0.5, iter - 1)
  
  # Set up descriptive statistics
  summary_table <- fast.data.frame(
    c(
      iter, median_boot, se_boot, ci_boot, 
      median_boot - ci_boot, median_boot + ci_boot,
      quantile(boot_n.dim, c(0.025, 0.975), na.rm = TRUE)
    ), ncol = 8,
    colnames = c(
      "n.Boots", "median.dim", "SE.dim", "CI.dim",
      "Lower.CI", "Upper.CI", "Lower.Quantile", "Upper.Quantile"
    ), row.names = ""
  )
  
  # Compute frequency
  frequencies <- count_table(boot_n.dim, proportion = TRUE)
  dimnames(frequencies)[[2]] <- c("# of Factors", "Frequency")
  
  # Return results (keep legacy names)
  return(
    list(
      iter = iter,
      bootGraphs = boot_networks,
      boot.wc = boot_memberships,
      boot.ndim = boot_n.dim,
      summary.table = summary_table,
      frequency = frequencies[
        order(frequencies[, "# of Factors"]),
      ],
      TEFI = nvapply(boot_object, function(x){x$TEFI})
    )
  )
  
}

#' @noRd
# Typical Walktrap `EGA.fit` ----
# Updated 07.07.2023
typical_walktrap_fit <- function(network, dimensions, ellipse)
{
  
  # Check for parameter search space
  if(!"steps" %in% names(ellipse)){
    steps <- 3:8 # default
  }else{
    
    # Set steps
    steps <- ellipse$steps
    
    # Remove objective function from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "steps"]
    
  }
  
  # Get the length of the steps
  step_length <- length(steps)
  
  # Obtain search matrix
  search_matrix <- t(nvapply(
    steps, function(step){
      do.call(
        what = community.detection,
        args = c(
          list( # Necessary call
            network = network,
            algorithm = "walktrap",
            steps = step
          ),
          ellipse # pass on ellipse
        )
      )
    }, LENGTH = dimensions[2]
  ))
  
  # Format names
  dimnames(search_matrix) <- list(
    steps, dimnames(network)[[2]]
  )
  
  # Return results
  return(
    list(
      search_matrix = search_matrix,
      parameters = steps
    )
  )
  
}

#' @noRd
# Typical Leiden `EGA.fit` ----
# Updated 07.08.2023
typical_leiden_fit <- function(network, dimensions, ellipse)
{
  
  # Check for objective function
  if(!"objective_function" %in% names(ellipse)){
    
    # Default for {EGAnet}
    objective_function <- "modularity"
    
    # Send warning
    warning(
      paste0(
        "{EGAnet} uses \"modularity\" as the default objective function in the Leiden algorithm. ",
        "In contrast, {igraph} uses \"CPM\". Set `objective_function = \"CPM\"` to use the Constant Potts ",
        "Model in {EGAnet}"
      ), call. = FALSE
    )
    
  }else{
    
    # Set objective function
    objective_function <- ellipse$objective_function
    
    # Remove objective function from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "objective_function"]
    
  }
  
  # Check for parameter search space
  if(!"resolution_parameter" %in% names(ellipse)){
    
    # Switch based on objective function
    if(objective_function == "CPM"){
      resolution_parameter <- seq.int(0, max(abs(network)), 0.01) # default 
    }else if(objective_function == "modularity"){
      resolution_parameter <- seq.int(0, 2, 0.05) # default
    }
    
  }else{
    
    # Set resolution parameter
    resolution_parameter <- ellipse$resolution_parameter
    
    # Remove resolution parameter from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "resolution_parameter"]
    
  }
  
  # Obtain search matrix
  search_matrix <- t(nvapply(
    resolution_parameter, function(resolution_parameter){
      do.call(
        what = community.detection,
        args = c(
          list( # Necessary call
            network = network,
            algorithm = "leiden",
            resolution_parameter = resolution_parameter,
            objective_function = objective_function
          ),
          ellipse # pass on ellipse
        )
      )
    }, LENGTH = dimensions[2]
  ))
  
  # Format names
  dimnames(search_matrix) <- list(
    resolution_parameter, dimnames(network)[[2]]
  )
  
  # Return results
  return(
    list(
      search_matrix = search_matrix,
      parameters = resolution_parameter
    )
  )
  
}

#' @noRd
# Typical Louvain `EGA.fit` ----
# Updated 07.08.2023
typical_louvain_fit <- function(network, dimensions, ellipse)
{
  
  # Check for parameter search space
  if(!"resolution_parameter" %in% names(ellipse)){
    resolution_parameter <- seq.int(0, 2, 0.05) # default
  }else{
    
    # Set resolution parameter
    resolution_parameter <- ellipse$resolution_parameter
    
    # Remove resolution parameter from ellipse to
    # make other arguments available in `do.call`
    ellipse <- ellipse[names(ellipse) != "resolution_parameter"]
    
  }
  
  # Obtain search matrix
  search_matrix <- t(nvapply(
    resolution_parameter, function(resolution_parameter){
      do.call(
        what = community.consensus,
        args = c(
          list( # Necessary call
            network = network,
            resolution = resolution_parameter
          ),
          ellipse # pass on ellipse
        )
      )
    }, LENGTH = dimensions[2]
  ))
  
  # Format names
  dimnames(search_matrix) <- list(
    resolution_parameter, dimnames(network)[[2]]
  )
  
  # Return results
  return(
    list(
      search_matrix = search_matrix,
      parameters = resolution_parameter
    )
  )
  
}

#' @noRd
# Typical `EGA.fit` memberships ----
# Needs to be handled consistently
# Updated 07.08.2023
estimate_typical_EGA.fit <- function(results, ellipse)
{
  
  # Get proper EGA object
  ega_object <- get_EGA_object(results)
  
  # Attributes from empirical EGA
  model_attributes <- attr(ega_object$network, "methods")
  algorithm_attributes <- attr(ega_object$wc, "methods")
  
  # Set model, algorithm and unidimensional method
  model <- tolower(model_attributes$model)
  algorithm <- tolower(algorithm_attributes$algorithm)
  
  # Check for {BGGM}
  if(model == "bggm"){
    stop("Due to CRAN check issues, `model = \"BGGM\"` is not available at the moment.")
  }
  
  # Get network
  network <- switch(
    model,
    # "bggm" = symmetric_matrix_lapply(results$bootGraphs, median),
    "glasso" = symmetric_matrix_lapply(results$bootGraphs, median),
    "tmfg" = symmetric_matrix_lapply(results$bootGraphs, mean)
  )
  
  # Make sure proper names are there
  dimnames(network) <- dimnames(ega_object$network)
  
  # Get network dimensions
  dimensions <- dim(network)
  
  # Branch for fit function
  fit_FUN <- switch(
    algorithm,
    "walktrap" = typical_walktrap_fit,
    "leiden" = typical_leiden_fit,
    "louvain" = typical_louvain_fit
  )
  
  # Get fit result
  fit_result <- fit_FUN(network, dimensions, ellipse)
  
  # Set objective function (will be NULL for Louvain and Walktrap)
  objective_function <- fit_result$objective_function
  
  # Obtain only unique solutions
  search_unique <- unique(fit_result$search_matrix, MARGIN = 1)
  
  # Determine best fitting solution
  fit_values <- apply(
    search_unique, 1, function(membership){
      tefi(ega_object$correlation, membership, verbose = FALSE)
    }$VN.Entropy.Fit
  )
  
  # Determine best fit index
  best_index <- which.min(fit_values)
  
  # Obtain best solution
  best_solution <- search_unique[best_index,]
  
  # Add methods to membership attributes
  attr(best_solution, "methods") <- list(
    algorithm = obtain_algorithm_name(algorithm),
    # `obtain_algorithm_name` is with `community.detection`
    objective_function = objective_function
  )
  
  # Set class (proper `EGA.estimate` printing)
  class(best_solution) <- "EGA.community"
  
  # Set up dimension variables data frame
  ## Mainly for legacy, redundant with named `wc`
  dim.variables <- fast.data.frame(
    c(dimnames(ega_object$network)[[2]], best_solution),
    nrow = dimensions[2], ncol = 2,
    colnames = c("items", "dimension")
  )
  
  # Return results
  return(
    list(
      graph = network,
      typical.dim.variables = dim.variables[order(dim.variables$dimension),],
      wc = best_solution, n.dim = unique_length(best_solution),
      EntropyFit = fit_values,
      Lowest.EntropyFit = fit_values[best_index],
      parameter.space = fit_result$parameters
    )
  )
  
}

#' @noRd
# Typical network and memberships ----
# Updated 24.10.2023
estimate_typicalStructure <- function(
    data, results, verbose, ...
)
{

  # Get ellipse arguments
  ellipse <- list(...)
  
  # If results are from `EGA.fit` or `hierEGA`, handle separately
  if(results$EGA.type == "ega.fit"){
    return(estimate_typical_EGA.fit(results, ellipse))
  }else if(results$EGA.type == "hierega"){
    # Get typical structure for lower order first, then higher
    ega_object <- get_EGA_object(results)$lower_order
  }else{ # Get proper EGA object
    ega_object <- get_EGA_object(results)
  }
  
  # Attributes from empirical EGA
  model_attributes <- attr(ega_object$network, "methods")
  algorithm_attributes <- attr(ega_object$wc, "methods")
  unidimensional_attributes <- attr(ega_object, "unidimensional")
  
  # Set model, algorithm and unidimensional method
  model <- tolower(model_attributes$model)
  algorithm <- tolower(algorithm_attributes$algorithm)
  uni.method <- tolower(unidimensional_attributes$uni.method)
  
  # Check for {BGGM}
  if(model == "bggm"){
    stop("Due to CRAN check issues, `model = \"BGGM\"` is not available at the moment.")
  }
  
  # Branch for hierarchical EGA
  if(results$EGA.type == "hierega"){
    
    # Get network
    network <- switch(
      model,
      # "bggm" = symmetric_matrix_lapply(results$lower_order$bootGraphs, median),
      "glasso" = symmetric_matrix_lapply(results$lower_order$bootGraphs, median),
      "tmfg" = symmetric_matrix_lapply(results$lower_order$bootGraphs, mean)
    )
    
  }else{
    
    # Get network
    network <- switch(
      model,
      # "bggm" = symmetric_matrix_lapply(results$bootGraphs, median),
      "glasso" = symmetric_matrix_lapply(results$bootGraphs, median),
      "tmfg" = symmetric_matrix_lapply(results$bootGraphs, mean)
    )
    
  }

  # Make sure proper names are there
  dimnames(network) <- dimnames(ega_object$network)
  
  # Check for function or non-Louvain method
  if(is.function(algorithm) || algorithm != "louvain"){
    
    # Apply non-Louvain method
    wc <- do.call(
      what = community.detection,
      args = c(
        list(
          network = network, algorithm = algorithm,
          membership.only = TRUE
        ),
        ellipse # pass on ellipse
      )
    )
    
  }else{ # for Louvain, use consensus clustering
    
    # Check for consensus method
    if(!"consensus.method" %in% names(ellipse)){
      ellipse$consensus.method <- "most_common" # default
    }
    
    # Check for consensus iterations
    if(!"consensus.iter" %in% names(ellipse)){
      ellipse$consensus.iter <- 1000 # default
    }
    
    # Force lower order for hierarchical EGA
    if(results$EGA.type == "hierega"){
      ellipse$order <- "lower"
    }

    # Apply consensus clustering
    wc <- do.call(
      what = community.consensus,
      args = c(
        list(
          network = network,
          membership.only = TRUE
        ),
        ellipse # pass on ellipse
      )
    )
    
  }

  # Obtain arguments for model
  model_ARGS <- switch(
    model,
    # "bggm" = c(
    #   obtain_arguments(BGGM::estimate, model_attributes),
    #   overwrite_arguments(
    #     # defaults for `BGGM:::select.estimate`
    #     list(cred = 0.95, alternative = "two.sided"), model_attributes
    #   )
    # ),
    "glasso" = obtain_arguments(EBICglasso.qgraph, model_attributes),
    "tmfg" = obtain_arguments(TMFG, model_attributes)
  )
  
  # Make adjustments for each model (removes extraneous arguments)
  model_ARGS <- adjust_model_arguments(model_ARGS)
  # Found in `EGA` internals
  
  # Set up arguments for unidimensional
  unidimensional_ARGS <- list( # standard arguments
    data = data, n = ega_object$n, 
    corr = model_attributes$corr, 
    na.data = model_attributes$na.data,
    model = model, uni.method = uni.method,
    verbose = FALSE
  )
  
  # `data` at this point will be data or correlation matrix
  # For non-BGGM network estimation, OK to use correlation matrix
  if(model_attributes$model != "bggm"){
    unidimensional_ARGS$data <- ega_object$correlation
  }
  
  # Additional arguments for model
  unidimensional_ARGS <- c(unidimensional_ARGS, model_ARGS)
  
  # Third, obtain the unidimensional result
  unidimensional_result <- do.call(
    what = community.unidimensional,
    args = unidimensional_ARGS
  )
  
  # Determine result
  if(unique_length(unidimensional_result) == 1){
    wc <- unidimensional_result
  }
  
  # Set up dimension variables data frame
  ## Mainly for legacy, redundant with named `wc`
  dim.variables <- fast.data.frame(
    c(dimnames(ega_object$network)[[2]], wc),
    nrow = length(wc), ncol = 2,
    colnames = c("items", "dimension")
  )
  
  # Branch for hierarchical EGA
  if(results$EGA.type == "hierega"){
    
    # Get hierarchical methods
    hierarchical_methods <- attr(results$EGA, "methods")
    
    # Check for scores
    if(hierarchical_methods$scores == "factor"){
      
      # Send warning
      warning(
        "Typical network structure for `hierEGA` is not supported for `scores = \"factor\"",
        call. = FALSE
      )
      
    }else if(hierarchical_methods$scores == "network"){
      
      # Compute network scores
      network_output <- net.scores(
        data = data, A = network, wc = wc,
        rotation = hierarchical_methods$rotation,
        loading.method = hierarchical_methods$loading.method,
        scoring.method = "network",
        ...
      )
      
      # Score estimates
      if(is.null(hierarchical_methods$rotation)){
        score_estimates <- network_output$scores$std.scores
        lower_loadings <- network_output$loadings$std
      }else{
        score_estimates <- network_output$scores$rot.scores
        lower_loadings <- network_output$loadings$rotated
      }
      
    }
    
    # Store higher order results
    higher_order <- EGA(
      data = score_estimates, corr = model_attributes$corr, 
      na.data = model_attributes$na.data,
      model = model, algorithm = ellipse$higher.algorithm,
      uni.method = unidimensional_attributes$uni.method,
      plot.EGA = FALSE, verbose = verbose, ...
    )
    
    # Return results
    return(
      list(
        lower_order = list(
          graph = network,
          typical.dim.variables = dim.variables,
          wc = wc, n.dim = unique_length(wc)
        ),
        higher_order = list(
          graph = higher_order$network,
          typical.dim.variables = higher_order$dim.variables,
          wc = higher_order$wc,
          n.dim = higher_order$n.dim
        ),
        parameters = list(
          lower_loadings = lower_loadings,
          lower_scores = score_estimates
        ) 
      )
    )
    
  }else{ 
    
    # Return results
    return(
      list(
        graph = network,
        typical.dim.variables = dim.variables[order(dim.variables$dimension),],
        wc = wc, n.dim = unique_length(wc)
      )
    )
    
  }

}

