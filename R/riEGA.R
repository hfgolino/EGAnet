#' @title Random-Intercept \code{\link[EGAnet]{EGA}}
#'
#' @description Estimates the number of substantive dimensions after controlling
#' for wording effects. EGA is applied to a residual correlation matrix after
#' subtracting and random intercept factor with equal unstandardized loadings
#' from all the regular and unrecoded reversed items in the database
#'
#' @param data Matrix or data frame.
#' Should consist only of variables to be used in the analysis.
#' \strong{Must} be raw data and not a correlation matrix
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
#' \item \code{"auto"} --- Automatically computes appropriate correlations for
#' the data using Pearson's for continuous, polychoric for ordinal,
#' tetrachoric for binary, and polyserial/biserial for ordinal/binary with
#' continuous. To change the number of categories that are considered
#' ordinal, use \code{ordinal.categories}
#' (see \code{\link[EGAnet]{polychoric.matrix}} for more details)
#'
#' \item \code{"cor_auto"} --- Uses \code{\link[qgraph]{cor_auto}} to compute correlations.
#' Arguments can be passed along to the function
#'
#' \item \code{"pearson"} --- Pearson's correlation is computed for all
#' variables regardless of categories
#'
#' \item \code{"spearman"} --- Spearman's rank-order correlation is computed
#' for all variables regardless of categories
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
#' \item \code{"pairwise"} --- Computes correlation for all available cases between
#' two variables
#'
#' \item \code{"listwise"} --- Computes correlation for all complete cases in the dataset
#'
#' }
#'
#' @param model Character (length = 1).
#' Defaults to \code{"glasso"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"BGGM"} --- Computes the Bayesian Gaussian Graphical Model.
#' Set argument \code{ordinal.categories} to determine
#' levels allowed for a variable to be considered ordinal.
#' See \code{?BGGM::estimate} for more details
#'
#' \item \code{"glasso"} --- Computes the GLASSO with EBIC model selection.
#' See \code{\link[EGAnet]{EBICglasso.qgraph}} for more details
#'
#' \item \code{"TMFG"} --- Computes the TMFG method.
#' See \code{\link[EGAnet]{TMFG}} for more details
#'
#' }
#'
#' @param algorithm Character or
#' \code{igraph} \code{cluster_*} function (length = 1).
#' Defaults to \code{"walktrap"}.
#' Three options are listed below but all are available
#' (see \code{\link[EGAnet]{community.detection}} for other options):
#'
#' \itemize{
#'
#' \item \code{"leiden"} --- See \code{\link[igraph]{cluster_leiden}} for more details
#'
#' \item \code{"louvain"} --- By default, \code{"louvain"} will implement the Louvain algorithm using
#' the consensus clustering method (see \code{\link[EGAnet]{community.consensus}}
#' for more information). This function will implement
#' \code{consensus.method = "most_common"} and \code{consensus.iter = 1000}
#' unless specified otherwise
#'
#' \item \code{"walktrap"} --- See \code{\link[igraph]{cluster_walktrap}} for more details
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
#' \item \code{"expand"} --- Expands the correlation matrix with four variables correlated 0.50.
#' If number of dimension returns 2 or less in check, then the data
#' are unidimensional; otherwise, regular EGA with no matrix
#' expansion is used. This method was used in the Golino et al.'s (2020)
#' \emph{Psychological Methods} simulation
#'
#' \item \code{"LE"} --- Applies the Leading Eigenvector algorithm
#' (\code{\link[igraph]{cluster_leading_eigen}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Leading Eigenvector solution is used; otherwise, regular EGA
#' is used. This method was used in the Christensen et al.'s (2023)
#' \emph{Behavior Research Methods} simulation
#'
#' \item \code{"louvain"} --- Applies the Louvain algorithm (\code{\link[igraph]{cluster_louvain}})
#' on the empirical correlation matrix. If the number of dimensions is 1,
#' then the Louvain solution is used; otherwise, regular EGA is used.
#' This method was validated Christensen's (2022) \emph{PsyArXiv} simulation.
#' Consensus clustering can be used by specifying either
#' \code{"consensus.method"} or \code{"consensus.iter"}
#'
#' }
#'
#' @param plot.EGA Boolean (length = 1).
#' If \code{TRUE}, returns a plot of the network and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#' @param verbose Boolean (length = 1).
#' Whether messages and (insignificant) warnings should be output.
#' Defaults to \code{FALSE} (silent calls).
#' Set to \code{TRUE} to see all messages and warnings for every function call
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[EGAnet]{auto.correlate}},
#' \code{\link[EGAnet]{network.estimation}},
#' \code{\link[EGAnet]{community.detection}},
#' \code{\link[EGAnet]{community.consensus}}, and
#' \code{\link[EGAnet]{EGA}}
#'
#' @return Returns a list containing:
#'
#' \item{EGA}{Results from \code{\link[EGAnet]{EGA}}}
#'
#' \item{RI}{A list containing information about the random-intercept
#' model (if the model converged):
#'
#' \itemize{
#'
#' \item \code{fit} --- The fit object for the random-intercept model using \code{\link[lavaan]{cfa}}
#'
#' \item \code{lavaan.args} --- The arguments used in \code{\link[lavaan]{cfa}}
#'
#' \item \code{loadings} --- Standardized loadings from the random-intercept model
#'
#' \item \code{correlation} --- Residual correlations after accounting for the random-intercept model
#'
#' }
#'
#' }
#'
#' \item{TEFI}{\code{link[EGAnet]{tefi}} for the estimated structure}
#'
#' \item{plot.EGA}{Plot output if \code{plot.EGA = TRUE}}
#'
#' @author
#' Alejandro Garcia-Pardina <alejandrogp97@gmail.com>,
#' Francisco J. Abad <fjose.abad@uam.es>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s at virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu.do>, and
#' Robert Moulder <rgm4fd@virginia.edu>
#'
#' @examples
#' # Obtain example data
#' wmt <- wmt2[,7:24]
#'
#' # riEGA example
#' riEGA(data = wmt, plot.EGA = FALSE)
#' # no plot for CRAN checks
#'
#' @references
#' \strong{Selection of CFA Estimator} \cr
#' Rhemtulla, M., Brosseau-Liard, P. E., & Savalei, V. (2012).
#' When can categorical variables be treated as continuous? A comparison of robust continuous and categorical SEM estimation methods under suboptimal conditions.
#' \emph{Psychological Methods}, \emph{17}, 354-373.
#'
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{EGAnet}
#'
#' @export
#'
# Random-Intercept EGA
# Superceded 'residualEGA.R' on 17.04.2022
# Updated 12.06.2024
riEGA <- function(
    data, n = NULL,
    corr = c("auto", "cor_auto", "pearson", "spearman"),
    na.data = c("pairwise", "listwise"),
    model = c("glasso", "TMFG"),
    algorithm = c("leiden", "louvain", "walktrap"),
    uni.method = c("expand", "LE", "louvain"),
    plot.EGA = TRUE, verbose = FALSE,
    ...
)
{

  # Argument errors (return data in case of tibble)
  data <- riEGA_errors(data, n, plot.EGA, verbose, ...)

  # Check for missing arguments (argument, default, function)
  corr <- set_default(corr, "auto", riEGA)
  na.data <- set_default(na.data, "pairwise", auto.correlate)
  model <- set_default(model, "glasso", network.estimation)
  algorithm <- set_default(algorithm, "walktrap", community.detection)
  uni.method <- set_default(uni.method, "louvain", community.unidimensional)

  # Obtain correlations
  output <- obtain_sample_correlations(data, n, corr, na.data, verbose)

  # Catch BGGM
  if(model == "bggm"){
    stop(
      "'model = BGGM' is not supported in `riEGA` because it requires the original data and the residual correlation matrix from the random-intercept model is used as input into `EGA`.",
      call. = FALSE
    )
  }

  # Ensure data has names
  data <- ensure_dimension_names(data)

  # Check that data is not symmetric or square matrix
  if(is_symmetric(data)){
    stop(
      "A symmetric matrix was input into 'data'. The original data needs to be used to properly estimate the random-intercept model.",
      call. = FALSE
    )
  }

  # Get data dimensions
  dimensions <- dim(data)

  # Get variable names
  variable_names <- dimnames(data)[[2]]

  # Get ellipse arguments
  ellipse <- list(...)

  # Get {lavaan}'s CFA function
  cfa_FUN <- silent_load(lavaan::cfa)

  # Get {lavaan} CFA arguments
  lavaan_ARGS <- obtain_arguments(cfa_FUN, ellipse)

  # Set random-intercept model
  lavaan_ARGS$model <- paste(
    "RI =~", paste0("1*", variable_names, collapse = " + ")
  )

  # Set {lavaan} arguments
  lavaan_ARGS[c(
    "sample.cov", "sample.nobs", "estimator",
    "std.lv", "se"
  )] <- list(
    output$correlation_matrix, output$n,
    swiftelse("estimator" %in% names(ellipse), ellipse$estimator, "ML"),
    FALSE, "none"
  )

  # Fit model
  fit <- try(
    do.call(
      what = "cfa_FUN",
      # not sure why {lavaan}'s `cfa` function
      # needs to be in quotes...but it does
      args = lavaan_ARGS
    ),
    silent = TRUE
  )

  # Check for non-convergence
  if(is(fit, "try-error")){

    # Let user know that random-intercept model did not converge
    message(
      "Random-intercept model did not converge suggesting that there is not substantial evidence for wording effects."
    )

    # Get correlations since they won't come from random-intercept model
    correlation_matrix <- output$correlation_matrix

    # Get correlations since they won't come from random-intercept model
    if(corr == "auto"){

      # Compute correlations
      correlation_matrix <- auto.correlate(
        data = data, corr = "pearson",
        na.data = na.data, verbose = verbose,
        needs_usable = FALSE, # skips usable data check
        ...
      )

    }else{

      # Obtain correlations using base R
      correlation_matrix <- cor(data, use = na.data, method = corr)

    }

  }else{ # Proceed with random-intercept model

    # Obtain random-intercept loadings
    ri_loadings <- lavaan::inspect(fit, what = "std")$lambda

    # Obtain residual correlation matrix
    correlation_matrix <- lavaan::residuals(fit, type = "cor")$cov

    # Change diagonal to 1
    diag(correlation_matrix) <- 1

    # Ensure positive definite
    if(!is_positive_definite(correlation_matrix)){

      correlation_matrix <- as.matrix(
        Matrix::nearPD(
          correlation_matrix,
          corr = TRUE,
          keepDiag = TRUE
        )$mat
      )

    }else{ # Remove {lavaan} class
      correlation_matrix <- unclass(correlation_matrix)
    }

  }


  # Estimate EGA
  ega_result <- EGA(
    data = correlation_matrix, n = dimensions[1],
    corr = corr, na.data = na.data, model = model,
    algorithm = algorithm, uni.method = uni.method,
    plot.EGA = FALSE, verbose = verbose,
    needs_usable = FALSE, # skips usable data check
    ...
  )

  # Set up results
  results <- list(
    EGA = ega_result,
    RI = list(fit = fit, lavaan.args = lavaan_ARGS)
  )

  # Add methods for random-intercept model
  attr(results$RI, "methods") <- list(
    model = lavaan_ARGS$model,
    estimator = lavaan_ARGS$estimator,
    ordered = lavaan_ARGS$ordered
  )

  # Check for random-intercept model convergence
  if(!is(fit, "try-error")){

    # Send loadings and correlations
    results$RI$loadings <- ri_loadings
    results$RI$correlation <- correlation_matrix

    # Add loadings attribute
    attr(results$RI, "methods")$loading <- unname(unique(ri_loadings))

    # Send message about recoding (regardless of 'verbose')
    message(
      paste0(
        "The random-intercept model converged. Wording effects likely. Results are only valid if data are ",
        styletext("unrecoded", defaults = "underline"), "."
      )
    )

  }

  # Make class "riEGA"
  class(results) <- "riEGA"

  # Add TEFI to the result
  results$TEFI <- results$EGA$TEFI

  # Check for plot
  if(plot.EGA && sum(results$EGA$network != 0)){

    # Set up plot
    results$Plot.EGA <- plot(results, ...)

    # Actually send the plot
    silent_plot(results$Plot.EGA)

  }

  # Return results
  return(results)

}

# Bug checking ----
## Basic input
# data = wmt2[,7:24]; n = NULL; corr = "auto"; estimator = "auto"
# na.data = "pairwise"; model = "glasso"; algorithm = "leiden"
# uni.method = "louvain"; plot.EGA = TRUE; verbose = FALSE;
# ellipse = list()

#' @noRd
# Errors ----
# Updated 07.09.2023
riEGA_errors <- function(data, n, plot.EGA, verbose, ...)
{

  # 'data' errors
  object_error(data, c("matrix", "data.frame", "tibble"), "riEGA")

  # Check for tibble
  if(get_object_type(data) == "tibble"){
    data <- as.data.frame(data)
  }

  # 'n' errors
  if(!is.null(n)){
    length_error(n, 1, "riEGA")
    typeof_error(n, "numeric", "riEGA")
  }

  # 'plot.EGA' errors
  length_error(plot.EGA, 1, "riEGA")
  typeof_error(plot.EGA, "logical", "riEGA")

  # 'verbose' errors
  length_error(verbose, 1, "riEGA")
  typeof_error(verbose, "logical", "riEGA")

  # Check for usable data
  if(needs_usable(list(...))){
    data <- usable_data(data, verbose)
  }

  # Return usable data in case of tibble
  return(data)

}

#' @exportS3Method
# S3 Print Method ----
# Updated 28.06.2023
print.riEGA <- function(x, ...)
{

  # Print regular EGA result
  print(x$EGA)

  # Add random-intercept model
  ri_attributes <- attr(x$RI, "methods")

  # Add break space
  cat("\n\n----\n\n")

  # Print estimator
  cat("Random-Intercept Estimator: ", ri_attributes$estimator, "\n")
  cat("Random-Intercept Loading: ", format_decimal(ri_attributes$loading, 3))

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 05.07.2023
summary.riEGA <- function(object, ...)
{
  print(object, ...) # same as print
}

#' @exportS3Method
# S3 Plot Method ----
# Updated 28.06.2023
plot.riEGA <- function(x, ...)
{

  # Return plot
  single_plot(
    network = x$EGA$network,
    wc = x$EGA$wc,
    ...
  )

}

