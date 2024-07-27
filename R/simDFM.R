#' Simulate data following a Dynamic Factor Model
#'
#' @description Function to simulate data following a dynamic factor model (DFM). Two DFMs are currently available:
#' the direct autoregressive factor score model (Engle & Watson, 1981; Nesselroade, McArdle, Aggen, and Meyers, 2002) and the
#' dynamic factor model with random walk factor scores.
#'
#' @param variab Number of variables per factor.
#'
#' @param timep Number of time points.
#'
#' @param nfact Number of factors.
#'
#' @param error Value to be used to construct a diagonal matrix Q. This matrix is p x p covariance matrix Q that will
#' generate random errors following a multivariate normal distribution with mean zeros.
#' The value provided is squared before constructing Q.
#'
#' @param dfm A string indicating the dynamical factor model to use.
#' Current options are:
#'
#' \itemize{
#'
#' \item \strong{\code{DAFS}} --- Simulates data using the direct autoregressive factor score model.
#' This is the default method
#'
#' \item \strong{\code{RandomWalk}} --- Simulates data using a dynamic factor model with random walk factor scores
#'
#'}
#'
#' @param loadings Magnitude of the loadings.
#'
#' @param autoreg Magnitude of the autoregression coefficients.
#'
#' @param crossreg Magnitude of the cross-regression coefficients.
#'
#' @param var.shock Magnitude of the random shock variance.
#'
#' @param cov.shock Magnitude of the random shock covariance
#'
#' @param burnin Number of n first samples to discard when computing the factor scores. Defaults to 1000.
#'
#' @examples
#' \dontrun{
#' # Estimate EGA network
#' data1 <- simDFM(variab = 5, timep = 50, nfact = 3, error = 0.05,
#' dfm = "DAFS", loadings = 0.7, autoreg = 0.8,
#' crossreg = 0.1, var.shock = 0.36,
#' cov.shock = 0.18, burnin = 1000)}
#'
#' @references
#' Engle, R., & Watson, M. (1981).
#' A one-factor multivariate time series model of metropolitan wage rates.
#' \emph{Journal of the American Statistical Association}, \emph{76}(376), 774-781.
#'
#' Nesselroade, J. R., McArdle, J. J., Aggen, S. H., & Meyers, J. M. (2002).
#' Dynamic factor analysis models for representing process in multivariate time-series. In D. S. Moskowitz & S. L. Hershberger (Eds.),
#' \emph{Multivariate applications book series. Modeling intraindividual variability with repeated measures data: Methods and applications}, 235-265.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @export
# Simulate dynamic factor model ----
# Updated 27.07.2024
simDFM <- function(
    variab, timep, nfact, error, dfm = c("DAFS","RandomWalk"),
    loadings, autoreg, crossreg, var.shock, cov.shock, burnin = 1000
)
{

  # Check for missing arguments (argument, default, function)
  dfm <- set_default(dfm, "DAFS", simDFM)

  # Check for errors
  simDFM_errors(
    variab, timep, nfact, error,
    loadings, autoreg, crossreg,
    var.shock, cov.shock, burnin
  )

  # Get N
  N <- burnin + timep

  # Compute scores
  if(dfm == "dafs"){ # Dynamic Factor Model

    # B = Matrix of Bl is a nfact x nfact matrix containing the
    # autoregressive and cross-regressive coefficients
    B <- matrix(crossreg, ncol = nfact, nrow = nfact)
    diag(B) <- autoreg

    # Shock = Random shock vectors following a multivariate normal
    # distribution with mean zeros and nfact x nfact q covariance matrix D
    D <- matrix(cov.shock, nfact, nfact)
    diag(D) <- var.shock

    # Set up shock
    Shock <- MASS_mvrnorm_quick(p = nfact, np = N * nfact, coV = D)

    # Get factor scores
    Fscores <- matrix(0, nrow = N, ncol = nfact)
    Fscores[1,] <- Shock[1,]

    # Loop over to get scores
    for(time in 2:N){
      Fscores[time,] <- Fscores[time-1,] %*% B + Shock[time,]
    }

    # Remove burn-in
    Fscores <- Fscores[-seq_len(burnin),]

  }else{ # Random Walk

    # Initialize factor scores
    Fscores <- matrix(
      rnorm_ziggurat(n = nfact * N),
      nrow = nfact, ncol = N
    )

    # Remove burn-in
    Fscores <- Fscores[,-seq_len(burnin)]

    # nfact x timep matrix of scaled latent trends
    Fscores <- scale(apply(Fscores, 1, cumsum))

  }

  # Get loading matrix
  LoadMat <- as.matrix(Matrix::bdiag(lapply(rep(loadings, nfact), rep, variab)))

  # Set number of variables
  p <- variab * nfact

  # Set up for generating data
  Q <- diag(error^2, nrow = p, ncol = p)
  e <- t(MASS_mvrnorm_quick(p = p, np = timep * p, coV = Q))

  # Simulate observed data
  observed <- t(tcrossprod(LoadMat, Fscores) + e)

  # Get factor sequence
  factor_sequence <- seq_len(nfact)

  # Set names
  dimnames(observed)[[2]] <- paste0("V", seq_len(p))
  dimnames(Fscores)[[2]] <- paste0("F", factor_sequence)

  # Return results
  return(
    list(
      data = observed,
      Fscores = Fscores,
      LoadMat = LoadMat,
      Structure = rep(factor_sequence, each = variab)
    )
  )

}

#' @noRd
# Errors ----
# Updated 27.07.2024
simDFM_errors <- function(
    variab, timep, nfact, error,
    loadings, autoreg, crossreg,
    var.shock, cov.shock, burnin
)
{

  # 'variab' errors
  length_error(variab, 1, "simDFM")
  typeof_error(variab, "numeric", "simDFM")
  range_error(variab, c(2, Inf), "simDFM")

  # 'timep' errors
  length_error(timep, 1, "simDFM")
  typeof_error(timep, "numeric", "simDFM")
  range_error(timep, c(2, Inf), "simDFM")

  # 'nfact' errors
  length_error(nfact, 1, "simDFM")
  typeof_error(nfact, "numeric", "simDFM")
  range_error(nfact, c(1, Inf), "simDFM")

  # 'error' errors
  length_error(error, 1, "simDFM")
  typeof_error(error, "numeric", "simDFM")
  range_error(error, c(0, 1), "simDFM")

  # 'loadings' errors
  length_error(loadings, 1, "simDFM")
  typeof_error(loadings, "numeric", "simDFM")
  range_error(loadings, c(-1, 1), "simDFM")

  # 'autoreg' errors
  length_error(autoreg, 1, "simDFM")
  typeof_error(autoreg, "numeric", "simDFM")
  range_error(autoreg, c(0, 1), "simDFM")

  # 'crossreg' errors
  length_error(crossreg, 1, "simDFM")
  typeof_error(crossreg, "numeric", "simDFM")
  range_error(crossreg, c(0, 1), "simDFM")

  # 'var.shock' errors
  length_error(var.shock, 1, "simDFM")
  typeof_error(var.shock, "numeric", "simDFM")
  range_error(var.shock, c(0, 1), "simDFM")

  # 'cov.shock' errors
  length_error(cov.shock, 1, "simDFM")
  typeof_error(cov.shock, "numeric", "simDFM")
  range_error(cov.shock, c(0, 1), "simDFM")

  # 'burnin' errors
  length_error(burnin, 1, "simDFM")
  typeof_error(burnin, "numeric", "simDFM")
  range_error(burnin, c(1, Inf), "simDFM")

}




