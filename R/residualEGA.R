#' Residualized \code{\link[EGAnet]{EGA}}
#'
#' \code{residualEGA} Estimates the number of dimensions after controlling for wording effects.
#' EGA is applied in the residual of a random intercept item factor model (RIIFA) with one method factor and one substantive factor.
#'
#' @param data Matrix or data frame.
#' Includes the variables to be used in the \code{residualEGA} analysis
#'
#' @param manifests Character vector.
#' Vector indicating the names of the variables (items) to be used in the analysis.
#'
#' @param lat Numeric integer.
#' Number of latent factors to be estimated.
#' Only one substantitve latent factor is recommended in the current version of the function.
#'
#' @param negative.items Numeric vector
#' A numeric vector indicating the column of the negative items.
#'
#' @param plot Boolean.
#' If \code{TRUE}, returns a plot of the residualized network and its estimated dimensions.
#' Defaults to \code{TRUE}
#'
#'
#' @return Returns a list containing:
#'
#' \item{openMx.model}{OpenMX model}
#'
#' \item{openMx.result}{OpenMX results}
#'
#' \item{openMx.std.par}{OpenMX standardized parameters}
#'
#' \item{ResidualMatrix}{Residual matrix}
#'
#' \item{EGA.Residuals}{Results of the residualized EGA}
#'
#' \item{Fit}{Fit metrics of the network structure, calculated using the ggmfit function of the
#' \code{\link[qgraph]{qgraph}} package}
#'
#' \item{WordLoads}{Loadings of the wording effects}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Robert Moulder <rgm4fd@virginia.edu>
#'
#' @examples
#'
#'
#'
#' \dontrun{
#'
#' # resEGA example
#' opt.res <- residualEGA(data = optmism, manifests = colnames(optmism),lat = 1, negative.items = c(3,7,9), plot = TRUE)
#'
#' # Fit:
#' opt.res$Fit
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA
#' and \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @importFrom OpenMx mxTryHard mxModel mxPath mxData mxRefModels mxStandardizeRAMpaths
#' @importFrom stats cov median sd qt
#' @importFrom utils combn
#'
#' @export
#'


### Residualized EGA ###

residualEGA <- function(data, manifests,lat, negative.items, plot = TRUE){
  factors <- paste0("F", seq(lat))
  wording <- paste0("W", 1)
  latents.2 <- c(factors, wording)

  #Substantive Latent Factors:

  substantive <- list()
  for (i in 1:length(factors)) {
    substantive <- append(substantive, OpenMx::mxPath(
      from = factors[i], to = manifests,
      arrows = 1,
      free   = c(FALSE, rep(TRUE, length(manifests)-1)),
      values = 1,
      labels=paste("L_",manifests,"_", factors[i],sep="")))
  }

  #Wording Factor:
  initial.values <- rep(1, length(manifests))
  initial.values[negative.items] <- -1
  length.word.factor <- length(manifests)/length(wording)

  wording.factors <- list()
  for (i in 1:length(wording)) {
    j <- i*length.word.factor
    wording.factors <- append(wording.factors, OpenMx::mxPath(
      from = wording[i], to = manifests[(j - (length.word.factor-1)):j],
      arrows = 1,
      free   = rep(FALSE, length.word.factor),
      values = initial.values[1:length.word.factor],
      labels=paste("W_",manifests[(j - (length.word.factor-1)):j],"_", wording[i],sep="")))
  }

  # Variance Factors:
  variances.factors <- OpenMx::mxPath(from = factors,
                              to=factors, arrows = 2, free = TRUE, values = .2,
                              labels = "Var_Factors")

  # Covariance Factors:
  if(lat>= 2){
    uniqueFactors <- utils::combn(factors,2)
    covariances.factors <- OpenMx::mxPath(from = uniqueFactors[1,],
                                  to=uniqueFactors[2,], arrows = 2, free = TRUE,
                                  labels = "Cov_Factors")
  } else{
    covariances.factors <- OpenMx::mxPath(from = factors[1],
                                  to=factors[1], arrows = 2, free = TRUE,
                                  labels = "Cov_Factors")
  }



  # Variance Wording:
  variances.wording <- OpenMx::mxPath(from = wording,
                              to=wording, arrows = 2, free = TRUE, values = 0.8,
                              labels="Var_Wording")

  # Variance Manifests:
  variances.manifests <- OpenMx::mxPath(from = manifests,
                                to=manifests, arrows = 2, free = TRUE, values = 0.8)


  # SEM Model:
  mod.ex <- OpenMx::mxModel('SEM',type='RAM',manifestVars=manifests,latentVars=latents.2,
                    substantive,
                    variances.factors,
                    covariances.factors,
                    wording.factors,
                    variances.wording,
                    variances.manifests,
                    OpenMx::mxPath(from="one",to=manifests,arrows=1,
                           values=.2,free=TRUE,
                           labels=paste("m_",manifests,sep="")),
                    OpenMx::mxData(observed=data,type='raw'))

  # SEM:
  modResult<-suppressWarnings(OpenMx::mxTryHard(mod.ex, 500))

  # Summary with Fit Indexes
  sum1 <- suppressWarnings(summary(modResult, refModels=OpenMx::mxRefModels(modResult, run = TRUE)))

  # Fit indexes:
  sum1$CFI
  sum1$RMSEA
  sum1$RMSEACI
  sum1$TLI
  sum1$Chi
  sum1$ChiDoF
  sum1$p

  #Standardized parameters:
  std.results <- OpenMx::mxStandardizeRAMpaths(modResult,SE=FALSE)

  # Calculating the residual matrix:


  # std. loadings methods (wording) factor
  n.FatLoad <- length(manifests)*lat

  word.loadings <- std.results[,8][(n.FatLoad+1):(n.FatLoad+length(manifests))]

  # Sample Cor:
  msamplecorr <- cor(data, use = "pairwise.complete.obs")

  # Residual Matrix:
  mresidual <- msamplecorr - (word.loadings %*% t(word.loadings))

  # EGA of the Residuals:
  res.ega <- EGA(mresidual, n = nrow(data), plot.EGA = FALSE, model = "glasso")
  cov.data <- cov(data, use="pairwise.complete.obs")

  fit.res <- qgraph::ggmFit(res.ega$network, cov.data, sampleSize = nrow(data))

  results <- list()
  results$openMx.model <- mod.ex
  results$openMx.result <- sum1
  results$openMx.std.par <- std.results
  results$ResidualMatrix <- mresidual
  results$EGA.Residuals <- res.ega
  results$Fit <- fit.res
  results$WordLoads <- word.loadings

  if (plot == TRUE) {
    plot.ega <- qgraph::qgraph(res.ega$network, layout = "spring",
                       vsize = 6, groups = as.factor(res.ega$wc))
  }
  class(results) <- "resEGA"
  return(results)
}
