#' CFA Fit of \code{\link[EGAnet]{EGA}} Structure
#'
#' Verifies the fit of the structure suggested by \code{\link[EGAnet]{EGA}} using confirmatory factor analysis
#'
#' @param ega.obj An \code{\link[EGAnet]{EGA}} object
#'
#' @param data A dataframe with the variables to be used in the analysis
#'
#' @param estimator The estimator used in the confirmatory factor analysis.
#' 'WLSMV' is the estimator of choice for ordinal variables.
#' 'ML' or 'WLS' for interval variables.
#' See \code{\link[lavaan]{lavOptions}} for more details
#'
#' @param plot.CFA Logical.
#' Should the CFA structure with its standardized loadings be plot?
#' Defaults to TRUE
#'
#' @param layout Layout of plot (see \code{\link[semPlot]{semPaths}}).
#' Defaults to "spring"
#'
#' @param ... Arguments passed to \code{\link[lavaan]{cfa}}
#'
#' @return Returns a list containing:
#'
#' \item{fit}{Output from \code{\link[lavaan]{cfa}}}
#'
#' \item{summary}{Summary output from \code{\link[lavaan]{lavaan-class}}}
#'
#' \item{fit.measures}{Fit measures: chi-squared,
#' degrees of freedom, p-value, CFI, RMSEA, GFI, and NFI.
#' Additional fit measures can be applied using the
#' \code{\link[lavaan]{fitMeasures}} function (see examples)}
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#'
#' @examples
#' # Load data
#' wmt <- wmt2[,7:24]
#' 
#' \dontshow{# Fast for CRAN
#' cor.wmt <- cor(wmt)
#' 
#' # Estimate EGA
#' ega.wmt <- EGA(data = wmt, uni = TRUE, n = nrow(wmt2), plot.EGA = FALSE)
#' }
#'
#' \donttest{
#' # Estimate EGA
#' ega.wmt <- EGA(data = wmt, uni = TRUE, plot.EGA = FALSE)
#' 
#' # Fit CFA model to EGA results
#' cfa.wmt <- CFA(ega.obj = ega.wmt, estimator = 'WLSMV', plot.CFA = TRUE, data = wmt)
#'
#' # Additional fit measures
#' lavaan::fitMeasures(cfa.wmt$fit, fit.measures = "all")
#' }
#'
#' # Load data
#' intel <- intelligenceBattery[,8:66]
#'
#' \donttest{
#' # Estimate EGA
#' ega.intel <- EGA(data = intel, plot.EGA = FALSE)
#'
#' # Fit CFA model to EGA results
#' cfa.intel <- CFA(ega.obj = ega.intel, estimator = 'WLSMV', plot.CFA = TRUE,
#' data = intel)
#' }
#' 
#' @references 
#' Christensen, A. P., Gross, G. M., Golino, H., Silvia, P. J., & Kwapil, T. R. (2019).
#' Exploratory graph analysis of the Multidimensional Schizotypy Scale.
#' \emph{Schizophrenia Research}, \emph{206}, 43-51.
#' doi: \href{https://doi.org/10.1016/j.schres.2018.12.018}{10.1016/j.schres.2018.12.018}
#' 
#' Golino, H., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PloS one}, \emph{12(6)}, e0174035.
#' doi: \href{https://doi.org/10.1371/journal.pone.0174035}{10.1371/journal.pone.0174035}
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{bootEGA}} to investigate the stability of EGA's estimation via bootstrap.
#'
#' @export
#'
# CFA model for EGA
# Updated 21.10.2020
CFA<- function(ega.obj, data, estimator, plot.CFA = TRUE, layout = "spring", ...) {

    strct <- split(ega.obj$dim.variables[, 1], list(ega.obj$dim.variables[, 2]))
    names(strct) <- paste("Fat", labels(strct))
    model.ega <- paste(names(strct), " =~ ", lapply(strct, function(x) paste(print(x), collapse = " + ")), collapse = " \n ")
    fit.mod.ega <- lavaan::cfa(model = model.ega, estimator = estimator, orthogonal = FALSE,
                               data = data, ...)
    summary.cfa <- summary(fit.mod.ega, fit.measures = TRUE)

    if (estimator == "WLSMV") {
        fit.measures.cfa <- lavaan::fitMeasures(fit.mod.ega, fit.measures = c("chisq.scaled", "df", "pvalue", "cfi.scaled", "rmsea.scaled"))
    } else {
        fit.measures.cfa <- lavaan::fitMeasures(fit.mod.ega, fit.measures = c("chisq", "df", "pvalue", "cfi", "rmsea", "gfi", "nfi"))
    }

    if (plot.CFA == TRUE) {
        plot.cfa <- semPlot::semPaths(fit.mod.ega, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
                                      sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
                                      "std", cut = 0.5)
    }

    cfa <- list()
    cfa$fit <- fit.mod.ega
    cfa$summary <- summary.cfa
    cfa$fit.measures <- fit.measures.cfa
    class(cfa) <- "CFA"
    return(cfa)
}
