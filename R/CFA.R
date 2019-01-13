#' CFA Fit of EGA Structure
#'
#' Verifies the fit of the structure suggested by EGA using confirmatory factor analysis
#'
#' @param ega.obj An \code{\link[EGA]{EGA}} object
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
#' #estimate EGA
#' ega.wmt <- EGA(data = wmt2[,7:24])
#' 
#' #fit CFA model to EGA results
#' cfa.wmt <- CFA(ega.obj = ega.wmt, estimator = 'WLSMV', plot.CFA = TRUE, data = wmt2)
#'
#' #additional fit measures
#' lavaan::fitMeasures(cfa.wmt$fit, fit.measures = "all")
#'
#' #estimate EGA
#' ega.intel <- EGA(data = intelligenceBattery[,8:66])
#' 
#' #fit CFA model to EGA results
#' cfa.intel <- CFA(ega.obj = ega.intel, estimator = 'WLSMV', plot.CFA = TRUE,
#' data = intelligenceBattery[,8:66])
#' 
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{bootEGA}} to investigate the stability of EGA's estimation via bootstrap.
#'
#' @export
#CFA model for EGA
CFA <- function(ega.obj, data, estimator, plot.CFA = TRUE, layout = "spring", ...) {
    
    strct <- split(ega.obj$dim.variables[, 1], list(ega.obj$dim.variables[, 2]))
    names(strct) <- paste("Fat", labels(strct))
    model.ega <- paste(names(strct), " =~ ", lapply(strct, function(x) paste(print(x), collapse = " + ")), collapse = " \n ")
    fit.mod.ega <- lavaan::cfa(model = model.ega, estimator = estimator, orthogonal = FALSE, se = "standard", test = "satorra-bentler",
                               data = data, ...)
    summary.cfa <- summary(fit.mod.ega, fit.measures = TRUE)
    fit.measures.cfa <- lavaan::fitMeasures(fit.mod.ega, fit.measures = c("chisq", "df", "pvalue", "cfi", "rmsea", "gfi", "nfi"))
    
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
