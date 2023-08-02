#' CFA Fit of \code{\link[EGAnet]{EGA}} or \code{\link[EGAnet]{hierEGA}} Structure
#'
#' Verifies the fit of the structure suggested by \code{\link[EGAnet]{EGA}} or by \code{\link[EGAnet]{hierEGA}} using confirmatory factor analysis
#'
#' @param ega.obj An \code{\link[EGAnet]{EGA}} object or an \code{\link[EGAnet]{hierEGA}}
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
#' \dontrun{
#' # Estimate EGA
#' ega.wmt <- EGA(
#'   data = wmt,
#'   plot.EGA = FALSE # No plot for CRAN checks
#' )
#' 
#' # Fit CFA model to EGA results
#' cfa.wmt <- CFA(
#'   ega.obj = ega.wmt, estimator = "WLSMV",
#'   plot.CFA = FALSE, # No plot for CRAN checks
#'   data = wmt
#' )
#'
#' # Additional fit measures
#' lavaan::fitMeasures(cfa.wmt$fit, fit.measures = "all")}
#' 
#' @references 
#' Christensen, A. P., Gross, G. M., Golino, H., Silvia, P. J., & Kwapil, T. R. (2019).
#' Exploratory graph analysis of the Multidimensional Schizotypy Scale.
#' \emph{Schizophrenia Research}, \emph{206}, 43-51.
#' 
#' Golino, H., & Epskamp, S. (2017).
#' Exploratory graph analysis: A new approach for estimating the number of dimensions in psychological research.
#' \emph{PLoS ONE}, \emph{12}, e0174035.
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{bootEGA}} to investigate the stability of EGA's estimation via bootstrap.
#'
#' @export
#'
# CFA model for EGA
# Updated 29.11.2022
CFA <- function(ega.obj, data, estimator, plot.CFA = TRUE, layout = "spring", ...) 
{
  
  # Set lavaan arguments
  lavaan.args <- list(...)
  
  ## lavaan.args
  if(length(lavaan.args) == 0){
    lavaan.args <- formals(lavaan::cfa)
    lavaan.args[length(lavaan.args)] <- NULL
    lavaan.args$std.lv <- TRUE
  }else{
    lavaan.default <- formals(lavaan::cfa)
    lavaan.default[length(lavaan.default)] <- NULL
    lavaan.default$std.lv <- TRUE
    
    if(any(names(lavaan.args) %in% names(lavaan.default))){
      lavaan.default[names(lavaan.args)] <- lavaan.args
    }
    
    lavaan.args <- lavaan.default
  }
  
  ## change key if NULL
  data <- lavaan.formula.names(data)
  
  ## Get default estimator
  categories <- apply(data, 2, function(x){
    length(na.omit(unique(x)))
  })
  
  ## Check categories
  if(any(categories < 6)){# Not all continuous
    lavaan.args$estimator <- "WLSMV"
    lavaan.args$missing <- "pairwise"
    lavaan.args$ordered <- TRUE
  }else{# All can be considered continuous
    lavaan.args$estimator <- "MLR"
    lavaan.args$missing <- "fiml"
  }
  
  
  if(class(ega.obj)=="EGA"){
    strct <- split(ega.obj$dim.variables[, 1], list(ega.obj$dim.variables[, 2]))
    names(strct) <- paste("Fat", labels(strct))
    model.ega <- paste(names(strct), " =~ ", lapply(strct, function(x) paste(print(x), collapse = " + ")), collapse = " \n ")
    lavaan.args$model <- model.ega
    lavaan.args$data <- data
    FUN <- lavaan::cfa
    fit.mod.ega <- do.call(
      what = "FUN",
      args = lavaan.args
    )
}
  if(class(ega.obj)=="hierEGA"){
     ## Organizing 1st level items and factors:
    hierEGA.obj <- ega.obj
    
    strct <- split(hierEGA.obj$lower_order$dim.variables[, 1], list(hierEGA.obj$lower_order$dim.variables[, 2]))
    names(strct) <- paste("FirstFact", labels(strct))
    
    ## Organizing 2nd level items and factors:
    
    dim.variables.second <- data.frame(items = colnames(hierEGA.obj$lower_order$network[]), 
                                       dimension = hierEGA.obj$higher_order$wc[as.character(hierEGA.obj$lower_order$wc)])
    
    strct.2nd <- split(dim.variables.second[, 1], list(paste("GenFact", dim.variables.second[, 2])))
    
    if(length(names(strct.2nd))>1){
      ## Pairing the names of the 2nd level factors:
      pair.names <- t(combn(names(strct.2nd),2))
      
      ## Lavaan Model based on hierEGA results:
      model.ega <- paste(paste(names(strct), " =~ ", lapply(strct, function(x) paste(print(x), collapse = " + ")), collapse = " \n "), 
                         paste(names(strct.2nd), " =~ ", lapply(strct.2nd, function(x) paste(print(x), collapse = " + ")), collapse = " \n "), 
                         paste(apply(pair.names, 1, function(x) paste(x[1],"~~", x[2], sep = " ")), collapse = " \n "), sep = " \n ")  
    }else{
      ## Lavaan Model based on hierEGA results:
      model.ega <- paste(paste(names(strct), " =~ ", lapply(strct, function(x) paste(print(x), collapse = " + ")), collapse = " \n "), 
                         paste(names(strct.2nd), " =~ ", lapply(strct.2nd, function(x) paste(print(x), collapse = " + ")), collapse = " \n "), sep = " \n ") 
    }
    
    lavaan.args$model <- model.ega
    lavaan.args$data <- data
    lavaan.args$orthogonal <- TRUE

    FUN <- lavaan::cfa
    fit.mod.ega <- do.call(
      what = "FUN",
      args = lavaan.args
    )
    }
  
  summary.cfa <- summary(fit.mod.ega, fit.measures = TRUE)
  
  if (estimator == "WLSMV") {
    fit.measures.cfa <- lavaan::fitMeasures(fit.mod.ega, fit.measures = c("chisq.scaled", "df", "pvalue", "cfi.scaled", "rmsea.scaled"))
  } else {
    fit.measures.cfa <- lavaan::fitMeasures(fit.mod.ega, fit.measures = c("chisq", "df", "pvalue", "cfi", "rmsea", "gfi", "nfi"))
  }
  
  if (plot.CFA == TRUE) {
    plot.cfa <- semPlot::semPaths(fit.mod.ega, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
                                  sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = layout,
                                  "std", cut = 0.5)
  }
  
    cfa <- list()
    cfa$model <- model.ega
    cfa$fit <- fit.mod.ega
    cfa$summary <- summary.cfa
    cfa$fit.measures <- fit.measures.cfa
    class(cfa) <- "CFA"
    return(cfa)
}
