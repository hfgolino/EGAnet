#Methods:

#Plot EGA:
plot.EGA <- function(object, layout = "spring", vsize = 6, ...) {
  groups = as.factor(object$wc)
  variable = object$glasso
  qgraph(variable, layout = layout, vsize = vsize, groups = groups, ...)
}


#Summary EGA:
summary.EGA <- function(object, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(object$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dim.variables)
}

#Print EGA:
print.EGA <- function(object, ...) {
  cat("EGA Results:\n")
  cat("\nNumber of Dimensions:\n")
  print(object$n.dim)
  cat("\nItems per Dimension:\n")
  print(object$dim.variables)
}

#Plot CFA:
plot.CFA <- function(object, layout = "spring", vsize = 6, ...) {
  semPaths(object$fit, title = FALSE, label.cex = 0.8, sizeLat = 8, sizeMan = 5, edge.label.cex = 0.6, minimum = 0.1,
           sizeInt = 0.8, mar = c(1, 1, 1, 1), residuals = FALSE, intercepts = FALSE, thresholds = FALSE, layout = "spring",
           "std", cut = 0.5)
}

#Summary CFA:
summary.CFA <- function(object, ...) {
  cat("Summary: Confirmatory Factor Analysis:\n")
  print(object$summary)
  cat("\n FIt Measures:\n")
  print(object$fit.measures)
}

