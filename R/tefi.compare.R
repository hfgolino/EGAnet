#' Compare Total Entropy Fit Index (TEFI) Between Two Structures
#'
#' This function computes the TEFI values for two different structures using 
#' bootstrapped correlation matrices from `bootEGA` and compares them using a non-parametric bootstrap test.
#' It also visualizes the distributions of TEFI values for both structures.
#' 
#' The null hypothesis being tested is that the TEFI values obtained in the bootstrapped correlation matrices for the `base` structure
#' are `higher` than the TEFI values obtained in the bootstrapped correlation matrices for the `comparison` structure.
#' 
#' Therefore, the `p-value` in this bootstrap test can be interpreted as follows:
#' 
#' \describe{
#'   \item{If the p-value is small (e.g., < 0.01):}{There is evidence to suggest that the TEFI values for the `base` structure tend to be lower 
#'   than those for `comparison`, indicating that the former might provide a better fit (lower entropy) compared to the latter. 
#'   In other words, if the null hypothesis is true, it is very unlikely to obtain a test statistic as extreme as the one we observe in the data using a bootstrap strategy.}
#'   
#'   \item{If the p-value is large (e.g., > 0.01):}{There is insufficient evidence to conclude that the `base` structure has lower TEFI than `comparison`, 
#'   suggesting that both structures may provide similar fits or that `comparison` might fit better. 
#'   In other words, if the null hypothesis is true, it is not unlikely to obtain a test statistic as extreme as the one we observe in the data using a bootstrap strategy.}
#' }
#'
#' @param bootega.obj A `bootEGA` object containing the bootstrapped correlation matrices.
#' @param base A numeric vector representing the first structure to be tested.
#' @param comparison A numeric vector representing the second structure to be tested.
#' @param base.name A character string representing the label for the first structure in the plot. Default is `"Base"`.
#' @param comparison.name A character string representing the label for the second structure in the plot. Default is `"Comparison"`.
#' @param base.color A character string specifying the color for the first structure in the plot. Default is `"blue"`.
#' @param comparison.color A character string specifying the color for the second structure in the plot. Default is `"red"`.
#' @param plot.TEFI Logical. Indicates if the plot of the TEFI comparison and the p-value should be plotted. Default is `TRUE`.
#' @return A list containing:
#' \describe{
#'   \item{TEFI.base}{A numeric vector of TEFI values for the first structure.}
#'   \item{TEFI.comparison}{A numeric vector of TEFI values for the second structure.}
#'   \item{p_value}{The p-value from the non-parametric bootstrap hypothesis test.}
#'   \item{df}{A data frame containing the TEFI values for both structures.}
#' }
#'
#' @examples
#' simulated <- simEGM(
#'   communities = 2, variables = 6,
#'   loadings = 0.7, correlations = 0.3,
#'   sample.size = 1000
#' )
#' 
#' boot.sim <- bootEGA(
#'   data = simulated$data, iter = 500,
#'   type = "parametric", ncores = 20
#' )
#' 
#' comparingTEFI <- tefi.compare(
#'   boot.sim,
#'   base = boot.sim$EGA$wc,
#'   comparison = c(1,1,1,2,2,2,3,3,3,4,4,4)
#' )
#' 
#' plot(comparingTEFI, base.name = "EGA Structure", comparison.name = "Alternative Structure", 
#'      base.color = "blue", comparison.color = "red")
#'
#' @export

tefi.compare <- function(bootega.obj, base, comparison, plot.TEFI = TRUE, base.color = "blue", comparison.color = "red") {
  # Calculate TEFI for each correlation matrix and structure
  ## First structure
  TEFI.base <- EGAnet:::nvapply(
    bootega.obj = bootEGA$bootCorrs, function(R){
      return(tefi(data = R, structure = base)$VN.Entropy.Fit)
    }
  )
  ## Second structure
  TEFI.comparison <- EGAnet:::nvapply(
    bootega.obj$bootCorrs, function(R){
      return(tefi(data = R, structure = comparison)$VN.Entropy.Fit)
    }
  )
  
  # Calculate p-value
  p <- mean(c(TRUE, TEFI.base >= TEFI.comparison), na.rm = TRUE)
  
  # Create a data frame for visualization
  df <- data.frame(
    TEFI = c(TEFI.base, TEFI.comparison),
    Structure = rep(c("Base", "Comparison"), each = bootega.obj$iter)
  )
  
  result <- list(
    TEFI.base = TEFI.base,
    TEFI.comparison = TEFI.comparison,
    p_value = p,
    df = df
  )
  
  # Set class
  class(result) <- "tefi.compare"
  
  # Plot
  if(plot.TEFI){
    EGAnet:::silent_plot(
      plot(
        result,
        base.color = base.color, 
        comparison.color = comparison.color
      )
    )
  }
  
  return(result)
}

# S3 method for plotting tefi.compare
#' @noRd
plot.tefi.compare <- function(x, base.name = "Base", comparison.name = "Comparison", 
                              base.color = "blue", comparison.color = "red", ...) {
  if (!inherits(x, "tefi.compare")) {
    stop("Input must be a tefi.compare object")
  }
  
  ggplot(x$df, aes(x = TEFI, fill = Structure)) +
    geom_density(alpha = 0.5) +
    geom_vline(aes(xintercept = mean(x$df$TEFI[x$df$Structure == base.name])), color = base.color, linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = mean(x$df$TEFI[x$df$Structure == comparison.name])), color = comparison.color, linetype = "dashed", size = 1) +
    labs(title = "TEFI Comparison",
         subtitle = paste("p-value:", round(x$p_value, 4)),
         x = "TEFI Value",
         y = "Density",
         fill = "Structure") +
    theme_minimal() +
    scale_fill_manual(values = c(base.color, comparison.color), labels = c(base.name, comparison.name)) +
    theme(legend.title = element_blank())
}
