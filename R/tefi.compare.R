#' @title Compare Total Entropy Fit Index (\code{\link[EGAnet]{tefi}}) Between Two Structures
#'
#' @description This function computes the \code{\link[EGAnet]{tefi}} values for two different structures using
#' bootstrapped correlation matrices from \code{\link[EGAnet]{bootEGA}} and compares them using a
#' non-parametric bootstrap test. It also visualizes the distributions of \code{\link[EGAnet]{tefi}} values
#' for both structures.
#'
#' @param bootega.obj A \code{\link[EGAnet]{bootEGA}} object
#'
#' @param base Numeric (length = columns in original dataset).
#' A vector representing the base structure to be tested
#'
#' @param comparison Numeric (length = columns in original dataset).
#' A vector representing the structure to be compared against the \code{base} structure
#'
#' @param plot.TEFI Boolean (length = 1).
#' Whether the TEFI comparison and the p-value should be plotted.
#' Defaults to \code{TRUE}
#'
#' @param ... Additional arguments that can be passed on to \code{\link[EGAnet]{plot.EGAnet}}.
#' See \code{Examples} for plotting arguments
#'
#' @return A list containing:
#'
#' \item{\code{TEFI.df}}{A data frame containing the TEFI values for both structures}
#'
#' \item{\code{p.value}}{The \emph{p}-value from the non-parametric bootstrap hypothesis test}
#'
#' @details The null hypothesis is that the TEFI values obtained in the bootstrapped correlation matrices for the \code{base}
#' structure are than the TEFI values obtained in the bootstrapped correlation matrices for the \code{comparison} structure.
#' Therefore, the \emph{p}-value in this bootstrap test can be interpreted as follows:
#'
#' \itemize{
#'
#' \item If the \emph{p}-value less than 0.05: TEFI values for the \code{base} structure tend to be lower
#' than the \code{comparison} structure, indicating that the former provides a better fit (lower entropy) than the latter
#'
#' \item If the \emph{p}-value is greater than 0.05: TEFI values for the \code{base} structure are not significantly lower than
#' the \code{comparison} structure, suggesting that both structures may provide similar fits or that \code{comparison} might fit better
#'
#' }
#'
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @examples
#' # Obtain data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Perform bootstrap EGA
#' boot.wmt <- bootEGA(
#'   data = simulated$data, iter = 500,
#'   type = "parametric", ncores = 2
#' )}
#'
#' # Perform comparison
#' comparing_tefi <- tefi.compare(
#'   boot.wmt,
#'   base = boot.wmt$EGA$wc, # Compare Walktrap
#'   comparison = community.detection(
#'    boot.wmt$EGA$network, algorithm = "louvain"
#'   ) # With Louvain
#' )
#'
#' # Plot options
#' plot(
#'   comparing_tefi,
#'   base.name = "Walktrap", base.color = "royalblue",
#'   comparison.name = "Louvain", comparison.color = "orange"
#' )
#'
#' @export
#'
# TEFI comparison ----
# Updated 28.01.2025
tefi.compare <- function(bootega.obj, base, comparison, plot.TEFI = TRUE, ...)
{

  # Argument errors
  tefi.compare_errors(bootega.obj, base, comparison, plot.TEFI, ...)

  # Calculate TEFI for the base structure
  TEFI.base <- nvapply(
    bootega.obj$bootCorrs, function(R){
      return(tefi(data = R, structure = base)$VN.Entropy.Fit)
    }
  )

  # Calculate TEFI for the comparison structure
  TEFI.comparison <- nvapply(
    bootega.obj$bootCorrs, function(R){
      return(tefi(data = R, structure = comparison)$VN.Entropy.Fit)
    }
  )

  # Create results
  result <- list(
    # Data frame
    TEFI.df = data.frame(
      TEFI = c(TEFI.base, TEFI.comparison),
      Structure = rep(c("Base", "Comparison"), each = bootega.obj$iter)
    ),
    # Calculate p-value
    p.value = mean(c(TRUE, TEFI.base >= TEFI.comparison), na.rm = TRUE)
  )

  # Set class
  class(result) <- "tefi.compare"

  # Plot
  if(plot.TEFI){
    silent_plot(plot(result, ...))
  }

  # Return result
  return(result)

}

#' @noRd
# Errors ----
# Updated 28.01.2025
tefi.compare_errors <- function(bootega.obj, base, comparison, plot.TEFI, ...)
{

  # 'bootega.obj' error
  if(!is(bootega.obj, "bootEGA")){

    # Throw error
    .handleSimpleError(
      h = stop,
      msg = paste0(
        "Object '", deparse(substitute(bootega.obj)), "' is not of ",
        "class 'bootEGA' (it is ", paste0("'", class(bootega.obj), "'", collapse = ", "),
        ")\n\'bootega.obj' input must be of class 'bootEGA'"
      ),
      call = "tefi.compare"
    )

  }

  # 'base' errors
  length_error(base, length(bootega.obj$EGA$wc), "tefi.compare")
  typeof_error(base, "numeric", "tefi.compare")

  # 'comparison' errors (at this point, `base` is OK)
  length_error(comparison, length(base), "tefi.compare")
  typeof_error(comparison, "numeric", "tefi.compare")

  # 'plot.TEFI' errors
  length_error(plot.TEFI, 1, "tefi.compare")
  typeof_error(plot.TEFI, "logical", "tefi.compare")

}

#' @exportS3Method
# S3 Print Method ----
# Updated 28.01.2025
print.tefi.compare <- function(x, ...)
{

  # Get values
  base_values <- x$TEFI.df$TEFI[x$TEFI.df$Structure == "Base"]
  comparison_values <- x$TEFI.df$TEFI[x$TEFI.df$Structure == "Comparison"]

  # Print summary
  cat(
    paste0(
      "Base TEFI: ", round(mean(base_values), 4),
      " (SD = ", round(sd(base_values), 4), ")",
      "\nComparison TEFI: ", round(mean(comparison_values), 4),
      " (SD = ", round(sd(comparison_values), 4), ")",
      "\n", styletext("p", defaults = "italics"),
      " (one-tailed, Base >= Comparison): ", round(x$p.value, 4)
    )
  )

}

#' @exportS3Method
# S3 Summary Method ----
# Updated 28.01.2025
summary.tefi.compare <- function(object, ...)
{

  # Same as print
  print(object, ...)

}

#' @exportS3Method
# S3 Plot Method ----
# Updated 28.01.2025
plot.tefi.compare <- function(
    x, base.name = "Base", base.color = "blue",
    comparison.name = "Comparison", comparison.color = "red", ...
)
{

  # Set up plot
  return(
    silent_plot(
      ggplot2::ggplot(x$TEFI.df, ggplot2::aes(x = TEFI, fill = Structure)) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::geom_vline(
          ggplot2::aes(xintercept = mean(x$TEFI.df$TEFI[x$TEFI.df$Structure == "Base"])),
          color = base.color, linetype = "dashed", size = 1
        ) +
        ggplot2::geom_vline(
          ggplot2::aes(xintercept = mean(x$TEFI.df$TEFI[x$TEFI.df$Structure == "Comparison"])),
          color = comparison.color, linetype = "dashed", size = 1
        ) +
        ggplot2::scale_fill_manual(
          labels = c(base.name, comparison.name),
          values = c(base.color, comparison.color)
        ) +
        ggplot2::labs(
          title = "TEFI Comparison",
          subtitle = paste("p-value:", round(x$p.value, 4)),
          x = "TEFI", y = "Density", fill = "Structure"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 16),
          plot.subtitle = ggplot2::element_text(size = 14),
          axis.title = ggplot2::element_text(size = 14),
          axis.text = ggplot2::element_text(size = 12),
          axis.text.y = ggplot2::element_blank(),
          legend.title = ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 12),
          legend.position = "bottom"
        )
    )
  )

}

#' @noRd
# Global variables needed for CRAN checks ----
# Updated 28.01.2025
utils::globalVariables(c("TEFI", "Structure"))
