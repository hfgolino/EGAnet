#'  Surrogate analysis for EGA.
#'
#' \code{surrogate} Compares the median number of dimensions (and it's 95% confidence intervals) estimated in n shuffled versions of the
#' original dataset and in n random subsamples (\code{\link{subsamples}}) of the data.
#' @param data A data.frame object.
#' @param n Numeric. Number of estimates.
#' @param ncores Numeric. Number of cores to use in parallel.
#' @param plot.surrogate Logic. If true, returns a plot with the surrogate analysis results.
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#'\dontrun{
#' surrogate.analysis <- surrogate(data = wmt2[,7:24], n = 100, ncores = 4)
#' }
#' @seealso \code{\link{subsamples}} to estimate the number of dimensions via EGA in n random subsamples of the original data and \code{\link{shuffle}}
#' to generate n estimations of the number of dimensions in shuffled versions of the original dataset
#' 
#' @importFrom stats rnorm
#' @importFrom rlang .data
#' 
#' @export
#'
## Surrogate EGA:

surrogate <- function(data, n, ncores, plot.surrogate = TRUE){

# Resampling, Surrogate and Normal Distribution:

  subsample <- subsamples(data, n, ncores)
  surrogate <- shuffle(data, n, ncores)
  normal.dist <- rnorm(n, mean = median(surrogate), sd = sd(surrogate))

# Computing the Median, SE and CI of the Median:

  distributions <- data.frame(Observed = subsample, Surrogate = surrogate, Normal.Dist = normal.dist)
  distributions.gather <- tidyr::gather(distributions, Dist)
  distributions.gather$Dist <- factor(distributions.gather$Dist, levels = c("Observed", "Surrogate", "Normal.Dist"))

  ciMult <- qt(0.95/2 + 0.5, n - 1)

  summary.distribution <- as.data.frame(distributions.gather %>%
                                          dplyr::group_by(Dist) %>%
                                          dplyr::summarise(
                                            Median <- median(value),
                                            SD <- sd(value),
                                            SE <- (1.253 * SD)/sqrt(1000)) %>%
                                          dplyr::mutate(
                                              CI <- SE * ciMult) %>%
                                          dplyr::mutate(
                                              Lower <- Median - CI,
                                              Upper <- Median + CI))

  sig1 <- ifelse(summary.distribution[1,7] >= summary.distribution[2,6] | summary.distribution[1,6] >= summary.distribution[2,7], "Non-Sig", "Sig")
  sig2 <- ifelse(summary.distribution[1,7] >= summary.distribution[3,6] | summary.distribution[1,6] >= summary.distribution[3,7], "Non-Sig", "Sig")

  summary.distribution$Significance <- NA
  summary.distribution$Significance <- c(" ", sig1, sig2)

  # Plot:
  if(plot.surrogate == TRUE) {
    g <- ggplot2::ggplot(summary.distribution, ggplot2::aes(x=Dist, y=Median, group=1)) +
      ggplot2::geom_point(ggplot2::aes(size=1000), alpha=0.52) +
      ggplot2::geom_errorbar(width=.1, ggplot2::aes(ymin=Lower, ymax=Upper), colour="darkred") +
      ggplot2::labs(x="Type",y= "Number of Dimensions", title="Mean and Confidence Intervals of the Number of Dimensions Assessed via EGA", subtitle="Observed (Random Subsample), Shuffling (Surrogate EGA) and Normal Distribution \nError Bar Using median as center with 95% Confidence Intervals", caption="Hudson Golino") +
      ggthemes::theme_fivethirtyeight()
    print(g)
  }
  results <- vector("list")
  results$summary.distribution <- summary.distribution
  results$subsample <- subsample
  results$shuffled <- surrogate
  results$normal.dist <- normal.dist
  return(summary.distribution)
}

