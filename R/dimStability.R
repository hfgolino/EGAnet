#' Dimension Stability Statistics from \code{\link[EGAnet]{bootEGA}}
#'
#' @description Based on the \code{\link[EGAnet]{bootEGA}} results, this function
#' computes the stability of dimensions. This is computed by assessing the proportion of
#' times the original dimension is exactly replicated in across bootstrap samples
#'
#' @param bootega.obj A \code{\link[EGAnet]{bootEGA}} object
#'
#' @param orig.wc Numeric or character.
#' A vector with community numbers or labels for each item.
#' Typically uses community results (\code{wc}) from \code{\link[EGAnet]{EGA}}
#'
#' @param item.stability Boolean.
#' Should the item stability statistics be computed
#' using \code{[EGAnet]{itemStability}}?
#' Defaults to \code{TRUE}
#'
#' @return When argument \code{item.stability = TRUE}, returns a list containing:
#'
#' \item{dimensions}{The dimensional stability of each dimension}
#'
#' \item{items}{The output from \code{[EGAnet]{itemStability}}}
#'
#' When argument \code{item.stability = FALSE}, returns a vector of the
#' dimensional stability of each dimension
#'
#' @examples
#'
#' # Load data
#' wmt <- wmt2[,7:24]
#'
#' \dontrun{
#' # Estimate EGA network
#' ega.wmt <- EGA(data = wmt, model = "glasso")
#'
#' # Estimate dimension stability
#' boot.wmt <- bootEGA(data = wmt, n = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "glasso",
#' type = "parametric", ncores = 4)
#'
#' # Estimate item stability statistics
#' dimStability(boot.wmt, orig.wc = ega.wmt$wc, item.stability = FALSE)
#' }
#'
#' @seealso \code{\link[EGAnet]{EGA}} to estimate the number of dimensions of an instrument using EGA and
#' \code{\link[EGAnet]{CFA}} to verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#'
#' @author Hudson F. Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#' 
# Dimension Stability function
# Updated 22.05.2020
dimStability <- function(bootega.obj, orig.wc, item.stability = TRUE)
{
  if(class(bootega.obj) != "bootEGA")
  {stop("Input for 'bootega.obj' is not a 'bootEGA' object")}

  # Compute item stability
  items <- itemStability(bootega.obj, orig.wc, item.freq = 0, plot.item.rep = item.stability)

  # Compute dimension stability
  ## Grab dimensions from itemStability output
  dims <- items$wc
  ## Identify unique dimensions
  uniq.dim <- items$uniq.name
  ## Idetify unique dimension numbers
  uniq.num <- items$uniq.num
  ## Number of dimensions
  dim.len <- length(uniq.dim)
  ## Initialize dimension stability vector
  dim.stab <- numeric(dim.len)
  names(dim.stab) <- uniq.dim

  # Loop through dims
  for(i in 1:dim.len)
  {
    # Target items
    target <- which(orig.wc == uniq.dim[i])

    # Initialize count vector
    dim.count <- numeric(length = ncol(dims))

    # Identify consistency across bootstraps
    for(j in 1:ncol(dims))
    {dim.count[j] <- all(dims[target,j] == uniq.num[i])}

    # Input mean of into vector
    dim.stab[i] <- mean(dim.count, na.rm = TRUE)
  }

  # Results list
  if(item.stability)
  {
    res <- list()
    res$dimensions <- round(dim.stab[order(names(dim.stab))],3)
    res$items <- items
  }else{res <- round(dim.stab[order(names(dim.stab))],3)}

  return(res)
}
#----
