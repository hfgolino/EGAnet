#'  Plot the replicability of EGA item-wise.
#'
#' \code{itemConfirm} Based on the bootEGA results, this function plots the number of times an item (variable) is estimated in the same
#' factor/dimension as originaly estimated by EGA.
#'
#' @param bootega.obj A bootEGA object
#' @author Hudson F. Golino <hfg9s at virginia.edu>
#' @examples
#' \dontrun{
#' ega.wmt <- EGA(data = wmt2[,7:24], model = "glasso")
#' boot.wmt <- bootEGA(data = wmt2[,7:24], n = 100, typicalStructure = TRUE,
#' plot.typicalStructure = TRUE, model = "GGM",
#' type = "parametric", ncores = 4, confirm = ega.wmt$wc)
#' itemConfirm(boot.wmt)
#' }
#'
#' \dontrun{
#' itemConfirm(bootEGA)
#' }
#' @seealso \code{\link{EGA}} to estimate the number of dimensions of an instrument using EGA and \code{\link{CFA}} to
#' verify the fit of the structure suggested by EGA using confirmatory factor analysis.
#' @export

itemConfirm <- function(bootega.obj){
    require(ggpubr)
    
    comm <- bootega.obj$orig.wc
    rain <- grDevices::rainbow(max(comm))
    
    item.rep <- data.frame(Item = names(bootega.obj$item.confirm),
                           Rep = bootega.obj$item.confirm,
                           Comm = factor(comm))
    
    
    ggdotchart(item.rep, x = "Item", y = "Rep",
               group = "Comm", color = "Comm",
               palette = rain,
               legend.title = "EGA Communities",
               sorting = NULL,
               add = "segments",
               rotate = TRUE,
               dot.size = 6,
               label = round(item.rep$Rep, 2),
               font.label = list(color = "black", size = 8,
                                 vjust = 0.5),
               ggtheme = theme_pubr()
    )
}
