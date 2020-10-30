#' @title Shiny App for \code{\link{EGAnet}}
#' 
#' @description An interactive Shiny application for running \code{\link{EGAnet}} analysis
#' 
#' @return A list called \code{resultShiny} containing:
#' 
#' \item{data}{The data imported into \code{\link[EGAnet]{shinyEGA}}}
#' 
#' \item{EGA}{The output generated from \code{\link[EGAnet]{EGA}}}
#' 
#' @examples
#' if(interactive())
#' {shinyEGA()}
#' 
#' @author Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
#' 
#' @export
# EGAnet Shiny App----
# Updated 30.10.2020
shinyEGA <- function(){
  shiny::runApp(appDir = system.file("Shiny", package="EGAnet"))
}