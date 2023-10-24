#' @title \code{\link[EGAnet]{EGA}} Color Palettes
#'
#' @description Color palettes for plotting \code{\link[GGally]{ggnet2}} 
#' \code{\link[EGAnet]{EGA}} network plots
#'
#' @param name Character.
#' Name of color scheme (see \code{\link[RColorBrewer]{RColorBrewer}}).
#' Defaults to \code{"polychrome"}.
#' \code{\link[EGAnet]{EGA}} palettes:
#'
#' \itemize{
#'
#' \item \code{"polychrome"} --- Default 40 color palette
#' 
#' \item \code{"grayscale"} --- "grayscale", "greyscale", or "colorblind" will produce
#' plots suitable for publication purposes
#'
#' \item \code{"blue.ridge1"} --- Palette inspired by the Blue Ridge Mountains
#'
#' \item \code{"blue.ridge2"} --- Second palette inspired by the Blue Ridge Mountains
#' 
#' \item \code{"rainbow"} --- Rainbow colors. Default for \code{\link[qgraph]{qgraph}}
#'
#' \item \code{"rio"} --- Palette inspired by Rio de Janiero, Brazil
#'
#'  \item \code{"itacare"} --- Palette inspired by Itacare, Brazil
#'
#' }
#'
#' For custom colors, enter HEX codes for each dimension in a vector
#'
#' @param wc Numeric vector.
#' A vector representing the community (dimension) membership
#' of each node in the network. \code{NA} values mean that the node
#' was disconnected from the network
#' 
#' @param sorted Boolean.
#' Should colors be sorted by \code{wc}?
#' Defaults to \code{FALSE}
#'
#' @author Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen at gmail.com>
#'
#' @return Vector of colors for community memberships
#'
#' @examples
#' # Default
#' color_palette_EGA(name = "polychrome", wc = ega.wmt$wc)
#'
#' # Blue Ridge Moutains 1
#' color_palette_EGA(name = "blue.ridge1", wc = ega.wmt$wc)
#'
#' # Custom
#' color_palette_EGA(name = c("#7FD1B9", "#24547e"), wc = ega.wmt$wc)
#' 
#' @seealso \code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
#'
#' @export
#' 
# Color palettes for EGA ----
# Updated 24.10.2023
color_palette_EGA <- function(
    name = c(
      "polychrome", "blue.ridge1", "blue.ridge2",
      "rainbow", "rio", "itacare", "grayscale"
    ), wc, sorted = FALSE
)
{

  # Set default for name
  if(missing(name)){name <- "polychrome"}
  
  # Get lowercase name
  lower_name <- tolower(name)
  
  # Get length of name
  length_name <- length(name)
  
  # Grayscale palette
  if(length_name == 1 && lower_name %in% c("greyscale", "grayscale", "colorblind")){
    name <- "grayscale"
  }

  # Set default for memberships
  if(missing(wc)){wc <- seq_along(name)}
  
  # Ensure that memberships are numeric
  wc <- as.numeric(factor(wc))

  # Get length of memberships
  length_wc <- length(wc)
  
  # Unique memberships
  unique_wc <- unique(wc)
  
  # Get number of communities
  communities <- length(unique_wc) - anyNA(unique_wc)
  
  # Send length error
  length_error(name, c(1, communities, length_wc))
  
  # Sort memberships
  if(isTRUE(sorted)){
    wc <- sort(wc, na.last = TRUE)
  }
  
  # Determine length of names
  if(length_name == 1){
    
    # Get sequence of memberships
    sequence_wc <- seq_len(communities)
    
    # Check if name is in color brewer
    if(name %in% row.names(RColorBrewer::brewer.pal.info)){
      return(silent_call(RColorBrewer::brewer.pal(communities, name)[wc]))
    }else if(
      lower_name %in% c(
        "polychrome", "blue.ridge1", "blue.ridge2",
        "rainbow", "rio", "itacare", "grayscale"
      ) 
    ){
      
      # Perform switch on palettes
      return(
        switch(
          lower_name,
          "polychrome" = polychrome(communities),
          "blue.ridge1" = blue_ridge1(),
          "blue.ridge2" = blue_ridge2(),
          "rainbow" = rainbow(communities),
          "rio" = rio(),
          "itacare" = itacare(),
          "grayscale" = grayscale(communities)
        )[wc]
      )
      
    }else{
      
      # Return single value for all
      return(rep(name, length_wc))
      
    }
    
  }
  
  # Otherwise, ensure hashtags and return
  name <- cvapply(name, function(x){
    swiftelse(substr(x, 1, 1) == "#", x, paste0("#", x))
  })
  
  # Check for whether colors are same length as memberships
  return(
    swiftelse(
      length_wc != length_name,
      name[wc], name
    )
  )

}

#' @noRd
# Polychrome color palette ----
# Updated 26.07.2023
polychrome <- function(communities)
{
  
  # Check for total colors
  if(communities <= 40){
    
    # Get {EGAnet} polychrome
    polychrome <- c(
      "#F03D2D", "#90DDF0", "#C8D96F", "#ef8a17", "#f5c900",
      "#ba42c0", "#17BEBB", "#9bafd9", "#f27a7d", "#f9c58d",
      "#f7f779", "#c5f9d7", "#a18dce", "#f492f0", "#919bff",
      "#c792df", "#ff4f79", "#f9a470", "#bc556f", "#f7a2a1",
      "#3a445d", "#5e5768", "#928779", "#d4d2a5", "#a11692",
      "#6d1a36", "#bce7fd", "#53917e", "#dd99bb", "#fcd0a1",
      "#a9fbd7", "#d81159", "#8f2d56", "#006ba6", "#39304a",
      "#ff470a", "#60b6f1", "#fcdebe", "#cff27e", "#b87d4b"
    )
    
  }else{
    
    # Get distinct colors from color brewer
    qual_col_pals <- RColorBrewer::brewer.pal.info[
      RColorBrewer::brewer.pal.info$category == "qual",
    ]
    # Obtain unique distinct colors
    polychrome <- unique(
      unlist(
        mapply(
          RColorBrewer::brewer.pal,
          qual_col_pals$maxcolors,
          rownames(qual_col_pals)
        )
      )
    )
    
  }
  
  # Return polychrome
  return(polychrome)
  
}

#' @noRd
# Blue ridge 1 color palette ----
# 7 colors
# Updated 26.07.2023
blue_ridge1 <- function()
{
  
  # Return colors
  return(
    c(
      "#272a39", "#24547e", "#4c6e98", "#7f616e",
      "#fdb184", "#fde8a9", "#fdcd9b"
    )
  )
  
}

#' @noRd
# Blue ridge 2 color palette ----
# 10 colors
# Updated 26.07.2023
blue_ridge2 <- function()
{
  
  # Return colors
  return(
    blue.ridge2 <- c(
      "#26405b", "#facf92", "#497397",
      "#8c7f8e", "#a5a9a9", "#68788b",
      "#e2a187", "#e3ccb5", "#c48480", "#fcac6c"
    )
  )
  
}

#' @noRd
# Rainbow color palette ----
# Updated 26.07.2023
rainbow <- function(communities)
{
  
  # Return colors
  return(
    grDevices::rainbow(communities)
  )
  
}

#' @noRd
# Rio color palette ----
# 10 colors
# Updated 26.07.2023
rio <- function()
{
  
  # Return colors
  return(
    c(
      "#fac9af", "#a95c5b", "#322a30",
      "#654145", "#f09b5f", "#985e36",
      "#ea897c", "#9c8062", "#524954", "#54544c"
    )
  )
  
}

#' @noRd
# Itacare color palette ----
# 10 colors
# Updated 26.07.2023
itacare <- function()
{
  
  # Return colors
  return(
    c(
      "#232b17", "#cbbda4", "#2888ab",
      "#0581c9", "#7e8056", "#d9e7e6",
      "#8ec0c5", "#a58a60", "#ad9342", "#a96c2e"
    )
  )
  
}

#' @noRd
# Grayscale color palette ----
# 16 colors
# Updated 26.07.2023
grayscale <- function(communities)
{
  
  # Grayscale colors
  colors <- c(
    "#F0F0F0", "#E9E9E9", "#E1E1E1", "#C2C2C2",
    "#A4A4A4", "#959595", "#8D8D8D", "#858585",
    "#767676", "#6F6F6F", "#606060", "#565656",
    "#505050", "#484848", "#454545", "#333333"
  )
  
  # Get distance
  distance <- round(16 / communities)
  
  # Check for zero distance
  if(distance == 0){
    distance <- 1
  }
  
  # Return colors
  return(colors[seq.int(1, 16, distance)])
  
}

