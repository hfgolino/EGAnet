% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color_palette_EGA.R
\name{color_palette_EGA}
\alias{color_palette_EGA}
\title{\code{\link[EGAnet]{EGA}} Color Palettes}
\usage{
color_palette_EGA(
  name = c("polychrome", "blue.ridge1", "blue.ridge2", "rainbow", "rio", "itacare",
    "grayscale"),
  wc,
  sorted = FALSE
)
}
\arguments{
\item{name}{Character.
Name of color scheme (see \code{\link[RColorBrewer]{RColorBrewer}}).
Defaults to \code{"polychrome"}.
\code{\link[EGAnet]{EGA}} palettes:

\itemize{

\item \code{"polychrome"} --- Default 40 color palette

\item \code{"grayscale"} --- "grayscale", "greyscale", or "colorblind" will produce
plots suitable for publication purposes

\item \code{"blue.ridge1"} --- Palette inspired by the Blue Ridge Mountains

\item \code{"blue.ridge2"} --- Second palette inspired by the Blue Ridge Mountains

\item \code{"rainbow"} --- Rainbow colors. Default for \code{\link[qgraph]{qgraph}}

\item \code{"rio"} --- Palette inspired by Rio de Janiero, Brazil

 \item \code{"itacare"} --- Palette inspired by Itacare, Brazil

}

For custom colors, enter HEX codes for each dimension in a vector}

\item{wc}{Numeric vector.
A vector representing the community (dimension) membership
of each node in the network. \code{NA} values mean that the node
was disconnected from the network}

\item{sorted}{Boolean.
Should colors be sorted by \code{wc}?
Defaults to \code{FALSE}}
}
\value{
Vector of colors for community memberships
}
\description{
Color palettes for plotting \code{\link[GGally]{ggnet2}}
\code{\link[EGAnet]{EGA}} network plots
}
\examples{
# Default
color_palette_EGA(name = "polychrome", wc = ega.wmt$wc)

# Blue Ridge Moutains 1
color_palette_EGA(name = "blue.ridge1", wc = ega.wmt$wc)

# Custom
color_palette_EGA(name = c("#7FD1B9", "#24547e"), wc = ega.wmt$wc)

}
\seealso{
\code{\link[EGAnet]{plot.EGAnet}} for plot usage in \code{\link{EGAnet}}
}
\author{
Hudson Golino <hfg9s at virginia.edu>, Alexander P. Christensen <alexpaulchristensen at gmail.com>
}
