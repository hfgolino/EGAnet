#' @title S3 Plot Methods for \code{\link{EGAnet}}
#' 
#' @name EGAnet-plot
#' 
#' @description General usage for plots created by \code{\link{EGAnet}}'s S3 methods.
#' Plots across the \code{\link{EGAnet}} package leverage \code{\link[GGally]{ggnet2}} 
#' and \code{\link[ggplot2]{ggplot}}.
#' 
#' Most plots allow the full usage of the \code{gg*} series functionality and therefore
#' plotting arguments should be referenced through those packages rather than here in
#' \code{\link{EGAnet}}.
#' 
#' The sections below list the functions and their usage for the S3 plot methods.
#' The plot methods are intended to be generic and without many arguments so that
#' nearly all arguments are passed to \code{\link[GGally]{ggnet2}} and \code{\link[ggplot2]{ggplot}}.
#' 
#' There are some constraints placed on certain plots to keep the \code{\link{EGAnet}} style
#' throughout the (network) plots in the package, so be aware that if some settings are
#' not changing your plot output, then these settings might be fixed 
#' to maintain the \code{\link{EGAnet}} style
#' 
#' \emph{Do \code{\link{EGAnet}} plots load slow? \cr
#' Check out our Wiki on faster plots in 
#' \href{https://tinyurl.com/EGAnet-plotting}{EGAnet (and R)}}
#' 
#' @usage 
#' plot(x, ...)
#' plot.dynEGA(x, base = 1, id = NULL, ...)
#' plot.dynEGA.Group(x, base = 1, ...)
#' plot.dynEGA.Individual(x, base = 1, id = NULL, ...)
#' plot.hierEGA(
#'   x, plot.type = c("multilevel", "separate"), 
#'   color.match = FALSE, ...
#' )
#' plot.invariance(
#'   x, p_type = c("p", "p_BH"), p_value = 0.05, ...
#' )
#' 
#' @param x \code{\link{EGAnet}} object with available S3 plot
#' method (see full list below)
#' 
#' @param color.palette Character (vector).
#' Either a character (length = 1) from the
#' pre-defined palettes in \code{\link[EGAnet]{color_palette_EGA}}
#' or character (length = total number of communities) using
#' HEX codes (see \strong{Color Palettes} and \strong{Examples} sections)
#' 
#' @param layout Character (length = 1).
#' Layouts can be set using \code{\link[sna]{gplot.layout}} and the
#' ending layout name; for example, \code{gplot.layout.circle} can be set
#' in these functions using \code{layout = "circle"} or \code{mode = "circle"}
#' (see \strong{Examples})
#' 
#' @param base Numeric (length = 1).
#' Plot to be used as the base for the configuration of the networks.
#' Uses the number of the order in which the plots are input.
#' Defaults to \code{1} or the first plot
#' 
#' @param id Numeric index(es) or character name(s).
#' IDs to use when plotting \code{\link[EGAnet]{dynEGA}} \code{level = "individual"}.
#' Defaults to \code{NULL} or 4 IDs drawn at random
#' 
#' @param plot.type Character (length = 1).
#' Whether \code{\link[EGAnet]{hierEGA}} networks should plotted in
#' a stacked, \code{"multilevel"} fashion or as \code{"separate"} plots.
#' Defaults to \code{"multilevel"}
#' 
#' @param color.match Boolean (length = 1).
#' Whether lower order community colors in the \code{\link[EGAnet]{hierEGA}} plot
#' should be "matched" and used as the border color for the higher order
#' communities.
#' Defaults to \code{FALSE}
#' 
#' @param p_type Character (length = 1).
#' Type of \emph{p}-value when plotting \code{\link[EGAnet]{invariance}}.
#' Defaults to \code{"p"} or uncorrected \emph{p}-value.
#' Set to \code{"p_BH"} for the Benjamini-Hochberg corrected \emph{p}-value
#' 
#' @param p_value Numeric (length = 1).
#' The \emph{p}-value to use alongside \code{p_type} when 
#' plotting \code{\link[EGAnet]{invariance}}.
#' Defaults to \code{0.05}
#' 
#' @param ... Additional arguments to pass on to
#' \code{\link[GGally]{ggnet2}} and \code{\link[sna]{gplot.layout}}
#' (see \strong{Examples})
#' 
#' @section \code{*EGA} Plots:
#' 
#' \code{\link[EGAnet]{bootEGA}}, \code{\link[EGAnet]{dynEGA}},
#' \code{\link[EGAnet]{EGA}}, \code{\link[EGAnet]{EGA.estimate}},
#' \code{\link[EGAnet]{EGA.fit}}, \code{\link[EGAnet]{hierEGA}},
#' \code{\link[EGAnet]{invariance}}, \code{\link[EGAnet]{riEGA}}
#' 
#' @section All Available S3 Plot Methods:
#' 
#' \code{\link[EGAnet]{boot.ergoInfo}}, \code{\link[EGAnet]{bootEGA}}, 
#' \code{\link[EGAnet]{dynEGA}}, \code{dynEGA.Group}, \code{dynEGA.Individual},
#' \code{dynEGA.Population}, \code{\link[EGAnet]{EGA}}, 
#' \code{\link[EGAnet]{EGA.estimate}}, \code{\link[EGAnet]{EGA.fit}}, 
#' \code{\link[EGAnet]{hierEGA}}, \code{\link[EGAnet]{infoCluster}}, 
#' \code{\link[EGAnet]{invariance}}, \code{\link[EGAnet]{itemStability}},
#' \code{\link[EGAnet]{riEGA}}
#' 
#' @section Color Palettes:
#' 
#' \code{\link[EGAnet]{color_palette_EGA}} will implement some color palettes in
#' \code{\link{EGAnet}}. The main \code{\link{EGAnet}} style palette is \code{"polychrome"}. 
#' This palette currently has 40 colors but there will likely be a need to expand it further
#' (e.g., \code{\link[EGAnet]{hierEGA}} demands a lot of colors).
#' 
#' The \code{color.palette} argument will also accept HEX code colors that 
#' are the same length as the number of communities in the plot.
#' 
#' In any network plots, the \code{color.palette} argument can be used to
#' select color palettes from \code{\link[EGAnet]{color_palette_EGA}} as well
#' as those in the color scheme of \code{\link[RColorBrewer]{RColorBrewer}}
#' 
#' @examples
#' \dontrun{
#' # Using different arguments in {GGally}'s `ggnet2`
#' plot(ega.wmt, node.size = 6, edge.size = 4)
#' 
#' # Using a different layout in {sna}'s `gplot.layout`
#' plot(ega.wmt, layout = "circle") # 'layout' argument
#' plot(ega.wmt, mode = "circle") # 'mode' argument
#' 
#' # Using different color palettes with `color_palette_EGA`
#' 
#' ## Pre-defined palette
#' plot(ega.wmt, color.palette = "blue.ridge2")
#' 
#' ## University of Virginia colors
#' plot(ega.wmt, color.palette = c("#232D4B", "#F84C1E"))
#' 
#' ## Vanderbilt University colors
#' ## (with additional {GGally} `ggnet2` argument)
#' plot(
#'   ega.wmt, color.palette = c("#FFFFFF", "#866D4B"), 
#'   label.color = "#000000"
#' )}
#' 
#' @aliases plot.EGAnet
#' 
#' @rdname EGAnet-plot
#' 
# Updated 04.08.2023
NULL