% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plots}
\alias{plots}
\alias{plot.bootEGA}
\alias{plot.CFA}
\alias{plot.dynEGA}
\alias{plot.dynEGA.Groups}
\alias{plot.dynEGA.Individuals}
\alias{plot.EGA}
\alias{plot.NetLoads}
\title{S3Methods for Plotting}
\usage{
\method{plot}{bootEGA}(x, plot.type = c("GGally", "qgraph"), plot.args = list(), ...)

\method{plot}{CFA}(x, layout = "spring", vsize = 6, ...)

\method{plot}{dynEGA}(x, title = "", plot.type = c("GGally", "qgraph"), plot.args = list(), ...)

\method{plot}{dynEGA.Groups}(x, ncol, nrow, title = "",
plot.type = c("GGally", "qgraph"), plot.args = list(),  ...)

\method{plot}{dynEGA.Individuals}(x, title = "",  id = NULL,
plot.type = c("GGally", "qgraph"), plot.args = list(), ...)

\method{plot}{EGA}(x, title = "", plot.type = c("GGally", "qgraph"), plot.args = list(), ...)

\method{plot}{NetLoads}(x, ...)
}
\arguments{
\item{x}{Object from \code{EGAnet} package}

\item{plot.type}{Character.
Plot system to use.
Current options are \code{\link[qgraph]{qgraph}} and \code{\link{GGally}}.
Defaults to \code{"GGally"}.}

\item{plot.args}{List.
A list of additional arguments for the network plot.
For \code{plot.type = "qgraph"}:

\itemize{

\item{\strong{\code{vsize}}}
{Size of the nodes. Defaults to 6.}

}
For \code{plot.type = "GGally"} (see \code{\link[GGally]{ggnet2}} for
full list of arguments):

\itemize{

\item{\strong{\code{vsize}}}
{Size of the nodes. Defaults to 6.}

\item{\strong{\code{label.size}}}
{Size of the labels. Defaults to 5.}

\item{\strong{\code{alpha}}}
{The level of transparency of the nodes, which might be a single value or a vector of values. Defaults to 0.7.}

\item{\strong{\code{edge.alpha}}}
{The level of transparency of the edges, which might be a single value or a vector of values. Defaults to 0.4.}

\item{\strong{\code{legend.names}}}
{A vector with names for each dimension}

\item{\strong{\code{color.palette}}}
{The color palette for the nodes. For custom colors,
enter HEX codes for each dimension in a vector.
See \code{\link[EGAnet]{color_palette_EGA}} for 
more details and examples}

}}

\item{...}{Arguments passed on to

\itemize{

\item{\code{\link[qgraph]{qgraph}}}
{Functions: bootEGA, dynEGA, dynEGA.Groups, dynEGA.Individuals, EGA, and net.loads}

\item{\code{\link[semPlot]{semPaths}}}
{Functions: CFA}

}}

\item{vsize}{Numeric.
Size of vertices in \code{\link[EGAnet]{CFA}} plots.
Defaults to \code{6}}

\item{layout}{Character.
Layout of plot (see \code{\link[semPlot]{semPaths}}).
Defaults to "spring"}

\item{ncol}{Numeric.
Number of columns}

\item{nrow}{Numeric.
Number of rows}

\item{title}{Character.
Title of the plot.
Defaults to \code{""}}

\item{id}{Numeric.
An integer or character indicating the ID of the individual to plot}
}
\value{
Plots of \code{EGAnet} object
}
\description{
Plots for \code{EGAnet} objects
}
\author{
Hudson Golino <hfg9s at virginia.edu> and Alexander P. Christensen <alexpaulchristensen@gmail.com>
}
