% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert2igraph.R
\name{convert2igraph}
\alias{convert2igraph}
\title{Convert networks to \code{igraph}}
\usage{
convert2igraph(A, diagonal = 0)
}
\arguments{
\item{A}{Matrix or data frame.
\emph{N} x \emph{N} matrix where \emph{N} is the number of nodes}

\item{diagonal}{Numeric.
Value to be placed on the diagonal of \code{A}.
Defaults to \code{0}}
}
\value{
Returns a network in the \code{igraph} format
}
\description{
Converts networks to \code{igraph} format
}
\examples{
convert2igraph(ega.wmt$network)

}
\author{
Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
}
