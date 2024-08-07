% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frobenius.R
\name{frobenius}
\alias{frobenius}
\title{Frobenius Norm (Similarity)}
\usage{
frobenius(network1, network2)
}
\arguments{
\item{network1}{Matrix or data frame.
Network to be compared}

\item{network2}{Matrix or data frame.
Second network to be compared}
}
\value{
Returns Frobenius Norm
}
\description{
Computes the Frobenius Norm (Ulitzsch et al., 2023)
}
\examples{
# Obtain wmt2 data
wmt <- wmt2[,7:24]

# Set seed (for reproducibility)
set.seed(1234)

# Split data
split1 <- sample(
  1:nrow(wmt), floor(nrow(wmt) / 2)
)
split2 <- setdiff(1:nrow(wmt), split1)

# Obtain split data
data1 <- wmt[split1,]
data2 <- wmt[split2,]

# Perform EBICglasso
glas1 <- EBICglasso.qgraph(data1)
glas2 <- EBICglasso.qgraph(data2)

# Frobenius norm
frobenius(glas1, glas2)
# 0.7070395

}
\references{
\strong{Simulation Study} \cr
Ulitzsch, E., Khanna, S., Rhemtulla, M., & Domingue, B. W. (2023).
A graph theory based similarity metric enables comparison of subpopulation psychometric networks
\emph{Psychological Methods}.
}
\author{
Hudson Golino <hfg9s at virginia.edu> & Alexander P. Christensen <alexander.christensen at Vanderbilt.Edu>
}
