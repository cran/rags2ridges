\name{GGMnetworkStats}
\alias{GGMnetworkStats}
\title{
Gaussian graphical model network statistics
}
\description{
Function that calculates various network statistics from a sparse precision matrix. The sparse precision matrix
is taken to represent the conditional indepence graph of a Gaussian graphical model. 
}
\usage{
GGMnetworkStats(sparseP)
}
\arguments{
  \item{sparseP}{
Sparse precision/partial correlation \code{matrix}.
}
}
\details{
The function calculates various network statistics from a sparse matrix. The input matrix \code{P} is assumed to 
be a sparse precision or partial correlation matrix. The sparse matrix is taken to represent a conditional independence graph.
In the Gaussian setting, conditional independence corresponds to zero entries in the (standardized) precision matrix. Each node in
the graph represents a Gaussian variable, and each undirected edge represents conditional dependence in the sense of a 
nonzero corresponding precision entry.

The function calculates various measures of centrality: node degree, betweenness centrality and closeness centrality. It also
calculates the number of positive and the number of negative edges for each node. In addition, for each variate the mutual information
(with all other variates), the variance, and the partial variance is represented. For more information on network measures, consult, e.g.,
Newman (2010).
}
\value{
An object of class \code{list}:
\item{degree}{A \code{numeric} vector with the node degree for each node.}
\item{betweenness}{A \code{numeric} vector representing the betweenness centrality for each node.}
\item{closeness}{A \code{numeric} vector representing the closeness centrality for each node.}
\item{nNeg}{An \code{integer} vector representing the number of negative edges for each node.}
\item{nPos}{An \code{integer} vector representing the number of positive edges for each node.}
\item{mutualInfo}{A \code{data.frame} with the mutual information (with all other nodes) for each node.}
\item{variance}{A \code{numeric} vector representing the variance of each node.}
\item{partialVariance}{A \code{numeric} vector representing the partial variance of each node.}
}
\references{
Newman, M.E.J. (2010). "Networks: an introduction", Oxford University Press.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeS}}, \code{\link{covML}}, \code{\link{sparsify}}, \code{\link{Ugraph}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Obtain sparsified partial correlation matrix
Pridge   <- ridgeS(Cx, 10, type = "Alt")
PCsparse <- sparsify(Pridge , threshold = "localFDR")$sparseParCor

## Represent the graph and calculate GGM network statistics
Ugraph(PCsparse, "fancy")
GGMnetworkStats(PCsparse)
}