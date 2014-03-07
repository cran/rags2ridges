\name{sparsify}
\alias{sparsify}
\title{
Determine the support of a partial correlation/precision matrix
}
\description{
Function that determines the support of a partial correlation/precision matrix and sparsifies it accordingly.
}
\usage{
sparsify(P, type = c("threshold", "localFDR"), threshold = 0.25, FDRcut = 0.8)
}
\arguments{
  \item{P}{
(Possibly regularized) precision \code{matrix}.
}
  \item{type}{
A \code{character} signifying type of sparsification. Must be one of: "threshold", "localFDR".
}
  \item{threshold}{
A \code{numeric} giving the cut-off for partial correlation element selection based on thresholding. 
}
  \item{FDRcut}{
A \code{numeric} giving the cut-off for partial correlation element selection based on local false discovery rates (FDR).
}
}
\details{
The function transforms the possibly regularized input precision matrix to a partial correlation matrix. Subsequently,
the support of this partial correlation matrix is determined. Support determination is performed either by simple
thresholding on the absolute values of matrix entries (\code{type = "threshold"}) or by usage of local FDR (\code{type = "localFDR"}). 
The argument \code{threshold} is only used when \code{type = "threshold"}. The argument \code{FDRcut} is only used when
\code{type = "localFDR"}. The function is to some extent a wrapper around certain \code{GeneNet} and \code{fdrtool} functions when
\code{type = "localFDR"}. In that case a mixture model is fitted to the nonredundant partial correlations by \code{fdrtool}. The decision to
retain elements is then based on the argument \code{FDRcut}. Elements with a posterior probability \eqn{\geq} FDRcut (equalling
1 - local FDR) are retained. See Schaefer and Strimmer (2005) for further details on usage of local FDR in graphical modeling.
}
\value{
A sparsified partial correlation \code{matrix}. 
}
\references{
Schaefer, J., and Strimmer, K. (2005). A shrinkage approach to large-scale covariance estimation
and implications for functional genomics. Statistical Applications in Genetics and Molecular Biology, 4:32.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeS}}, \code{\link{optPenaltyCV}}, \code{\link{symm}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- cov(X)

## Obtain regularized precision under optimal penalty
OPT <- optPenaltyCV(X, 15, 30, 20, output = "light")

## Determine support regularized (standardized) precision under optimal penalty
sparsify(symm(OPT$optPrec), type = "localFDR")
}