\name{sparsify}
\alias{sparsify}
\title{
Determine the support of a partial correlation/precision matrix
}
\description{
Function that determines the support of a partial correlation/precision matrix by thresholding and sparsifies it accordingly.
}
\usage{
sparsify(P, threshold = c("absValue", "localFDR"), absValueCut = 0.25, 
FDRcut = 0.8, verbose = TRUE)
}
\arguments{
  \item{P}{
(Possibly regularized) precision \code{matrix}.
}
  \item{threshold}{
A \code{character} signifying type of sparsification by thresholding. Must be one of: "absValue", "localFDR".
}
  \item{absValueCut}{
A \code{numeric} giving the cut-off for partial correlation element selection based on absolute value thresholding. 
}
  \item{FDRcut}{
A \code{numeric} giving the cut-off for partial correlation element selection based on local false discovery rate (FDR) thresholding.
}
  \item{verbose}{
A \code{logical} indicating if intermediate output should be printed on screen.
}
}
\details{
The function transforms the possibly regularized input precision matrix to a partial correlation matrix. Subsequently,
the support of this partial correlation matrix is determined. Support determination is performed either by simple
thresholding on the absolute values of matrix entries (\code{threshold = "absValue"}) or by usage of local FDR (\code{threshold = "localFDR"}). 
The argument \code{absValueCut} is only used when \code{threshold = "absValue"}. The argument \code{FDRcut} is only used when
\code{threshold = "localFDR"}. The function is to some extent a wrapper around 
certain \href{http://cran.r-project.org/web/packages/GeneNet/index.html}{GeneNet}
and \href{http://cran.r-project.org/web/packages/fdrtool/index.html}{fdrtool} functions when
\code{threshold = "localFDR"}. In that case a mixture model is fitted to the nonredundant partial correlations by 
\href{http://cran.r-project.org/web/packages/fdrtool/index.html}{fdrtool}. The decision to
retain elements is then based on the argument \code{FDRcut}. Elements with a posterior probability \eqn{\geq} FDRcut (equalling
1 - local FDR) are retained. See Schaefer and Strimmer (2005) for further details on usage of local FDR in graphical modeling.
}
\value{
If the input \code{P} is a standardized precision matrix the function returns a sparsified partial correlation \code{matrix}.
If the input \code{P} is an unstandardized precision matrix the function returns an object of class \code{list}:
\item{sparseParCor}{A \code{matrix} representing the sparsified partial correlation/standardized precision matrix}
\item{sparsePrecision}{A \code{matrix} representing the sparsified precision matrix}
}
\references{
Schaefer, J., and Strimmer, K. (2005). A shrinkage approach to large-scale covariance estimation
and implications for functional genomics. Statistical Applications in Genetics and Molecular Biology, 4:32.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{ridgeS}}, \code{\link{optPenalty.aLOOCV}}, \code{\link{optPenalty.LOOCV}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]

## Obtain regularized precision under optimal penalty
OPT <- optPenalty.LOOCV(X, lambdaMin = .5, lambdaMax = 30, step = 100)

## Determine support regularized (standardized) precision under optimal penalty
sparsify(OPT$optPrec, threshold = "localFDR")
}