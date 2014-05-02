\name{conditionNumber}
\alias{conditionNumber}
\title{
Visualize the spectral condition number against the regularization parameter
}
\description{
Function that visualizes the spectral condition number of the regularized precision matrix against the regularization parameter. 
The function can be used to heuristically determine the (minimal) value of the penalty parameter.
}
\usage{
conditionNumber(S, lambdaMin, lambdaMax, step, type = "Alt", target = diag(1/diag(S)), 
verticle = FALSE, value = 1, verbose = TRUE)
}
\arguments{
  \item{S}{
Sample covariance \code{matrix}.
}
  \item{lambdaMin}{
A \code{numeric} giving the minimum value for the penalty parameter.
}
  \item{lambdaMax}{
A \code{numeric} giving the maximum value for the penalty parameter.
}
  \item{step}{
An \code{integer} determining the number of steps in moving through the grid [\code{lambdaMin}, \code{lambdaMax}].
}
  \item{type}{
A \code{character} indicating the type of ridge estimator to be used. Must be one of: "Alt", "ArchI", "ArchII".
}
  \item{target}{
A target \code{matrix} (in precision terms) for Type I ridge estimators.
}
  \item{verticle}{
A \code{logical} indicating if output graph should come with a verticle line at a pre-specified value for the penalty parameter.
}
  \item{value}{
A \code{numeric} indicating a pre-specified value for the penalty parameter.
}
  \item{verbose}{
A \code{logical} indicating if intermediate output should be printed on screen.
}
}
\details{
The proposed ridge estimators (see \code{\link{ridgeS}}) are rotation equivariant: The effect of the ridge penalty on the precision estimate is equivalent to
shrinkage of the eigenvalues of the unpenalized precision estimate \eqn{\mathbf{S}^{-1}}. Maximum shrinkage implies that all eigenvalues are 
forced to be equal. The spectral condition number (ratio of maximum to minimum eigenvalue) of the regularized precision matrix may
function as a heuristic in determining the (minimal) value of the penalty parameter. A matrix with a high condition number is near-singular
and its inversion is numerically unstable. Such a matrix is said to be ill-conditioned. Numerically, ill-conditioning will mean that
small changes in the penalty parameter lead to dramatic changes in the condition number. From a numerical point of view one can thus track
the range of the penalty parameter for which the regularized precision matrix is ill-conditioned. When plotting the condition number against
the penalty parameter, the point of stabilization of the former will show up as an elbow (when working in the \eqn{p > n} situation). This
suggest the following fast heuristic (in analogy to the use of the scree plot in factor analysis) for determining the (minimal) value
of the penalty parameter: The value of the penalty parameter for which the spectral condition number starts to stabilize
may be termed an acceptable (minimal) value. 

The function outputs a graph of the
spectral condition number over the range [\code{lambdaMin}, \code{lambdaMax}]. When \code{verticle = TRUE} a verticle line is added 
at the constant \code{value}. This option can be used to assess if the optimal penalty obtained by, e.g., the routines \code{\link{optPenaltyCV}}
or \code{\link{optPenalty.aLOOCV}}, has led to a precision estimate that is well-conditioned.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>, Wessel N. van Wieringen
}
\seealso{
\code{\link{covML}}, \code{\link{optPenaltyCV}}, \code{\link{optPenalty.aLOOCV}}, \code{\link{ridgeS}}
}
\examples{
## Obtain some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
Cx <- covML(X)

## Assess spectral condition number across grid of penalty parameter
conditionNumber(Cx, lambdaMin = .001, lambdaMax = 25, step = 100)
}