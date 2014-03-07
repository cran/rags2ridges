\name{rags2ridges-package}
\alias{rags2ridges-package}
\alias{rags2ridges}
\docType{package}
\title{
Ridge estimation for high-dimensional precision matrices
}
\description{
Package contains proper L2-penalized ML estimators for the precision matrix as well as supporting functions to employ these estimators in a graphical modeling setting.
}
\details{
The main function of the package is \code{\link{ridgeS}} which enables archetypal and proper alternative ML ridge estimation of the precision matrix. The alternative 
ridge estimators can be found in van Wieringen and Peeters (2014) and encapsulate both target and non-target shrinkage for the multivariate normal precision 
matrix. The estimators are analytic and enable estimation in large \eqn{p} small \eqn{n} settings. Supporting functions to employ these estimators in a graphical
modeling setting are also given. These supporting functions enable, a.o., the determination of the optimal value of the penalty parameter, the determination of 
the support of a shrunken precision estimate, as well as various visualization options.
}
\author{
Carel F.W. Peeters, Wessel, N. van Wieringen
\cr
Maintainer: Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\references{
van Wieringen, W.N. and Peeters, C.F.W. (2014). Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data. arXiv:1403.0904 [stat.ME].
}