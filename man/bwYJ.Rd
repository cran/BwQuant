\name{bwYJ}
\alias{bwYJ}
\title{Computing the plug-in bandwidth proposed by Yu and Jones (1998)}
\description{Function to compute a bandwidth selector for local linear quantile regression following the plug-in rule proposed by Yu and Jones (1998).}
\usage{bwYJ(x, y, tau)}
\arguments{
  \item{x}{numeric vector of \code{x} data.}
  \item{y}{numeric vector of \code{y} data. This must be the same length as \code{x}.}
  \item{tau}{the quantile order where the regression function is to be estimated. It must be a number strictly between 0 and 1.}
}
\value{Returns a bandwidth for a local linear estimate of the \code{tau}-quantile regression function.}
\author{Mercedes Conde-Amboage and Cesar Sanchez-Sellero.}
\references{
Ruppert, D., Sheather, S. J. and Wand, M. P. (1995). An efective bandwidth selector for local least squares regression. Journal of the American Statistical Association. 90, 1257-1270.

Yu, K. and Jones, M. C. (1998). Local linear quantile regression. Journal of the American Statistical Association, 93, 228-237.
}
\seealso{The obtained bandwidth can be used in the function \code{\link{llqr}} to produce a local linear estimate of the \code{tau}-quantile regression function.}
\examples{
set.seed(1234)
x=runif(100)
y=10*(x^4+x^2-x)+rexp(100)
tau=0.25
bwYJ(x,y,tau)
}
\keyword{local linear quantile regression }
\keyword{bandwidth selection}
