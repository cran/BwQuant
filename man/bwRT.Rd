\name{bwRT}
\alias{bwRT}
\title{Computing a bandwidth using a rule of thumb}
\description{Function to compute a bandwidth selector for local linear quantile regression following the rule of thumb presented in Section 2.1 of Conde-Amboage and Sanchez-Sellero (2018).}
\usage{bwRT(x, y, tau)}
\arguments{
  \item{x}{numeric vector of \code{x} data.}
  \item{y}{numeric vector of \code{y} data. This must be the same length as \code{x}.}
  \item{tau}{the quantile order where the regression function is to be estimated. It must be a number strictly between 0 and 1.}
}
\value{Returns a bandwidth for a local linear estimate of the \code{tau}-quantile regression function.}
\author{Mercedes Conde-Amboage and Cesar Sanchez-Sellero.}
\references{Conde-Amboage, M. and Sanchez-Sellero, C. (2018). A plug-in bandwidth selector for nonparametric quantile regression. TEST, 28, 423-450. <doi:10.1007/s11749-018-0582-6>.}
\seealso{The obtained bandwidth can be used in the function \code{\link{llqr}} to produce a local linear estimate of the \code{tau}-quantile regression function.}
\examples{
set.seed(1234)
x=runif(100)
y=10*(x^4+x^2-x)+rexp(100)
tau=0.25
bwRT(x,y,tau)
}
\keyword{local linear quantile regression}
\keyword{bandwidth selection}
