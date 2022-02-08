\name{bwCV}
\alias{bwCV}
\title{Computing the cross-validation bandwidth proposed by Abberger (1998)}
\description{Function to compute a bandwidth for local linear quantile regression following the cross-validation criteria presented by Abberger (1998).}
\usage{bwCV(x, y, hseq, tau)}
\arguments{
  \item{x}{numeric vector of \code{x} data.}
  \item{y}{numeric vector of \code{y} data. This must be the same length as \code{x}.}
  \item{hseq}{sequence of values where the cross-validation function will be evaluated.}
  \item{tau}{the quantile order where the regression function is to be estimated. It must be a number strictly between 0 and 1.}
}
\details{The cross-validation function is evaluated at each element of \code{hseq}. Then, the cross-validation selector will be the element of \code{hseq} that minimizes the cross-validation function.}
\value{Returns a number with the chosen bandwidth.}
\author{Mercedes Conde-Amboage and Cesar Sanchez-Sellero.}
\references{Abberger, K. (1998). Cross-validation in nonparametric quantile regression. Allgemeines Statistisches Archiv, 82, 149-161.

Abberger, K. (2002). Variable data driven bandwidth choice in nonparametric quantile regression. Technical Report.}
\seealso{The obtained bandwidth can be used in the function \code{\link{llqr}} to produce a local linear estimate of the \code{tau}-quantile regression function.}
\examples{
set.seed(1234)
x=runif(100)
y=10*(x^4+x^2-x)+rexp(100)
hseq=seq(0.05,0.8,length=21)
tau=0.25
bwCV(x,y,hseq,tau)
}
\keyword{local linear quantile regression}
\keyword{bandwidth selection}
