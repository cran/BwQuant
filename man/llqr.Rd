\name{llqr}
\alias{llqr}
\title{Fitting a local linear quantile regression model}
\description{Function that estimates the quantile regression function using a local linear kernel smoother.}
\usage{llqr(x, y, tau, t, h)}
\arguments{
  \item{x}{numeric vector of \code{x} data.}
  \item{y}{numeric vector of \code{y} data. This must be the same length as \code{x}.}
  \item{tau}{the quantile order where the regression function is to be estimated. It must be a number strictly between 0 and 1.}
  \item{t}{the values of \code{x} at which the quantile regression model is to be estimated.}
    \item{h}{the bandwidth parameter.}
}
\value{
A list with the following components:
  \item{x.values}{the given points at which the evaluation occurs.}
  \item{y.values}{the estimated values of the quantile regression function at the given \code{x.values}.}
}
\author{Mercedes Conde-Amboage and Cesar Sanchez-Sellero.}
\references{
Fan, J., Hu, T. C. and Truong, Y. K. (1994). Robust nonparametric function estimation. Scandinavian Journal of Statistics, 21, 433-446.

Yu, K. and Jones, M. C. (1998). Local linear quantile regression. Journal of the American Statistical Association, 93, 228-237.

}
\seealso{The argument \code{h} with the bandwidth parameter can be fixed to some arbitrary value or chosen by one of the procedures implemented in the functions \code{\link{bwCV}}, \code{\link{bwPI}}, \code{\link{bwRT}} or \code{\link{bwYJ}}.}
\examples{
set.seed(1234)
x=runif(100)
y=10*(x^4+x^2-x)+rexp(100)
tau=0.25
h=bwPI(x,y,tau)
t=seq(0,1,length=101)
m=llqr(x,y,tau,t,h)
plot(x,y)
lines(m$x.values,m$y.values)
}
\keyword{local linear quantile regression }
\keyword{bandwidth selection}
