bwB <- function (x,tau) {

  if(!is.numeric(tau)|!is.vector(tau)|any(!is.finite(tau))) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(length(tau)!=1) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(!is.numeric(tau)|(tau>=1)|(tau<=0)) stop("The parameter 'tau' must be a single number between 0 and 1")
  
  if (!is.numeric(x)){stop("'x' must be numeric")}
  if(!length(x)>1 ){stop("'x' must be a numeric vector")}
 
  iaux=complete.cases(x)
  if(sum(!iaux)!=0){ii=which(iaux==FALSE);x=x[-ii];warning("Missing values have been removed from 'x'")}
 
	x=sort(x);n=length(x)
	m=(0.25*n^(8/9))%/%1;n2=(n*tau)%/%1

	h0=n^(-0.2) * ((4.5 * dnorm(qnorm(tau))^4)/(2 * qnorm(tau)^2 + 1)^2)^0.2
	j1=ceiling(min(tau+h0,1)*n);j2=max(ceiling(max(tau-h0,0)*n),1)

	s=(x[j1]-x[j2])/(2*h0)
	Z=0.5*(n/m)^3*(x[min(n2+2*m,n)]-2*x[min(n2+m,n)]+2*x[max(n2-m,1)]-x[max(n2-2*m,1)])
	return(((4.5*(s/Z)^2)/n)^(1/5))
}
