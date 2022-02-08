llqr <- function(x,y,tau,t,h){

  if(!is.numeric(tau)|!is.vector(tau)|any(!is.finite(tau))) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(length(tau)!=1) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(!is.numeric(tau)|(tau>=1)|(tau<=0)) stop("The parameter 'tau' must be a single number between 0 and 1")
  
  if(!is.numeric(tau)|!is.vector(tau)|any(!is.finite(tau))) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(length(tau)!=1) stop("The parameter 'tau' must be a single number between 0 and 1")
  if(!is.numeric(tau)|(tau>=1)|(tau<=0)) stop("The parameter 'tau' must be a single number between 0 and 1")

  if (!is.numeric(x)){stop("'x' must be numeric")}
  if (!is.numeric(y)){stop("'y' must be numeric")}
  if(!length(x)>1 | !length(y)>1){stop("'x' and 'y' must be numeric vectors with the same length")}
  if(length(x)!=length(y)){stop("'x' and 'y' must be numeric vectors with the same length")}
  
  iaux=complete.cases(x,y)
  if(sum(!iaux)!=0){ii=which(iaux==FALSE);x=x[-ii];y=y[-ii];warning("Missing values have been removed from 'x' and 'y'")}

  nt=length(t);n=length(x); beta0=1;beta1=0;z=c();weight=c()
  storage.mode(n)="integer"; storage.mode(z)="double"; storage.mode(y)="double"
  storage.mode(weight)="double"; storage.mode(beta0)="double"; storage.mode(beta1)="double"
  y.estimated<-numeric(nt)
  for (i in 1:nt){
    z<-x-t[i]
    aux<-dnorm(z/h);weight=aux/sum(aux)
    if(max(weight)>=0.9999){
      y.estimated[i]<-y[which.max(weight)]
    }else{
      result<-.Fortran("barro", n, z, y, weight, beta0, beta1,tau)
      y.estimated[i]<-result[[5]]
    }
  }
  return(list("x.values"=t,"y.values"=y.estimated))
}
